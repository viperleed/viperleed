"""Module badpixels of viperleed.gui.measure.camera.

Defines the BadPixelsFinder class that acquires multiple videos
from a camera and computes a pixel-badness array.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-10-09'
__license__ = 'GPLv3+'

from ast import literal_eval
from collections.abc import MutableSequence
from contextlib import nullcontext
from datetime import datetime
import functools
from pathlib import Path

import numpy as np
from scipy.signal import convolve2d
from scipy import ndimage

from PyQt5 import QtCore as qtc

from viperleed import NUMPY2_OR_LATER
from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.camera import tifffile
from viperleed.gui.measure.camera import cameracalibration as _calib
from viperleed.gui.measure.classes.calibrationtask import (
    CalibrationTaskOperation,
    )


N_DARK = 100    # Number of dark frames used for finding bad pixels
N_NOISE = 3000  # Number of dark frames used for finding noisy pixels
_INVOKE = qtc.QMetaObject.invokeMethod
_UNIQUE = qtc.Qt.UniqueConnection
_unique_connect = functools.partial(base.safe_connect, type=_UNIQUE)


def _report_progress(func):
    """Decorate calculation functions."""

    @functools.wraps(func)
    def _wrapper(self, *args, **kwargs):
        """Wrap decorated function."""
        # Names are of the form "find_WHAT_other_stuff"
        _, what, *_ = func.__name__.split('_')
        section = _FinderSection.from_short_name(what)
        steps = section.n_steps
        self.progress_occurred.emit(self.progress_name, *section, 0, steps)
        ret_val = func(self, *args, **kwargs)
        self.progress_occurred.emit(self.progress_name, *section, steps, steps)
        return ret_val

    return _wrapper


class BadPixelsFinderErrors(base.ViPErLEEDErrorEnum):
    """Class for bad-pixel-finder errors."""

    FLAT_FRAME_WRONG_LIGHT = (
        211,
        "Flat frame is too {}. Cannot automatically "
        "adjust exposure time and gain.")


class _FinderSection(CalibrationTaskOperation):
    """Enumeration class for bad-pixel finder sections."""

    ACQUIRE_DARK_SHORT, ACQUIRE_DARK_LONG, ACQUIRE_DARK_MEDIUM = 1, 2, 3
    ACQUIRE_FLAT, CALCULATE_FLICKERY, CALCULATE_HOT = 4, 5, 6
    CALCULATE_DEAD, CALCULATE_BAD, DONE = 7, 8, 9

    @property
    def long_name(self):
        """Return a long description for this section."""
        if self is _FinderSection.DONE:
            return 'Done.'
        operation, what, *_ = self.name.split('_')

        if operation == 'ACQUIRE':
            txt = f'Acquiring {what.lower()} frame'
            txt += 's ' if what == 'DARK' else ' '
            return self._format_long_name(txt + 'with {:.1f} ms exposure...')
        # Calculating
        return f'Calculating coordinates of {what.lower()} pixels...'

    @property
    def n_steps(self):
        """Return the number of steps in this section."""
        if self is _FinderSection.ACQUIRE_DARK_MEDIUM:
            return N_NOISE
        if 'DARK' in self.name:
            return N_DARK
        if self is _FinderSection.CALCULATE_FLICKERY:
            return 2  # short- and long-term telegraph noise
        return 1

    @classmethod
    def from_short_name(cls, name):
        """Return an instance from a short name."""
        if name.lower() in ('dark-short', 'dark-long', 'dark-medium', 'flat'):
            elem = f"ACQUIRE_{name.replace('-', '_').upper()}"
        elif name.lower() == 'done':
            elem = "DONE"
        else:
            elem = f"CALCULATE_{name.upper()}"
        return getattr(cls, elem)


class BadPixelsFinder(_calib.CameraCalibrationTask):
    """Class for finding bad pixels."""

    name = "Bad-pixel detection"
    progress_name = "Finding bad pixels..."
    description = (
        "Detect bad pixels in a camera. Bad pixels are especially "
        "bright, especially dark, or especially noisy pixels. They "
        "are corrected by replacing with a suitable average of "
        "neighbours. This makes your LEED-IV movies significantly "
        "cleaner. It should be performed regularly (e.g., every 6 "
        "months), as camera sensors may degrade over time."
        )

    def __init__(self, camera):
        """Initialize finder task."""
        super().__init__(camera)

        # __adjustments is used for the exposure/gain of flat;
        # __force_exposure is a flag keeping track of whether we
        # are resorting to forcing the exposure/gain to averages.
        # See __adjust_exposure_and_gain for more info. When we
        # force the exposure, we will use a looser constraint on
        # whether __frame_acceptable.
        self.__adjustments = _Adjustments()
        self.__force_exposure = False

        # Adjust camera settings: full sensor, no binning, no averaging
        _, (max_roi_w, max_roi_h), *_ = camera.get_roi_size_limits()
        full_roi = f"(0, 0, {max_roi_w}, {max_roi_h})"
        self._task_settings.set('camera_settings', 'roi', full_roi)
        self._task_settings.set('camera_settings', 'binning', '1')
        self._task_settings.set('measurement_settings', 'n_frames', '1')
        self._task_settings.set('camera_settings', 'bad_pixels_path', '')

        # Load settings in the camera. This is needed such
        # that self.frame_info returns the correct values
        # NB: THIS MAY NOT BE TRUE IF CAMERA IS IN ANOTHER THREAD.
        # In that case, we would have to ensure there's a way to
        # distinguish between using a queued and blocking-queued
        # connection types for the underlying _INVOKE.
        self.update_device_settings()

        # _imgs contains: the sum and sum of squares of short- and
        # long-exposure movies for the 'dark' frame (i.e., camera
        # with a cap on) and one frame for a 'flat field' (e.g.,
        # white paper right in front of the lens).
        height, width, dtype = self.frame_info
        self._imgs = {
            _FinderSection.ACQUIRE_DARK_SHORT:
                BadPixelsSumStorage(height, width),
            _FinderSection.ACQUIRE_DARK_LONG:
                BadPixelsSumStorage(height, width),
            _FinderSection.ACQUIRE_DARK_MEDIUM:
                BadPixelsMaxMinStorage(height, width, dtype),
            _FinderSection.ACQUIRE_FLAT:
                BadPixelsSumStorage(height, width),
            }
        self._current_section = _FinderSection.ACQUIRE_DARK_SHORT

        # _badness contains a measure of the badness of
        # each pixel. 'Normal' pixels have badness 0, bad
        # pixels have badness >= 6, very good pixels < 0.
        self._badness = np.zeros((height, width), dtype=float)

        # __bad_pixels will contain info about bad
        # pixel coordinates and their replacements
        self.__bad_pixels = BadPixels(self.camera)

    @property
    def large_gain(self):
        """Return a large gain in decibel."""
        _, max_gain = self._limits['gain']
        return min(20, max_gain)

    @property
    def long_exposure(self):
        """Return a long exposure time in milliseconds."""
        _, max_exposure = self._limits['exposure']
        return min(1000, max_exposure)

    @property
    def medium_exposure(self):
        """Return a medium exposure time in milliseconds."""
        _, exposure, _ = sorted((*self._limits['exposure'], 200))
        return exposure

    @property
    def missing_frames(self):
        """Return whether frames are missing for the current section."""
        if self._current_section is _FinderSection.ACQUIRE_DARK_MEDIUM:
            todo = N_NOISE
        elif 'DARK' in self._current_section.name:
            todo = N_DARK
        else:
            todo = 1
        return self._frames_done < todo

    @property
    def short_exposure(self):
        """Return a short exposure time in milliseconds."""
        min_exposure, _ = self._limits['exposure']
        return max(20, min_exposure)

    @qtc.pyqtSlot()
    def abort(self):
        """Abort self and any running tasks."""
        if self.is_aborted:
            return
        for task in self.camera.calibration_tasks['bad_pixels']:
            if task.is_running:
                _INVOKE(task, 'abort')
        super().abort()

    def begin_acquiring(self):
        """Start acquiring images."""
        if self._current_section is _FinderSection.ACQUIRE_DARK_SHORT:
            self.set_info_text(
                "A 'dark' movie will be acquired.<br>Make sure no light "
                "enters the camera.<br>For example, you can <b>close the "
                "lens with its cap.</b>"
                )
        elif self._current_section is _FinderSection.ACQUIRE_FLAT:
            self.set_info_text(
                "A 'flat' frame will be acquired.<br>Make sure the camera "
                "is viewing a featureless background with uniform "
                "illumination.<br>For example, you can place a <b>white "
                "piece of paper right in front of the lens</b>. Normal room "
                "light should be uniform enough."
                )
        if self._current_section not in (_FinderSection.ACQUIRE_DARK_LONG,
                                         _FinderSection.ACQUIRE_DARK_MEDIUM):
            # Show dialog. Its .accepted signal will trigger
            # .continue_(), its .rejected will trigger .abort()
            self.show_info()
            return
        self.continue_()

    @qtc.pyqtSlot(np.ndarray)
    def _check_and_store_frame(self, frame):
        """Store a new frame if it is acceptable."""
        sec = self._current_section
        if not self.__frame_acceptable(frame):
            if 'DARK' in sec.name:
                self.emit_error(
                    _calib.CameraCalibrationErrors.DARK_FRAME_TOO_BRIGHT
                    )
            else:
                self.__adjust_exposure_and_gain(frame)
            return

        self._imgs[sec].add_frame(frame)
        self._frames_done += 1
        self.__report_acquisition_progress()
        self._trigger_next_frame()

    @qtc.pyqtSlot()
    def continue_(self, *_):
        """Continue acquiring after user confirmation."""
        # In all cases start from the same gain and pick exposure
        gain = self.large_gain
        if self._current_section is _FinderSection.ACQUIRE_DARK_SHORT:
            exposure = self.short_exposure
        elif self._current_section is _FinderSection.ACQUIRE_DARK_LONG:
            exposure = self.long_exposure
        elif self._current_section is _FinderSection.ACQUIRE_DARK_MEDIUM:
            exposure = self.medium_exposure
        else:
            # Here we use the large gain and a somewhat short
            # exposure to begin with, but these are adjusted when
            # frames come back to have an average image intensity
            exposure = self.short_exposure
        self._task_settings.set("measurement_settings", "gain",
                                str(gain))
        self._task_settings.set("measurement_settings", "exposure",
                                str(exposure))

        # Apply the new settings and start the camera. When the camera
        # is done preparing it will emit .started(), which will trigger
        # acquisition of the first frame. Each time a frame is done
        # processing (image_processed) it is checked and stored. Then a
        # new acquisition is triggered, if needed.
        self.update_device_settings()
        self.__report_acquisition_progress(exposure)
        self.start_camera()

    @qtc.pyqtSlot()
    def find(self):
        """Find bad pixels in the camera."""
        self.start()

    @_report_progress
    def find_bad_and_replacements(self):  # too-many-locals
        """Find bad pixels and their optimal replacements.

        The best replacements are chosen among opposing
        neighbors such that their badness is minimal.

        Returns
        -------
        None.
        """
        # Each bad pixel will be replaced by the average of
        # two best neighbors, chosen among the opposing ones.
        # The badness of neighbors is multiplied by weights
        # that scale with the square of the distance:
        #
        #                    5 4 5
        #                  5 2 1 2 5
        #                  4 1 x 1 4
        #                  5 2 1 2 5
        #                    5 4 5

        offsets = np.array(((1, 0), (0, 1), (1, 1), (1, -1), (2, 0),
                            (0, 2), (1, 2), (2, 1), (1, -2), (2, -1)))
        badness = self._badness.copy()

        bad_coords = np.asarray(np.where(badness >= 6))
        bad_y, bad_x = bad_coords
        badness[bad_y, bad_x] = np.inf
        width, height = badness.shape

        # For each offset, calculate the total badness
        # of bad+offset and bad-offset, 'skipping'
        # those offsets that wrap around the image.
        total_badness = np.ones((len(offsets), *badness.shape)) * np.inf
        for i, offset in enumerate(offsets):
            repl_plus = bad_coords.T + offset
            repl_minus = bad_coords.T - offset
            weighted_badness = badness * (offset[0]**2 + offset[1]**2)

            # Pick only those bad pixels whose replacements
            # are on the inside of the image.
            inside = np.all((np.all(repl_plus >= 0, axis=1),
                             np.all(repl_minus >= 0, axis=1),
                             np.all(repl_plus < (width, height), axis=1),
                             np.all(repl_minus < (width, height), axis=1)),
                            axis=0)

            plus_y, plus_x = repl_plus[inside].T
            minus_y, minus_x = repl_minus[inside].T
            bad_inside_y, bad_inside_x = np.transpose(bad_coords.T[inside])

            total_badness[i][bad_inside_y, bad_inside_x] = (
                weighted_badness[minus_y, minus_x]
                + weighted_badness[plus_y, plus_x]
                )

        # Finally, pick as replacements for the bad pixels those
        # offsets that realize the minimum overall badness
        best_badness = total_badness.min(axis=0)

        # Eliminate from the bad pixels (and replacements) those
        # that are uncorrectable. They are such that the 'best'
        # replacement is still a bad pixel. Keep also track of them,
        # for information purposes.
        correctable = np.isfinite(best_badness[bad_y, bad_x])
        uncorrectable = np.asarray((bad_y[~correctable],
                                    bad_x[~correctable])).T

        bad_y, bad_x = bad_y[correctable], bad_x[correctable]

        # Finally, get the offsets of the replacement pixels,
        # and store all the info in a BadPixels object.
        best_offset_indices = total_badness.argmin(axis=0)[bad_y, bad_x]
        self.__bad_pixels = BadPixels(
            self.camera,
            bad_coordinates=np.asarray((bad_y, bad_x)).T,
            replacement_offsets=offsets[best_offset_indices],
            uncorrectable=uncorrectable
            )

    @_report_progress
    def find_dead_pixels(self):
        """Detect dead pixels from flat frame.

        Dead pixels are those whose intensity is much smaller
        than the average of the neighbours (excluding flickery
        and hot pixels). This method adds to the badness the
        ratio I_neighbors/I_pixel - 1.

        Return
        ------
        None.
        """
        pix_min, _, intensity_range = self._limits['intensity']
        # We can use .mean() because there is only a single flat frame.
        flat = self._imgs[_FinderSection.ACQUIRE_FLAT].mean() - pix_min

        # We have to exclude the already-detected bad pixels.
        # For this to happen, we will (1) add a bit to all
        # pixels (so that no pixel is zero), (2) set the bad
        # ones to zero, (3) convolve with an appropriate
        # averaging kernel, also counting nonzero pixels (i.e.,
        # the good ones). This way, bad ones do not factor in
        # the neighbour average.

        offset = 0.1*intensity_range
        flat_bad = flat + offset
        flat_bad_ones = np.ones(flat_bad.shape)  # for counting

        bad_px_mask = self._badness >= 6

        flat_bad[bad_px_mask] = 0
        flat_bad_ones[bad_px_mask] = 0

        # 5x5 kernel for averaging a diamond-shaped patch of neighbours
        kernel = np.array(((0, 0, 1, 0, 0),
                           (0, 1, 1, 1, 0),
                           (1, 1, 0, 1, 1),
                           (0, 1, 1, 1, 0),
                           (0, 0, 1, 0, 0)))

        neighbor_sum = convolve2d(flat_bad, kernel, mode='same')
        neighbor_count = convolve2d(flat_bad_ones, kernel, mode='same')
        neighbor_ave = neighbor_sum / neighbor_count - offset

        # Isolated pixels in a patch of bad neighbours
        # will also be marked as bad.
        neighbor_ave[neighbor_count < 0.5] = np.inf

        # As will those that have both zero neighbour average
        # and zero intensity (dead in a cluster of bad ones)
        delta_badness = neighbor_ave / flat - 1
        delta_badness[np.isnan(delta_badness)] = np.inf

        self._badness += delta_badness

    def find_flickery_pixels(self):
        """Prepare pixel badness based on short- and long-term noise."""
        _section = _FinderSection.CALCULATE_FLICKERY
        self.progress_occurred.emit(self.progress_name, *_section,
                                    0, _section.n_steps)
        self._find_fast_flickery_pixels()
        self.progress_occurred.emit(self.progress_name, *_section,
                                    1, _section.n_steps)
        self._find_long_term_flickery_pixels()
        self.progress_occurred.emit(self.progress_name, *_section,
                                    _section.n_steps, _section.n_steps)

    @_report_progress
    def find_hot_pixels(self):
        """Set badness of hot pixels to infinity."""
        pix_min, _, intensity_range = self._limits['intensity']
        _long = self._imgs[_FinderSection.ACQUIRE_DARK_LONG]
        dark_ave = _long.mean() - pix_min
        hot_pixels = dark_ave >= 0.25*intensity_range

        self._badness[hot_pixels] = np.inf

    @qtc.pyqtSlot(tuple)
    def _on_device_error(self, error):
        """Silence all bad-pixels-related errors."""
        # The difference with the base-class implementation is that
        # self IS NOT part of the list of calibration tasks, but
        # we should nonetheless skip all bad-pixels-related errors
        if self.camera.is_bad_pixels_error(error):
            return True
        return super()._on_device_error(error)

    def restore_device(self):
        """Restore settings and other device attributes."""
        # Load either old or newly saved bad pixels file
        _INVOKE(self.camera, 'update_bad_pixels')
        for task in self.camera.calibration_tasks['bad_pixels']:
            task.mark_as_done(False)
        super().restore_device()

    def save_and_cleanup(self):
        """Save a bad pixels file and finish."""
        bp_path = self.original_settings.get("camera_settings",
                                             "bad_pixels_path",
                                             fallback='')
        self.__bad_pixels.write(bp_path)
        self.restore_device()
        done = _FinderSection.DONE
        self.progress_occurred.emit(self.progress_name, *done,
                                    done.n_steps, done.n_steps)
        self.mark_as_done(True)

    @qtc.pyqtSlot()
    def start(self):
        """Start this task."""
        if not super().start():
            return False

        # See if there is any preliminary task to be performed
        task = self.__get_preliminary_task()
        if task:
            # Make sure all tasks have the same original settings
            task.original_settings.read_dict(self.original_settings)
            _unique_connect(task.done, self.__on_preliminary_task_done)
            _unique_connect(task.aborted, self.abort)
            _unique_connect(task.error_occurred, self.error_occurred)
            _unique_connect(task.progress_occurred, self.progress_occurred)
                            # self.__on_preliminary_task_progress)              # TODO: would be nice to have only one "overall" progress bar
            _INVOKE(task, 'start')
            return False

        # No more tasks to do. Proceed to finding bad pixels.
        # Clear any bad pixel info in the camera (frames should be raw)
        bad_old = self.camera.bad_pixels
        if bad_old:
            bad_old.clear()
        self._current_section = _FinderSection.first()
        self.begin_acquiring()
        return True

    def _find_fast_flickery_pixels(self):
        """Prepare badness based on how flickery pixels are.

        Badness will be set according to how much each pixel
        fluctuates in the 'dark' frames. We take as a measure
        of badness the ratio between the variance of a pixel
        and the overall average variance: the intensity of
        flickery pixels varies much more than the average.

        Two contributions are taken into account: long-exposed
        dark frames are used to detect pixels that flicker due
        to higher photon sensitivity, whereas short-exposed
        frames detect pixels that are flickery due to readout
        noise.

        Returns
        -------
        None.
        """
        long_flicker = self._imgs[_FinderSection.ACQUIRE_DARK_LONG].var()
        short_flicker = self._imgs[_FinderSection.ACQUIRE_DARK_SHORT].var()
        long_mean, short_mean = long_flicker.mean(), short_flicker.mean()
        self._badness = np.zeros_like(long_flicker)

        # The conditional checks prevent division by zero for very good
        # cameras. '-1' --> badness == 0 for normally flickery pixels
        if short_mean:
            self._badness += short_flicker / short_mean - 1
        if long_mean:
            self._badness += long_flicker / long_mean - 1

    def _find_long_term_flickery_pixels(self):
        """Calculate badness for pixels with sporadic burst noise.

        Badness will be set according to how much each pixel
        fluctuates in the 'dark' frames. As a measure of badness we
        take the difference between maximum and minimum intensity
        picked up by each pixel. The threshold when a pixel is
        considered to have too much telegraph noise is when the
        flicker of the pixel is about 3 times the average flicker.
        (3-1)^2 * 1.5 = 6,
        where 6 is the threshold used for the total badness.

        To determine the long-term flickering we use a medium exposure
        time with a very high frame count. The overall measurement for
        burst noise detection alone should not take less than ten
        minutes.
        """
        flicker = self._imgs[_FinderSection.ACQUIRE_DARK_MEDIUM].range_
        flicker_mean = flicker.mean()
        self._badness += ((flicker / flicker_mean - 1)**2) * 1.5

    @qtc.pyqtSlot()
    def _trigger_next_frame(self, *_):
        """Trigger acquisition of a new frame if necessary."""
        if not self.can_trigger_camera:
            return

        if self.missing_frames:
            self.trigger_camera()
            return

        # One section over, go to the next one
        self._current_section = self._current_section.next_()
        if self._current_section is _FinderSection.CALCULATE_FLICKERY:
            # Done with all sections. Can proceed to calculations.
            self.find_flickery_pixels()
            self.find_hot_pixels()
            self.find_dead_pixels()
            self.find_bad_and_replacements()
            self.save_and_cleanup()
            return
        # Proceed with the next acquisition
        self._frames_done = 0
        self.__adjustments.clear()
        self.begin_acquiring()

    def __adjust_exposure_and_gain(self, frame):
        """Use frame to estimate optimal exposure/gain for flat frames."""
        old_exposure = self.camera.exposure
        old_gain = self.camera.gain
        pix_min, _, intensity_range = self._limits['intensity']

        intensity = frame.mean() - pix_min
        exposure = old_exposure * 0.5 * intensity_range / intensity
        exposure, gain = self.__correct_exposure_gain(exposure, old_gain)
        if exposure is None:
            return
        self.__adjustments.add(exposure * 10**(gain / 20))

        # Now decide whether we really want to set these new values. We
        # do not want this process to take forever, in case there are
        # oscillations. Use the estimates only in case we have already
        # done a lot of attempts,  and only if we can provide a reasonable
        # estimate for the new exposure/gain.
        self.__force_exposure = False
        estimate = self.__adjustments.stable_value
        if len(self.__adjustments) >= 20 and estimate is not None:              # TODO: should we also have a hard limit with a 'retry' button?
            exposure = estimate * 10**(-gain / 20)
            exposure, gain = self.__correct_exposure_gain(exposure, gain)
            self.__force_exposure = True

        self._task_settings.set('measurement_settings', 'exposure',
                                f'{exposure:.3f}')
        self._task_settings.set('measurement_settings', 'gain', f'{gain:.1f}')
        self.update_device_settings()
        self.__report_acquisition_progress(exposure)
        self.start_camera()

    def __correct_exposure_gain(self, exposure, gain):
        """Return in-range exposure and gain."""
        exposure_min, exposure_max = self._limits['exposure']
        if exposure < exposure_min:
            # Decrease the gain
            gain += 20*np.log10(exposure/exposure_min)
            exposure = exposure_min
        elif exposure > exposure_max:
            # Increase the gain (this is very unlikely to happen)
            gain += 20*np.log10(exposure/exposure_max)
            exposure = exposure_max

        # Now check that the gain is OK: the only cases in which
        # it cannot be OK is if we need very short or very long
        # exposures (otherwise the gain stayed the same).
        min_gain, max_gain = self._limits['gain']
        if gain < min_gain:
            # Too much intensity
            self.emit_error(BadPixelsFinderErrors.FLAT_FRAME_WRONG_LIGHT,
                            'bright')
            return None, 0
        if gain > max_gain:
            # Too little intensity
            self.emit_error(BadPixelsFinderErrors.FLAT_FRAME_WRONG_LIGHT,
                            'dark')
            return None, 0
        return exposure, gain

    def __frame_acceptable(self, frame):
        """Return whether frame has acceptable intensity.

        Whether a frame is acceptable depends on the current
        'section': if we're acquiring a dark frame, the average
        intensity should be low. For a flat one, instead, we
        would like to end up with an average intensity in the
        middle of the pixel value range.

        Parameters
        ----------
        frame : numpy.ndarray
            The frame to be checked

        Returns
        -------
        acceptable : bool
            True if frame is OK.
        """
        pixel_min, _, intensity_range = self._limits['intensity']
        relative_intensity = (frame.mean() - pixel_min) / intensity_range
        name = self._current_section.name
        if 'DARK' in name:
            return relative_intensity < 0.2
        _min, _max = (0.45, 0.55) if not self.__force_exposure else (0.3, 0.7)
        return _min < relative_intensity < _max                                 # TODO: do we want to complain if there are very bright areas (gradients)?

    def __get_preliminary_task(self):
        """Return a CameraCalibrationTask to run before self, or None."""
        tasks_to_do = (t for t in self.camera.calibration_tasks['bad_pixels']
                       if t.to_be_done)
        return next(tasks_to_do, None)

    @qtc.pyqtSlot()
    def __on_preliminary_task_done(self):
        """React to a preliminary task being finished."""
        task = self.sender()

        # Make sure to update the original settings to the ones
        # that may have been edited by the task itself
        self.original_settings.read_dict(task.original_settings)
        self.start()

    def __report_acquisition_progress(self, exposure=None):
        """Report progress of acquisition."""
        if exposure is None:
            exposure = self.camera.exposure
        section = self._current_section
        section.set_format(exposure)
        n_adjustments = len(self.__adjustments)
        self.progress_occurred.emit(self.progress_name, *section,
                                    self._frames_done + n_adjustments,
                                    section.n_steps + n_adjustments)


class _Adjustments(MutableSequence):
    """Class for handling Holt-Winters smoothing of camera exposure/gain.

    This is used to assess whether a series of exposures/gains has
    settled, as well as determining the average value in the settled
    region. Holt-Winters exponential smoothing is very good to detect
    'ends of transients'
    """

    def __init__(self, *args, **kwargs):
        """Initialize instance."""
        self.__params = {'alpha': 0, 'beta': 0,
                         'thresholds': kwargs.get('thresholds', (7/100, 0.1))}
        self.alpha = kwargs.get('alpha', 0.6)
        self.beta = kwargs.get('alpha', 0.2)

        self.__list = []
        self.__holt_winters = [], []
        self.__last_settled = None
        for value in args:
            self.add(value)

    @property
    def alpha(self):
        """Return the smoothing factor for levels."""
        return self.__params['alpha']

    @alpha.setter
    def alpha(self, new_alpha):
        """Set a value for the smoothing factor for levels."""
        self.__params['alpha'] = new_alpha
        self.__params['alphac'] = 1 - new_alpha

    @property
    def beta(self):
        """Return the smoothing factor for trends."""
        return self.__params['beta']

    @beta.setter
    def beta(self, new_beta):
        """Set a value for the smoothing factor for trends."""
        self.__params['beta'] = new_beta
        self.__params['betac'] = 1 - new_beta

    @property
    def is_settled(self):
        """Return whether self has settled."""
        return self.__last_settled is not None

    @property
    def settled_values(self):
        """Return the values in the most recent 'settled' region."""
        return [] if not self.is_settled else self[self.__last_settled:]

    @property
    def stable_value(self):
        """Return the mean stable value, or None if not stable."""
        # Consider stable if the relative standard error of
        # the mean in settled_values is less than 10%
        values = np.array(self.settled_values)
        n_values = len(values)
        if n_values < 2:  # Need 2 values for s.err calculation
            return None
        mean_x, mean_x2 = values.mean(), (values**2).mean()
        rel_serr = mean_x2 / mean_x**2 - 1
        rel_serr = np.sqrt(rel_serr / (n_values - 1))
        if 3 * rel_serr >= self.__stable_thresh:
            return None
        return mean_x

    def __delitem__(self, index):
        """Prevent removal of items at index."""
        # In principle I could allow removal from the end as
        # this does not require recomputing anything, but I'm
        # not going to use it anyway
        raise NotImplementedError(
            f"Cannot remove item at positon {index}. Deletion is forbidden."
            )

    def __getitem__(self, index):
        """Return item at index."""
        return self.__list[index]

    def __len__(self):
        """Return length of self."""
        return len(self.__list)

    def __repr__(self):
        """Return a string representation of self."""
        txt = f"_Adjustments({repr(self.__list)}, "
        txt += f"alpha={self.alpha}, beta={self.beta}, "
        txt += f"thresholds={self.__params['thresholds']})"
        return txt

    def __setitem__(self, index, value):
        """Prevent setting elements other than at the end."""
        _invalid = (isinstance(index, slice)
                    and (index.step != 1 or index.start != len(self)))
        _invalid |= (isinstance(index, int) and index != len(self))
        if _invalid:
            raise NotImplementedError("Can only add elements to the end!")
        if isinstance(index, int):
            value = (value,)
        for val in value:
            self.add(val)

    def __str__(self):
        """Return a string representation of self."""
        return str(self.__list)

    def add(self, value):
        """Add one value, computing median smoothing, level and trend."""
        self.__list.append(value)

        # 5-element-window median can still change the i-th element
        # until the (i+3)-th arrives. Makes no sense to calculate
        # anything till the 4th element is here.
        if len(self) < 4:
            return

        level, trend = self.__holt_winters

        idx = len(self) - 4   # This is the stable value
        self.__median = ndimage.median_filter(self, 5, mode='nearest')
        median = self.__median[idx]
        if not idx:           # First stable median value
            level.append(median)
            return

        old_level = level[-1]
        if not trend:  # Second stable median value
            trend.append(median - old_level)                     # mul: /
        old_trend = trend[-1]
        new_level = (self.alpha * median
                     + self.__alphac * (old_level + old_trend))  # mul: *
        new_trend = (self.beta * (new_level - old_level)         # mul: /
                     + self.__betac * old_trend)
        level.append(new_level)
        trend.append(new_trend)

        # Now we decide if we are settled, and store the
        # index if we are. Notice that min(idx) == 1.
        settled = abs(new_trend) < self.__settle_thresh * median
        if not self.is_settled and settled:
            self.__last_settled = idx
        elif not settled:
            self.__last_settled = None

    def clear(self):
        """Clear list."""
        # No need to clear __median as it is recalculated each time
        for _list in self.__holt_winters:
            _list.clear()
        self.__last_settled = None
        super().clear()

    def insert(self, index, value):
        """Prevent setting elements other than at the end."""
        if index != len(self):
            raise NotImplementedError("Can only add elements to the end!")
        self.add(value)

    @property
    def __alphac(self):
        """Return 1 - self.alpha."""
        return self.__params['alphac']

    @property
    def __betac(self):
        """Return 1 - self.alpha."""
        return self.__params['betac']

    @property
    def __settle_thresh(self):
        """Return the threshold for settling."""
        return self.__params['thresholds'][0]

    @property
    def __stable_thresh(self):
        """Return the threshold for settling."""
        return self.__params['thresholds'][1]


# pylint: disable=too-many-instance-attributes
# Nine seems alright here, especially as they
# all are private attributes.
class BadPixels:
    """Class for reading, writing and manipulating bad pixels."""

    def __init__(self, camera, bad_coordinates=None,
                 replacement_offsets=None, uncorrectable=None):
        """Initialize BadPixels object.

        Parameters
        ----------
        bad_coordinates : Sequence or None, optional
            Coordinates of bad pixels. Should be a (N,2)-shaped
            sequence with N elements, each of which corresponds
            to the row_no and col_no coordinate of a bad pixel.
            If None, methods that manipulate or write bad pixels
            will raise RuntimeError. Bad pixels can be read from a
            file with the .read(filename) method. Default is None.
        replacement_offsets : Sequence or None, optional
            Offsets of pixels to be used as replacements for
            the bad ones. Should be a (N, 2)-shaped sequence
            with the same number of elements N as bad_pixels.
            Each element is a (row_no, col_no) offset (in pixels)
            of pixels to be used for replacement of the bad ones.
            Default is None.
        uncorrectable : Sequence or None, optional
            Coordinates of bad pixels that cannot be corrected
            (because they have only bad pixels as neighbors).

        Returns
        -------
        None.

        Raises
        ------
        TypeError
            If the shape of bad_pixels or replacement_offsets
            is not (N, 2).
        TypeError
            If camera is not a CameraABC subclass.
        ValueError
            If bad_pixels is given without any replacement_offsets
            or the other way around.
        ValueError
            If the number of bad_pixels and replacement_offsets
            is different.
        """
        if ((bad_coordinates is not None and replacement_offsets is None)
            or (replacement_offsets is not None
                and bad_coordinates is None)):
            raise ValueError(
                f"{self.__class__.__name__}: both bad_coordinates and "
                "replacement_offsets should be given when either is."
                )

        if bad_coordinates is None:
            bad_coordinates = np.ones((1, 2))*np.nan
        else:
            bad_coordinates = np.asarray(bad_coordinates)

        if replacement_offsets is None:
            replacement_offsets = np.ones((1, 2))*np.nan
        else:
            replacement_offsets = np.asarray(replacement_offsets)

        if len(bad_coordinates.shape) != 2 or bad_coordinates.shape[1] != 2:
            raise TypeError(f"{self.__class__.__name__}: bad_coordinates "
                            "must be a (N,2)-shaped sequence of coordinates. "
                            f"Found invalid shape {bad_coordinates.shape}.")
        if (len(replacement_offsets.shape) != 2
                or replacement_offsets.shape[1] != 2):
            raise TypeError(
                f"{self.__class__.__name__}: replacement_offsets "
                "must be a (N,2)-shaped sequence of coordinates. "
                f"Found invalid shape {replacement_offsets.shape}."
                )
        if len(bad_coordinates) != len(replacement_offsets):
            raise ValueError(
                f"{self.__class__.__name__}: inconsisten number of "
                f"bad_coordinates ({len(bad_coordinates)}) and of "
                f"replacement_offsets ({len(replacement_offsets)})."
                )

        self.__camera = camera
        self.__bad_coords = bad_coordinates
        self.__replacements = replacement_offsets
        self.__uncorrectable = uncorrectable

        self.__bad_coords_roi = bad_coordinates
        self.__replacements_roi = replacement_offsets
        self.__uncorrectable_roi = uncorrectable
        self.__datetime = datetime.now()

        # Prepare a base name used for saving bad-pixel info to disk
        now = self.__datetime.strftime("%Y%m%d_%H%M%S")
        self.__base_name = f"{self.__camera.name_clean}_{now}"

    def __str__(self):
        """Return a string representation of self."""
        formatter = "{}->{}"
        txt_lst = [formatter.format(b, r)
                   for b, r in zip(self.__bad_coords, self.__replacements)]
        return f"BadPixels({', '.join(txt_lst)})"

    @property
    def bad_pixel_coordinates(self):
        """Return coordinates of bad pixels.

        The coordinates always refer to a region of interest. Only
        correctable bad pixels within the region of interest are
        returned.

        The array returned should NOT be edited in place. Make a
        .copy() if any editing is needed.

        Returns
        -------
        bad_pixel_coordinates : numpy.ndarray
            Shape (N, 2), with N number of bad pixels. Each entry
            corresponds to the (row_no, col_no) coordinates of a
            bad pixel.
        """
        return self.__bad_coords_roi

    @property
    def file_name(self):
        """Return the currently read file name."""
        if not self.__has_info:
            return ''
        return self.__base_name

    @property
    def n_bad_pixels_sensor(self):
        """Return the total number of bad pixels.

        The count includes the uncorrectable ones. Use
        n_bad_pixels_roi for those that are within the
        current region of interest.

        Returns
        -------
        n_bad_pixels_sensor : int
            Total number of bad pixels for this camera.
        """
        return len(self.__bad_coords) + self.n_uncorrectable_sensor

    @property
    def n_bad_pixels_roi(self):
        """Return the number of bad pixels in the current ROI.

        The count includes the uncorrectable ones.

        Returns
        -------
        n_bad_pixels_roi : int
            Number of bad pixels in the current region of interest.
        """
        return len(self.__bad_coords_roi) + self.n_uncorrectable_roi

    @property
    def n_uncorrectable_sensor(self):
        """Return the total number of uncorrectable pixels."""
        if self.__uncorrectable is None:
            return 0
        return len(self.__uncorrectable)

    @property
    def n_uncorrectable_roi(self):
        """Return the number of uncorrectable pixels in the current ROI."""
        if self.__uncorrectable_roi is None:
            return 0
        return len(self.__uncorrectable_roi)

    @property
    def replacement_offsets(self):
        """Return the offsets to pixels that should replace bad ones.

        Bad pixels can be replaced by the average of
        .bad_pixel_coordinates + .replacement_offsets and
        .bad_pixel_coordinates - .replacement_offsets.

        The array returned should NOT be edited in place. Make a
        .copy() if any editing is needed.

        Returns
        -------
        replacement_offsets : numpy.ndarray
            Shape (N, 2), with N number of bad pixels. Entries
            are ordered as in .bad_pixel_coordinates.
        """
        return self.__replacements_roi

    @property
    def uncorrectable(self):
        """Return the coordinates of uncorrectable bad pixels.

        The coordinates always refer to a region of interest.
        Bad pixels are uncorrectable if they only have bad
        neighbours (in the allowed replacement directions).

        The array returned should NOT be edited in place. Make a
        .copy() if any editing is needed.

        Returns
        -------
        uncorrectable : numpy.ndarray
            Shape (M, 2), with M number of uncorrectable bad pixels.
        """
        if self.__uncorrectable_roi is None:
            return np.empty(0)
        return self.__uncorrectable_roi

    def apply_roi(self, no_roi=False):
        """Pick bad pixels and replacements inside the current camera ROI.

        This method must be called before the .bad_pixel_coordinates
        and .replacement_offsets properties return values that are up
        to date with a region of interest.

        Once this method runs, only correctable bad pixels and their
        replacements will be returned by the properties above.

        Parameters
        ----------
        no_roi : bool, optional
            If True, deactivate region of interest completely.
            Default is False.

        Raises
        ------
        RuntimeError
            If called before any bad pixel information was loaded.
        ValueError
            If no_roi=False and some of the other arguments is missing.
        """
        if not self.__has_info:
            raise RuntimeError(f"{self.__class__.__name__}: No bad pixels "
                               "to apply a region of interest to.")

        # Reset all ROI-related info
        self.__bad_coords_roi = self.__bad_coords
        self.__replacements_roi = self.__replacements
        self.__uncorrectable_roi = self.__uncorrectable

        if no_roi:
            return

        top_x, top_y, width, height = self.__camera.roi

        # Recalculate coordinates based on the new origin
        # NB: the bad coordinates are (row#, col#)!
        if self.__bad_coords.size:
            bad_coords = self.__bad_coords - (top_y, top_x)
            repl_minus = bad_coords - self.__replacements
            repl_plus = bad_coords + self.__replacements

            # Mask all bad pixels and replacements outside the ROI:
            mask = np.all((np.all(bad_coords >= (0, 0), axis=1),
                           np.all(repl_minus >= (0, 0), axis=1),
                           np.all(repl_plus >= (0, 0), axis=1),
                           np.all(bad_coords < (height, width), axis=1),
                           np.all(repl_minus < (height, width), axis=1),
                           np.all(repl_plus < (height, width), axis=1)),
                          axis=0)

            self.__bad_coords_roi = bad_coords[mask]
            self.__replacements_roi = self.__replacements[mask]

        # Do the same for uncorrectable pixels (if any)
        if self.__uncorrectable is not None and self.__uncorrectable.size:
            uncorrectable = self.__uncorrectable - (top_y, top_x)
            mask = np.all((np.all(uncorrectable >= (0, 0), axis=1),
                           np.all(uncorrectable < (height, width), axis=1)),
                          axis=0)
            self.__uncorrectable_roi = uncorrectable[mask]

    def clear(self):
        """Clear any bad-pixel information in the object.

        After this method is called, bool(self) returns False.

        Returns
        -------
        None.
        """
        # pylint: disable=unnecessary-dunder-call
        # Looks funky, but it is the fastest way to clear every
        # single attribute of self and remove all bad-pixels info
        self.__init__(self.__camera)

    def get_uncorrectable_clusters_sizes(self):
        """Return a sorted list of sizes for uncorrectable clusters.

        Only the uncorrectable pixels in the current
        region of interest will be considered.

        Returns
        -------
        sizes : list
            Each element is the number of pixels in a
            cluster of uncorrectable ones.

        Raises
        ------
        RuntimeError
            If this method is called before any bad-pixel
            information has been loaded.
        """
        if not self:
            raise RuntimeError("No bad pixel information present")

        if not self.n_uncorrectable_roi:
            return []

        # Prepare an image of uncorrectable pixels
        width, height, *_ = self.__camera.image_info
        uncorrectable = np.zeros((height, width), dtype='uint8')
        bad_y, bad_x = self.uncorrectable.T
        uncorrectable[bad_y, bad_x] = 255

        # Then label 4-connected features, and find their sizes
        labelled, _ = ndimage.label(uncorrectable)
        clusters = ndimage.find_objects(labelled)
        return sorted(
            (np.count_nonzero(labelled[c]) for c in clusters)
            )

    def __get_bad_px_file(self, filepath, most_recent):
        """Return the path to a bad-pixels file for self.camera."""
        filepath = Path(filepath)
        cam_name = self.__camera.name_clean

        # Search in filepath for bad pixels files of the current camera
        bad_px_files = list(filepath.glob(f"{cam_name}*.badpx"))
        if not bad_px_files:
            # Invalidate the information, then complain
            self.__bad_coords = np.ones((1, 2)) * np.nan
            raise FileNotFoundError(
                f"{self.__class__.__name__}: could not find a "
                f"bad pixels file for camera {self.__camera.name} "
                f"in directory {filepath.resolve()}"
                )

        idx = 0      # the most recent
        if not most_recent:
            idx = 1  # the second most recent

        # Get the correct bad pixels file
        try:
            return sorted(((f.stem, f) for f in bad_px_files),
                          reverse=True)[idx]
        except IndexError:
            self.__bad_coords = np.ones((1, 2)) * np.nan
            raise FileNotFoundError(
                f"{self.__class__.__name__}: found only one "
                f"bad pixels file for camera {self.__camera.name} "
                f"in directory {filepath.resolve()}. Cannot load "
                "the second-most recent file."
                ) from None

    def read(self, filepath='', most_recent=True):
        """Read bad pixel information from a file path.

        Parameters
        ----------
        filepath : str or Path, optional
            Path to the folder containing a bad pixel file.
            If an empty string, the current directory will
            be searched. Default is an empty string.
        most_recent : bool, optional
            If True, the most recent bad pixel file for the
            current camera will be read. If False, the second-
            most-recent file will be loaded. If False and there
            was no second-most-recent file, a FileNotFoundError
            is raised. Default is True.

        Raises
        ------
        FileNotFoundError
            If filepath does not contain a bad pixels file
            for the current camera.
        ValueError
            If filepath contains a corrupt bad pixels file
            for the current camera.
        """
        base_name, filename = self.__get_bad_px_file(filepath, most_recent)

        bad_pixels = []
        replacements = []
        uncorrectable = []
        with open(filename, 'r', encoding='utf-8') as bad_px_file:
            for line in bad_px_file:
                if line.strip().startswith('#'):  # Comment
                    continue
                try:
                    bad, replacement = line.split('),')
                except ValueError as err:
                    # Not a valid file
                    raise ValueError("Could not read bad pixel information "
                                     f"from file {filename}") from err
                bad += ')'
                if 'NaN' not in replacement:
                    bad_pixels.append(literal_eval(bad.strip()))
                    replacements.append(literal_eval(replacement.strip()))
                else:
                    # Uncorrectable
                    uncorrectable.append(literal_eval(bad.strip()))

        if not bad_pixels and not uncorrectable:
            # Invalidate the information, then complain
            self.__bad_coords = np.ones((1, 2)) * np.nan
            raise ValueError("Could not read bad pixel information "
                             f"from file {filename}")
        self.__bad_coords = np.array(bad_pixels)
        self.__replacements = np.array(replacements)
        self.__uncorrectable = np.array(uncorrectable)
        self.__base_name = base_name
        self.apply_roi()

    def save_uncorrectable_mask(self, filepath=''):
        """Save an image with uncorrectable pixels as white.

        The image will always contain only the uncorrectable
        pixels within the current ROI. Use .apply_roi(no_roi=True)
        if an image of the uncorrectable pixels on the whole camera
        sensor is needed.

        Parameters
        ----------
        filepath : str or Path, optional
            Base directory where the image should be saved.
            Default is '', i.e., current directory.

        Raises
        ------
        RuntimeError
            If called before any bad pixel info is present.
        """
        self.__save_mask_base(self.uncorrectable, filepath, 'uncorr',
                              "Uncorrectable bad pixels")

    def save_bad_px_mask(self, filepath=''):
        """Save an image with bad pixels as white.

        The image will always contain only the bad pixels within the
        current ROI. Use .apply_roi(no_roi=True) if an image of the
        bad pixels on the whole camera sensor is needed.

        Parameters
        ----------
        filepath : str or Path, optional
            Base directory where the image should be saved.
            Default is '', i.e., current directory.

        Raises
        ------
        RuntimeError
            If called before any bad pixel info is present.
        """
        self.__save_mask_base(self.bad_pixel_coordinates, filepath, 'bad',
                              "Bad pixels")

    def __save_mask_base(self, coordinates, filepath, extra_name, comment):
        """Save coordinates as a binary mask in filepath."""
        if not self.__has_info:
            raise RuntimeError(f"{self.__class__.__name__}: No "
                               "bad pixel information present.")

        # Prepare the image
        width, height, *_ = self.__camera.image_info
        mask = np.zeros((height, width), dtype='uint8')

        if coordinates.size:
            bad_y, bad_x = coordinates.T
            mask[bad_y, bad_x] = 255

        date_time = self.__datetime.strftime("%Y:%m:%d %H:%M:%S")
        filename = Path(filepath, f"{self.__base_name}_{extra_name}_mask.tiff")
        tiff = tifffile.TiffFile.from_array(
            mask,
            comment=f"{comment} mask for {self.__camera.name}",
            model=self.__camera.name,
            date_time=date_time
            )
        tiff.write(filename)

    def write(self, filepath=''):
        """Write bad pixel information to file.

        Parameters
        ----------
        filepath : str or Path, optional
            Path to the folder containing a bad pixel file.
            If an empty string, the file will be saved in the
            current directory. Default is an empty string.

        Raises
        ------
        RuntimeError
            If self has no bad pixel coordinates.
        """
        if not self.__has_info:
            raise RuntimeError(f"{self.__class__.__name__}: No "
                               "bad pixel coordinates to write.")

        # Since numpy 2.0 (numpy/numpy/pull/22449), numpy data types
        # are printed out together with scalars. This means that
        # str(np.int(x)) == 'np.int(x)' rather than 'x'. Reverting the
        # behavior to the v1.25 one ensures that we can then read the
        # information back using ast.
        print_np_scalar_as_number = (
            np.printoptions(legacy='1.25') if NUMPY2_OR_LATER
            else nullcontext()
            )
        with print_np_scalar_as_number:
            # Prepare the column contents
            columns = [
                ['Bad pixels',
                 *(f"{tuple(b)}" for b in self.__bad_coords),
                 *(f"{tuple(b)}" for b in self.__uncorrectable)],
                ['Replacements',
                 *(f"{tuple(r)}" for r in self.__replacements),
                 *(['NaN']*len(self.__uncorrectable))]
                ]

        # Find padding lengths:
        bad_width = max(len('# Bad pixels'), *(len(row) for row in columns[0]))
        repl_width = max(len(row) for row in columns[1])

        # And format the lines
        columns[0][0] = f"#{columns[0][0]:>{bad_width-1}}"
        columns[0][1:] = [f"{row:>{bad_width}}" for row in columns[0][1:]]
        columns[1] = [f"{row:>{repl_width}}" for row in columns[1]]
        lines = '\n'.join(f'{bad}, {repl}' for bad, repl in zip(*columns))

        filepath = Path(filepath)
        if not filepath.exists():
            filepath.mkdir(parents=True)
        filename = filepath / f"{self.__base_name}.badpx"

        width, height, *_ = self.__camera.image_info
        n_bad = self.n_bad_pixels_sensor
        n_uncorr = self.n_uncorrectable_sensor
        comment = (f"# Bad pixels file for camera {self.__camera.name}. \n"
                   f"# Total number of bad pixels: {n_bad}, i.e., "
                   f"{100*n_bad/(width*height):.2f}% of the sensor.\n")
        if n_uncorr:
            comment += f"# Of these, {n_uncorr} cannot be corrected.\n"
        comment += ("# Bad pixel coordinates are (row number, column "
                    "number), i.e., (y, x).\n#\n")

        with open(filename, 'w', encoding='utf-8') as bad_px_file:
            bad_px_file.write(comment)
            bad_px_file.write(lines)

    @property
    def __has_info(self):
        """Return whether bad pixel information is present."""
        # This means correctly __init__ed, or read from a valid file
        return np.all(np.isfinite(self.__bad_coords))

    def __bool__(self):
        """Return the truth value of self."""
        return bool(self.__has_info and self.n_bad_pixels_sensor)


class BadPixelsMaxMinStorage:
    """Container for maximum and minimum of bad-pixel frames."""

    def __init__(self, height, width, dtype):
        """Initialize bad pixel peak-to-peak storage class."""
        self._count = 0
        if dtype not in (np.uint8, np.uint16, np.uint32, np.uint64):
            raise TypeError('frame pixel value type is not a np.uint.')
        self._max = np.zeros((height, width), dtype=dtype)
        max_value = np.iinfo(dtype).max
        self._min = np.full((height, width), max_value, dtype=dtype)

    @property
    def max(self):
        """Return pixel maxima."""
        return self._max

    @property
    def min(self):
        """Return pixel minima."""
        return self._min

    @property
    def range_(self):
        """Return pixel maximum minus pixel minimum."""
        return self._max - self._min

    def add_frame(self, frame):
        """Add frame to frame storage."""
        self._count += 1
        self._max = np.maximum(self._max, frame)
        self._min = np.minimum(self._min, frame)


class BadPixelsSumStorage:
    """Class for containing the sum of bad pixel frames."""

    def __init__(self, height, width):
        """Initialize bad pixel sum storage class."""
        self._count = 0
        self._sum = np.zeros((height, width), dtype=np.uint64)
        self._sum_squares = np.zeros((height, width), dtype=np.uint64)

    @property
    def frame_sum(self):
        """Return the summed up frames."""
        return self._sum

    def add_frame(self, frame):
        """Add frame to frame sum and sum of squares."""
        frame = frame.astype(np.uint64, copy=False)
        self._count += 1
        self._sum += frame
        self._sum_squares += frame**2

    def mean(self):
        """Return the mean of the stored frames."""
        return self._sum / self._count

    def var(self):
        """Return the variance calculated from the stored frames."""
        mean_of_squares = self._sum_squares / self._count
        return mean_of_squares - self.mean()**2
