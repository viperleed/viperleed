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
from copy import deepcopy
from datetime import datetime
from enum import Enum
from pathlib import Path

import numpy as np
from scipy.signal import convolve2d
from scipy import ndimage

from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw)

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.camera import abc
from viperleed.gui.measure.camera import tifffile
from viperleed.gui.measure.widgets.camerawidgets import CameraViewer


N_DARK = 100  # Number of dark frames used for finding bad pixels


def _report_progress(func):
    """Decorate calculation functions."""

    def _wrapper(self, *args, **kwargs):
        """Wrap decorated function."""
        # Names are of the form "find_WHAT_other_stuff"
        _, what, *_ = func.__name__.split('_')
        section = _FinderSection.from_short_name(what)
        tasks = section.n_tasks
        self.progress_occurred.emit(section.long_name, section.value,
                                    section.n_sections, 0, tasks)
        ret_val = func(self, *args, **kwargs)
        self.progress_occurred.emit(section.long_name, section.value,
                                    section.n_sections, tasks, tasks)
        return ret_val

    return _wrapper


class BadPixelsFinderErrors(base.ViPErLEEDErrorEnum):
    """Class for bad-pixel-finder errors."""
    DARK_FRAME_TOO_BRIGHT = (210,
                             "Dark frame has too much intensity. Camera "
                             "optics is not appropriately covered.")
    FLAT_FRAME_WRONG_LIGHT = (211,
                              "Flat frame is too {}. Cannot automatically "
                              "adjust exposure time and gain.")


class _FinderSection(Enum):
    """Enumeration class for bad-pixel finder sections."""

    ACQUIRE_DARK_SHORT, ACQUIRE_DARK_LONG, ACQUIRE_FLAT = 1, 2, 3
    CALCULATE_FLICKERY, CALCULATE_HOT, CALCULATE_DEAD = 4, 5, 6
    CALCULATE_BAD = 7
    DONE = 8

    @property
    def long_name(self):
        """Return a long description for this section."""
        if self.name == 'DONE':
            return "Done."
        operation, what, *_ = self.name.split('_')

        if operation == "ACQUIRE":
            txt = f"Acquiring {what.lower()} frame"
            txt += "s " if what == "DARK" else " "
            return txt + "with {:.1f} ms exposure..."
        # Calculating
        return f"Calculating coordinates of {what.lower()} pixels..."

    @property
    def n_sections(self):
        """Return the total number of sections."""
        return 8

    @property
    def n_tasks(self):
        """Return the number of tasks in this section."""
        if 'DARK' in self.name:
            return N_DARK
        return 1

    @classmethod
    def from_short_name(cls, name):
        """Return an instance from a short name."""
        if name.lower() in ('dark-short', 'dark-long', 'flat'):
            elem = f"ACQUIRE_{name.replace('-', '_').upper()}"
        elif name.lower() == 'done':
            elem = "DONE"
        else:
            elem = f"CALCULATE_{name.upper()}"
        return getattr(cls, elem)


class BadPixelsFinder(qtc.QObject):
    """Class for finding bad pixels."""

    error_occurred = qtc.pyqtSignal(tuple)
    done = qtc.pyqtSignal()
    aborted = qtc.pyqtSignal()
    progress_occurred = qtc.pyqtSignal(
        str,  # current operation
        int,  # current operation index (1-based)
        int,  # total no. operations
        int,  # current progress
        int   # total tasks in this section
        )

    def __init__(self, camera, parent=None):
        """Initialize object."""
        if not isinstance(camera, abc.CameraABC):
            raise TypeError(
                f"{self.__class__.__name__}: invalid type "
                f"'{type(camera).__name__}' for camera argument. "
                "Must be a concrete subclass of CameraABC."
                )
        super().__init__(parent=parent)

        self.__camera = camera
        self.__viewer = CameraViewer(self.__camera, stop_on_close=False,
                                     visible=False)

        self.__camera.error_occurred.connect(self.error_occurred)
        self.__camera.image_processed.connect(self.__check_and_store_frame)
        self.__camera.camera_busy.connect(self.__trigger_next_frame)

        self.__frames_done = 0
        self.__adjustments = 0  # used for progress reporting

        # Keep track of the old camera settings, as we will need to
        # change them on the fly later. Also, store the previous
        # visibility state of the viewer. This may not match the
        # False above, as there can only be one camera viewer per
        # each camera object.
        # Notice the deep copies of camera settings, which make sure
        # that camera.settings, self.__original['settings'] and
        # self.__new_settings are distinct objects. This is necessary
        # because the camera does some checks for quantities changed
        # to decide whether to recalculate stuff.
        self.__original = {'settings': deepcopy(self.__camera.settings),
                           'process': self.__camera.process_info.copy(),
                           'show_auto': self.__viewer.show_auto,
                           'was_visible': self.__viewer.isVisible()}

        # Now set new settings that will be used
        # throughout the rest of the processing.
        self.__viewer.show_auto = False
        self.__viewer.hide()
        self.__new_settings = deepcopy(self.__camera.settings)

        _, (max_roi_w, max_roi_h), *_ = self.__camera.get_roi_size_limits()
        full_roi = f"(0, 0, {max_roi_w}, {max_roi_h})"
        self.__new_settings.set('camera_settings', 'roi', full_roi)
        self.__new_settings.set('camera_settings', 'binning', '1')
        self.__new_settings.set('camera_settings', 'mode', 'triggered')
        self.__new_settings.set('measurement_settings', 'n_frames', '1')
        self.__new_settings.set('camera_settings', 'bad_pixels_path', '')

        # And set them into the camera. This stops
        # the camera if it is currently running.
        self.__camera.settings = deepcopy(self.__new_settings)
        self.__camera.process_info.filename = ''

        # Clear any bad pixel info in the camera (frames should be raw)
        bad_old = self.__camera.bad_pixels
        if bad_old:
            bad_old.clear()

        width, height, n_bytes, _ = self.__camera.image_info
        dtype = np.uint8 if n_bytes == 1 else np.uint16

        # __imgs contains: the short- and long-exposure movies
        # (N_DARK frames each) for the 'dark' frame (i.e., camera
        # with a cap on) and one frame for a 'flat field' (e.g.,
        # white paper right in front of the lens).
        self.__imgs = {
            "dark-short": np.zeros((N_DARK, height, width), dtype=dtype),
            "dark-long": np.zeros((N_DARK, height, width), dtype=dtype),
            "flat": np.zeros((height, width), dtype=dtype)
            }
        self.__current_section = "dark-short"

        # __badness contains a measure of the badness of
        # each pixel. 'Normal' pixels have badness 0, bad
        # pixels have badness >= 6, very good pixels < 0.
        self.__badness = np.zeros((height, width), dtype=float)

        # __bad_pixels will contain info about bad
        # pixel coordinates and their replacements
        self.__bad_pixels = BadPixels(self.__camera)

        name = camera.name

        self.__info = qtw.QMessageBox(parent=parent)
        self.__info.setWindowTitle(f"Bad pixel detection for {name}")
        self.__info.setIcon(self.__info.Information)
        self.__info.setText(
            f"Frames will be acquired from {name} camera."
            )

        self.__calc_thread = _BadPixCalculationThread(self)

    @property
    def short_exposure(self):
        """Return a short exposure time in milliseconds."""
        min_exposure, _ = self.__camera.get_exposure_limits()
        return max(20, min_exposure)

    @property
    def long_exposure(self):
        """Return a long exposure time in milliseconds."""
        _, max_exposure = self.__camera.get_exposure_limits()
        return min(1000, max_exposure)

    @property
    def large_gain(self):
        """Return a large gain in decibel."""
        _, max_gain = self.__camera.get_gain_limits()
        return min(20, max_gain)

    @property
    def missing_frames(self):
        """Return whether frames are missing for the current section."""
        if 'dark' in self.__current_section:
            todo = self.__imgs[self.__current_section].shape[0]
        else:
            todo = 1
        return self.__frames_done < todo

    def begin_acquiring(self):
        """Start acquiring images."""
        if "dark" in self.__current_section:
            self.__info.setInformativeText(
                "A 'dark' movie will be acquired.\nMake sure no light "
                "enters the camera.\nFor example, you can close the "
                "lens with its cap.\n\nPress 'OK' when ready."
                )
        else:
            self.__info.setInformativeText(
                "A 'flat' frame will be acquired.\nMake sure the camera "
                "is viewing a featureless background with uniform "
                "illumination.\nFor example, you can place a white "
                "piece of paper right in front of the lens. Normal room "
                "light should be uniform enough.\n\nPress 'OK' when ready."
                )
        if 'long' not in self.__current_section:
            self.__info.exec_()

        if 'dark' in self.__current_section:
            # Here we set fixed gain and exposure time
            if 'short' in self.__current_section:
                exposure = self.short_exposure
            else:
                exposure = self.long_exposure

            self.__new_settings.set("measurement_settings", "exposure",
                                    str(exposure))
            self.__new_settings.set("measurement_settings", "gain",
                                    str(self.large_gain))

        else:
            # Here we use the large gain and a somewhat short
            # exposure to begin with, but these are adjusted when
            # frames come back to have an average image intensity
            self.__new_settings.set("measurement_settings", "gain",
                                    str(self.large_gain))
            self.__new_settings.set("measurement_settings", "exposure",
                                    str(self.short_exposure))

        # Apply the new settings and start the camera.
        # When the camera is done preparing it will emit
        # image_processed, which will trigger storage of
        # frames, and a camera_busy, which will trigger
        # acquisition of the next frame, if needed.
        self.__camera.settings = deepcopy(self.__new_settings)
        self.__report_acquisition_progress()
        self.__camera.start()

    def abort(self):
        """Abort finding bad pixels."""
        self.__restore_settings()  # Also stops
        self.__camera.update_bad_pixels()  # Restore old file
        self.__frames_done = 0
        self.aborted.emit()

    def find(self):
        """Find bad pixels in the camera."""
        if self.__camera.busy:
            # Cannot start detecting bad pixels as long as the
            # camera is busy (esp. right after it is started)
            base.emit_error(self, abc.CameraErrors.UNSUPPORTED_WHILE_BUSY,
                            'find bad pixels')
            return

        self.__current_section = "dark-short"
        self.begin_acquiring()

    def __adjust_exposure_and_gain(self, frame):
        """Use a frame to estimate optimal exposure for flat frames."""
        old_exposure = self.__camera.exposure
        old_gain = self.__camera.gain
        pix_min, pix_max = self.__camera.intensity_limits
        intensity_range = pix_max - pix_min
        exposure_min, exposure_max = self.__camera.get_exposure_limits()

        intensity = frame.mean() - pix_min
        new_exposure = old_exposure * 0.5*intensity_range / intensity
        if new_exposure < exposure_min:
            # Decrease the gain
            new_gain = old_gain + 20*np.log10(new_exposure/exposure_min)
            new_exposure = exposure_min
        elif new_exposure > exposure_max:
            # Increase the gain (this is very unlikely to happen)
            new_gain = old_gain + 20*np.log10(new_exposure/exposure_max)
            new_exposure = exposure_max
        else:
            new_gain = old_gain

        # Now check that the gain is OK: the only cases in which
        # it cannot be OK is if we need very short or very long
        # exposures (otherwise the gain stayed the same).
        min_gain, max_gain = self.__camera.get_gain_limits()
        if new_gain < min_gain:
            # Too much intensity
            base.emit_error(self, BadPixelsFinderErrors.FLAT_FRAME_WRONG_LIGHT,
                            'bright')
            return
        if new_gain > max_gain:
            # Too little intensity
            base.emit_error(self, BadPixelsFinderErrors.FLAT_FRAME_WRONG_LIGHT,
                            'bright')
            return
        self.__new_settings.set('measurement_settings', 'exposure',
                                f'{new_exposure:.3f}')
        self.__new_settings.set('measurement_settings', 'gain',
                                f'{new_gain:.1f}')
        self.__camera.settings = deepcopy(self.__new_settings)
        self.__adjustments += 1
        self.__report_acquisition_progress()
        self.__camera.start()

    def __check_and_store_frame(self, frame):
        """Store a new frame."""
        sec = self.__current_section
        if not self.__frame_acceptable(frame):
            if 'dark' in sec:
                base.emit_error(self,
                                BadPixelsFinderErrors.DARK_FRAME_TOO_BRIGHT)
            else:
                self.__adjust_exposure_and_gain(frame)
            return

        if 'dark' in sec:
            self.__imgs[sec][self.__frames_done, :, :] = frame
        else:
            self.__imgs[sec][:, :] = frame
        self.__frames_done += 1
        self.__report_acquisition_progress()

    def __report_acquisition_progress(self):
        """Report progress of acquisition."""
        section = _FinderSection.from_short_name(self.__current_section)
        name = section.long_name.format(self.__camera.exposure)
        self.progress_occurred.emit(name, section.value, section.n_sections,
                                    self.__frames_done + self.__adjustments,
                                    section.n_tasks + self.__adjustments)

    @_report_progress
    def find_bad_and_replacements(self):
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
        badness = self.__badness.copy()

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
            self.__camera,
            bad_coordinates=np.asarray((bad_y, bad_x)).T,
            replacement_offsets=offsets[best_offset_indices],
            uncorrectable=uncorrectable
            )

    @_report_progress
    def find_dead_pixels(self):
        """Detect dead pixels from flat frame.

        Dead pixels are those whose intensity is much smaller
        than the average of the neighbors (excluding flickery
        and hot pixels). This method adds to the badness the
        ratio I_neighbors/I_pixel - 1.

        Return
        ------
        None.
        """
        pix_min, pix_max = self.__camera.intensity_limits
        intensity_range = pix_max - pix_min

        flat = self.__imgs['flat'] - pix_min

        # We have to exclude the already-detected bad pixels.
        # For this to happen, we will (1) add a bit to all
        # pixels (so that no pixel is zero), (2) set the bad
        # ones to zero, (3) convolve with an appropriate
        # averaging kernel, also counting nonzero pixels (i.e.,
        # the good ones). This way, bad ones do not factor in
        # the neighbor average.

        offset = 0.1*intensity_range
        flat_bad = flat + offset
        flat_bad_ones = np.ones(flat_bad.shape)  # for counting

        bad_px_mask = self.__badness >= 6

        flat_bad[bad_px_mask] = 0
        flat_bad_ones[bad_px_mask] = 0

        # 5x5 kernel for averaging a diamond-shaped patch of neighbors
        kernel = np.array(((0, 0, 1, 0, 0),
                           (0, 1, 1, 1, 0),
                           (1, 1, 0, 1, 1),
                           (0, 1, 1, 1, 0),
                           (0, 0, 1, 0, 0)))

        neighbor_sum = convolve2d(flat_bad, kernel, mode='same')
        neighbor_count = convolve2d(flat_bad_ones, kernel, mode='same')
        neighbor_ave = neighbor_sum / neighbor_count - offset

        # Isolated pixels in a patch of bad neighbors
        # will also be marked as bad.
        neighbor_ave[neighbor_count < 0.5] = np.inf

        # As will those that have both zero neighbor average
        # and zero intensity (dead in a cluster of bad ones)
        delta_badness = neighbor_ave / flat - 1
        delta_badness[np.isnan(delta_badness)] = np.inf

        self.__badness += delta_badness

    @_report_progress
    def find_flickery_pixels(self):
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
        long_flicker = self.__imgs['dark-long'].var(axis=0)
        short_flicker = self.__imgs['dark-short'].var(axis=0)

        # '-1' is such that normally flickery pixels have badness == 0
        self.__badness = (long_flicker / long_flicker.mean() - 1
                          + short_flicker / short_flicker.mean() - 1)

    @_report_progress
    def find_hot_pixels(self):
        """Set badness of hot pixels to infinity."""
        pix_min, pix_max = self.__camera.intensity_limits
        intensity_range = pix_max - pix_min

        dark_ave = self.__imgs['dark-long'].mean(axis=0) - pix_min
        hot_pixels = dark_ave >= 0.25*intensity_range

        self.__badness[hot_pixels] = np.inf

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
        pixel_min, pixel_max = self.__camera.intensity_limits
        intensity_range = pixel_max - pixel_min
        intensity = frame.mean() - pixel_min
        if 'dark' in self.__current_section:
            return intensity < 0.2*intensity_range
        return 0.45*intensity_range < intensity < 0.55*intensity_range

    def __restore_settings(self):
        """Restore the original settings of the camera."""
        self.__camera.settings = self.__original['settings']
        self.__camera.process_info.restore_from(self.__original['process'])
        self.__viewer.show_auto = self.__original['show_auto']
        if self.__original['was_visible']:
            self.__viewer.show()

    def save_and_cleanup(self):
        """Save a bad pixels file and finish."""
        bp_path = self.__original['settings'].get("camera_settings",
                                                  "bad_pixels_path",
                                                  fallback='')
        self.__bad_pixels.write(bp_path)
        self.__restore_settings()
        done = _FinderSection.DONE
        self.progress_occurred.emit(done.long_name, done.value,
                                    done.n_sections,
                                    done.n_tasks, done.n_tasks)
        self.done.emit()

    def __trigger_next_frame(self, *_):
        """Trigger acquisition of a new frame if necessary."""
        if self.__camera.busy or not self.__camera.is_running:
            return

        if self.missing_frames:
            self.__camera.trigger_now()
            return

        # One section over, go to the next one
        sections = tuple(self.__imgs.keys())
        next_idx = sections.index(self.__current_section) + 1
        if next_idx >= len(sections):
            # Done with all sections. Can proceed to calculations.
            # Do it in a separate threat to keep the UI responsive.
            self.__calc_thread.start()
        else:
            self.__current_section = sections[next_idx]
            self.__frames_done = 0
            self.__adjustments = 0
            self.begin_acquiring()


class _BadPixCalculationThread(qtc.QThread):

    def __init__(self, finder, parent=None):
        super().__init__(parent=parent)
        self.__finder = finder

    def run(self):
        """Run calculations."""
        self.__finder.find_flickery_pixels()
        self.__finder.find_hot_pixels()
        self.__finder.find_dead_pixels()
        self.__finder.find_bad_and_replacements()
        self.__finder.save_and_cleanup()


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

        if not isinstance(camera, abc.CameraABC):
            raise TypeError(f"{self.__class__.__name__}: invalid type "
                            f"{type(camera).__name__} for camera argument. "
                            "Expected a subclass of CameraABC.")

        self.__camera = camera
        self.__bad_coords = bad_coordinates
        self.__replacements = replacement_offsets
        self.__uncorrectable = uncorrectable

        self.__bad_coords_roi = bad_coordinates
        self.__replacements_roi = replacement_offsets
        self.__uncorrectable_roi = uncorrectable
        self.__datetime = datetime.now()

        # Prepare a base name used for saving bad-pixel info to disk
        cam_name = self.__camera.name.replace(' ', '_')
        now = self.__datetime.strftime("%Y%m%d_%H%M%S")
        self.__base_name = f"{cam_name}_{now}"

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
        neighbors (in the allowed replacement directions).

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

        if no_roi:
            self.__bad_coords_roi = self.__bad_coords
            self.__replacements_roi = self.__replacements
            self.__uncorrectable_roi = self.__uncorrectable
            return

        top_x, top_y, width, height = self.__camera.roi

        # Recalculate coordinates based on the new origin
        # NB: the bad coordinates are (row#, col#)!
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
        if self.__uncorrectable is None or not self.__uncorrectable.size:
            return
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
        filepath = Path(filepath)
        cam_name = self.__camera.name.replace(' ', '_')

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
            base_name, filename = sorted(((f.stem, f) for f in bad_px_files),
                                         reverse=True)[idx]
        except IndexError:
            self.__bad_coords = np.ones((1, 2)) * np.nan
            raise FileNotFoundError(
                f"{self.__class__.__name__}: found only one "
                f"bad pixels file for camera {self.__camera.name} "
                f"in directory {filepath.resolve()}. Cannot load "
                "the second-most recent file."
                ) from None

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
        if not self.__has_info:
            raise RuntimeError(f"{self.__class__.__name__}: No "
                               "bad pixel information present.")

        # Prepare the image
        width, height, *_ = self.__camera.image_info
        mask = np.zeros((height, width), dtype='uint8')

        if len(self.uncorrectable):
            bad_y, bad_x = self.uncorrectable.T
            mask[bad_y, bad_x] = 255

        date_time = self.__datetime.strftime("%Y:%m:%d %H:%M:%S")
        filename = Path(filepath, f"{self.__base_name}_uncorr_mask.tiff")
        tiff = tifffile.TiffFile.from_array(
            mask,
            comment=f"Uncorrectable bad pixels mask for {self.__camera.name}",
            model=self.__camera.name,
            date_time=date_time
            )
        tiff.write(filename)

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
        if not self.__has_info:
            raise RuntimeError(f"{self.__class__.__name__}: No "
                               "bad pixel information present.")

        # Prepare the image
        width, height, *_ = self.__camera.image_info
        mask = np.zeros((height, width), dtype='uint8')

        if len(self.bad_pixel_coordinates):
            bad_y, bad_x = self.bad_pixel_coordinates.T
            mask[bad_y, bad_x] = 255

        date_time = self.__datetime.strftime("%Y:%m:%d %H:%M:%S")
        filename = Path(filepath, f"{self.__base_name}_bad_mask.tiff")
        tiff = tifffile.TiffFile.from_array(
            mask,
            comment=f"Bad pixels mask for {self.__camera.name}",
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
        return np.all(np.isfinite(self.__bad_coords))

    def __bool__(self):
        """Return the truth value of self."""
        if not self.__has_info:
            return False
        if not self.bad_pixel_coordinates.size:
            return False
        return True
