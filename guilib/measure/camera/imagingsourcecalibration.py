"""Module imagingsourcecalibration of viperleed.guilib.measure.camera.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2022-10-12
Author: Michele Riva

Defines CalibrationTask subclasses to be used with cameras from
The Imaging Source.
"""

from PyQt5 import QtCore as qtc
import numpy as np

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.camera import cameracalibration as _calib
from viperleed.guilib.measure.classes.calibrationtask import (
    CalibrationTaskOperation
    )


_INVOKE = qtc.QMetaObject.invokeMethod


class ImagingSourceCalibrationError(base.ViPErLEEDErrorEnum):
    """Class for bad-pixel-finder errors."""
    MAXIMUM_DARK_LEVEL_REACHED = (
        530,
        "Failed to calibrate dark level for camera {}. Camera intensity "
        "histogram is cut even at the largest dark level settings."
        )


_N_DARK_LEVEL = 3  # No. of changes of dark level


class _DarkLevelOperation(CalibrationTaskOperation):
    """Enumeration of dark-level calibration operations."""
    ACQUIRE_IMAGES, VERIFY_DARK_LEVEL, DONE = 1, 2, 3

    @property
    def long_name(self):
        """Return a descriptive name for this operation."""
        if self is _DarkLevelOperation.DONE:
            return "Done."
        if self is _DarkLevelOperation.ACQUIRE_IMAGES:
            txt = "Acquiring dark frames with {exposure:.0f} ms exposure "
            txt += "(dark level = {dark})..."
        elif self is _DarkLevelOperation.VERIFY_DARK_LEVEL:
            txt = "Verifying new dark level ({dark})..."
        return self._format_long_name(txt)

    @property
    def n_steps(self):
        """Return the number of steps in this operation."""
        if self is _DarkLevelOperation.ACQUIRE_IMAGES:
            return _N_DARK_LEVEL
        return 1


class DarkLevelCalibration(_calib.CameraCalibrationTask):
    """Calibration task that determines an appropriate dark level.

    The dark level is considered appropriate if the intensity
    histogram of an dark-frame image is somewhat centred. This
    means, in practice, that there is a reasonable amount of
    noise in the camera.

    The algorithm is as follows:
    - acquire images with a relatively long exposure time (~1s)
      and a mid-range gain (~20dB)
    - Take one image at dark = min + delta, with delta = (bk_range)//4.
    - Calculate img_mean - px_min - 6*img_std. If this is positive,
      keep the value and use next_dark = min + delta dark // 2,
      otherwise scrap the value and try again with
      next_dark = min + delta * 2; decrease scaling factor (==2)
      by 1.25; repeat. Worst-case scenario: reach max and still
      negative --> camera problem.
    - Once we have two non-negative values and their dark levels,
      find dark level at which value == 0 by linear fit. The
      intercept is at (<x><xy> - <y><x**2>)/(<xy> - <x><y>).
    - Pick dark_level = 1.05*ceil(intercept) and check if the new
      level works (with the 6*sigma criterion). If not, increase
      fraction to 1.1, then 1.2, then 1.4, etc...
    """

    name = "Dark-level calibration"
    progress_name = 'Finding optimal dark level...'
    description = (
        "Calibrate the dark level (also known as black level, or "
        "brightness) of this camera. This means finding an appropriate "
        "value of photon intensity to be used as the minimum level for "
        "the camera sensor. Pixels receiving less than this intensity "
        "will have a value equal to the minimum value that the ADC can "
        "provide. Running this calibration will normalize the level of "
        "noise in images."
        )

    def __init__(self, camera):
        """Initialize this task."""
        super().__init__(camera)

        try:
            drk_min, drk_max = camera.get_black_level_limits()
        except camera.exceptions:
            # Probably camera is not connected. Try
            # once more, then maybe error out again
            camera.connect_()
            drk_min, drk_max = camera.get_black_level_limits()

        self._limits['dark'] = drk_min, drk_max - drk_min
        self.__img_data = {      # Info to compute intercept
            'sum_dark': 0,       # sum(xi)
            'sum_int': 0,        # sum(yi)
            'sum_dark_int': 0,   # sum(xi yi)
            'sum_dark_sq': 0,    # sum(xi**2)
            }
        self.__scales = {     # scaling factors
            'dark': 2.0,      # for __dark_delta
            'intercept': .05  # for picking the final dark level
            }
        _, dark_range = self._limits['dark']
        self.__dark_delta = dark_range // 4

        _, max_exposure = self._limits['exposure']
        _, max_gain = self._limits['gain']
        _, (max_roi_w, max_roi_h), *_ = camera.get_roi_size_limits()
        full_roi = f"(0, 0, {max_roi_w}, {max_roi_h})"

        self._task_settings.set('camera_settings', 'roi', full_roi)
        self._task_settings.set('camera_settings', 'binning', '1')
        self._task_settings.set('camera_settings', 'black_level',
                                str(self.dark_level))
        self._task_settings.set('measurement_settings', 'n_frames', '1')
        self._task_settings.set("measurement_settings", "exposure",
                                str(min(500, max_exposure)))
        self._task_settings.set("measurement_settings", "gain",
                                str(min(20, max_gain)))

    @property
    def dark_level(self):
        """Return the dark level currently used."""
        _min, _range = self._limits['dark']
        return _min + min(self.__dark_delta, _range)

    @property
    def missing_frames(self):
        """Return whether frames are missing to complete acquisition."""
        return self._frames_done < _N_DARK_LEVEL

    @qtc.pyqtSlot(np.ndarray)
    def _check_and_store_frame(self, frame):
        """Check that mean - 6*sigma > pix_min."""
        mean_, stdev = frame.mean(), frame.std(ddof=1)
        int_min, _, int_range = self._limits['intensity']
        mean_ -= int_min

        if mean_ > 0.5*int_range:  # Too much for a dark frame
            base.emit_error(
                self, _calib.CameraCalibrationErrors.DARK_FRAME_TOO_BRIGHT
                )
            return

        # Check if we just have done enough frames (see also below)
        just_finished = self.missing_frames

        intensity = mean_ - 6*stdev
        if intensity > 0:
            # Store the current results
            _dark = self.__dark_delta
            self.__img_data['sum_dark'] += _dark
            self.__img_data['sum_int'] += intensity
            self.__img_data['sum_dark_int'] += _dark*intensity
            self.__img_data['sum_dark_sq'] += _dark**2
            self._frames_done += 1

        try:
            self.__calculate_next_dark(intensity)
        except RuntimeError:
            # We reached the maximum. Error already reported.
            return

        # just_finished: We were missing frames before, but we're not
        # missing frames now. Will check the prospective dark level
        # by acquiring one more image.
        just_finished &= not self.missing_frames

        if self.missing_frames or just_finished or intensity <= 0:
            self._trigger_next_frame()
            self.__report_progress()
        else:
            self.__finish()

    @qtc.pyqtSlot()
    def continue_(self, *_):
        """Start camera after user confirmation."""
        self.__report_progress()
        self.start_camera()

    @qtc.pyqtSlot()
    def start(self):
        """Start calibrating the dark level."""
        ok_to_start = super().start()
        if ok_to_start:
            self.set_info_text(
                "A 'dark' movie will be acquired.<br>Make sure no light "
                "enters the camera.<br>For example, you can <b>close the "
                "lens with its cap.</b>"
                )
            self.show_info()
        return ok_to_start

    @qtc.pyqtSlot()
    def _trigger_next_frame(self, *_):
        """Trigger acquisition of a new frame if necessary."""
        if not self.can_trigger_camera:
            return
        self.camera.settings.set('camera_settings', 'black_level',
                                 str(self.dark_level))
        _INVOKE(self.camera, 'set_black_level')
        self.trigger_camera()

    def __calculate_next_dark(self, intensity):
        """Appropriately change dark level and scaling."""
        if not self.missing_frames:
            # We have enough frames, try the new dark level
            sum_x, sum_y, sum_xy, sum_xx = self.__img_data.values()
            intercept = round(
                (sum_x * sum_xy - sum_y * sum_xx)
                / (self._frames_done * sum_xy - sum_x * sum_y)
                )
            self.__dark_delta = round(intercept
                                      * (1 + self.__scales['intercept']))
            if intensity < 0:
                self.__scales['intercept'] *= 2
            return

        # More frames to be acquired
        if intensity > 0:
            # Reduce dark level
            self.__dark_delta /= self.__scales['dark']
        else:
            # Increase dark level, and reduce scaling
            # so next time we scale less (up or down)
            # and we never end up at the same exact value
            self.__dark_delta *= self.__scales['dark']
            self.__scales['dark'] *= 0.8

        self.__dark_delta = round(self.__dark_delta)
        _, _drk_range = self._limits['dark']
        if self.__dark_delta > _drk_range:
            # We have increased dark level all the way up
            # but we're still not getting anything reasonable
            base.emit_error(
                self,
                ImagingSourceCalibrationError.MAXIMUM_DARK_LEVEL_REACHED,
                self.camera.name
                )
            raise RuntimeError

    def __finish(self):
        """Store new dark level and stop camera."""
        _settings = self._original['settings']
        _settings.set('camera_settings', 'black_level',
                      str(self.camera.black_level))
        _settings.update_file()
        done = _DarkLevelOperation.DONE
        self.progress_occurred.emit(self.progress_name, *done,
                                    done.n_steps, done.n_steps)
        self.mark_as_done(True)  # Also stops camera

    def __report_progress(self):
        """Report info about the current state."""
        steps_done = self._frames_done
        if self.missing_frames:
            section = _DarkLevelOperation.ACQUIRE_IMAGES
            n_steps = section.n_steps
        else:
            section = _DarkLevelOperation.VERIFY_DARK_LEVEL
            steps_done -= _N_DARK_LEVEL
            n_steps = section.n_steps + steps_done
        section.set_format(exposure=self.camera.exposure,
                           dark=self.camera.black_level)
        self.progress_occurred.emit(self.progress_name,
                                    *section, steps_done, n_steps)
