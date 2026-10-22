"""Module imagingsourcecalibration of viperleed.gui.measure.camera.

Defines CalibrationTask subclasses to be used with cameras from
The Imaging Source.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-10-12'
__license__ = 'GPLv3+'

from copy import deepcopy

from PyQt5 import QtCore as qtc
import numpy as np

from viperleed.gui.measure.camera.cameracalibration import (
    CameraCalibrationErrors,
    CameraCalibrationTask,
    )
from viperleed.gui.measure.classes.calibrationtask import (
    CalibrationTaskOperation,
    )
from viperleed.gui.measure.hardwarebase import ViPErLEEDErrorEnum


_INVOKE = qtc.QMetaObject.invokeMethod


class ImagingSourceCalibrationError(ViPErLEEDErrorEnum):
    """Class for bad-pixel-finder errors."""

    MAXIMUM_DARK_LEVEL_REACHED = (
        530,
        "Failed to calibrate dark level for camera {}. Camera intensity "
        "histogram is cut even at the largest dark level settings."
        )


_N_DARK_LEVEL = 3  # No. of changes of dark level


class _DarkLevelOperation(CalibrationTaskOperation):
    """Enumeration of dark-level calibration operations."""

    ENSURE_CORRECT_MINIMUM = 1
    ACQUIRE_IMAGES, VERIFY_DARK_LEVEL, DONE = 2, 3, 4

    @property
    def long_name(self):
        """Return a descriptive name for this operation."""
        if self is _DarkLevelOperation.DONE:
            return "Done."
        if self is _DarkLevelOperation.ENSURE_CORRECT_MINIMUM:
            txt = "Acquiring dark frame with {exposure:.2f} ms exposure..."
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


class DarkLevelCalibration(CameraCalibrationTask):
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

        self.__current_section = _DarkLevelOperation.first()
        self.__section_settings_dict = self.__prepare_task_settings()

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
        intensity = mean_ - 6*stdev

        if mean_ > 0.5*int_range:  # Too much for a dark frame
            self.emit_error(CameraCalibrationErrors.DARK_FRAME_TOO_BRIGHT)
            return

        this_section = self.__current_section
        assert this_section <= _DarkLevelOperation.VERIFY_DARK_LEVEL
        if this_section is _DarkLevelOperation.ENSURE_CORRECT_MINIMUM:
            # Will soon start acquiring the actual frames
            self.__update_limits_and_continue()
            return

        if (this_section is _DarkLevelOperation.VERIFY_DARK_LEVEL
                and intensity > 0):
            # Dark level successfully verified
            self.__finish()
            return

        if intensity > 0:
            # Store the current results
            _dark = self.__dark_delta
            self.__img_data['sum_dark'] += _dark
            self.__img_data['sum_int'] += intensity
            self.__img_data['sum_dark_int'] += _dark*intensity
            self.__img_data['sum_dark_sq'] += _dark**2
            self._frames_done += 1

        if (this_section is _DarkLevelOperation.ACQUIRE_IMAGES
                and not self.missing_frames):
            self.__current_section = this_section.next_()

        try:
            self.__calculate_next_dark(intensity)
        except RuntimeError:
            # We reached the maximum. Error already reported.
            return

        self._trigger_next_frame()
        self.__report_progress()

    @qtc.pyqtSlot()
    def continue_(self, *_):
        """Start camera after user confirmation."""
        self._task_settings.read_dict(self.__section_settings)
        self.__report_progress()

        self.update_device_settings()
        self.start_camera()

    @qtc.pyqtSlot()
    def start(self):
        """Start calibrating the dark level."""
        ok_to_start = super().start()
        if ok_to_start:
            self.__current_section = _DarkLevelOperation.first()
            if self.camera.has_zero_minimum:
                # If we already know, skip the first segment
                self.__current_section = self.__current_section.next_()
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

        if self.__current_section >= _DarkLevelOperation.ACQUIRE_IMAGES:
            self.camera.settings.set('camera_settings', 'black_level',
                                     str(self.dark_level))
            _INVOKE(self.camera, 'set_black_level')

        self.trigger_camera()

    @property
    def __section_settings(self):
        """Return a dictionary of settings for the current section."""
        return self.__section_settings_dict.get(self.__current_section, {})

    def __calculate_next_dark(self, intensity):
        """Appropriately change dark level and scaling."""
        if self.__current_section is _DarkLevelOperation.VERIFY_DARK_LEVEL:
            # We have enough frames, try the new dark level
            sum_x, sum_y, sum_xy, sum_xx = self.__img_data.values()
            intercept = round(
                (sum_x * sum_xy - sum_y * sum_xx)
                / (self._frames_done * sum_xy - sum_x * sum_y)
                )
            self.__dark_delta = round(intercept                                 # TODO: here there's potential for delta < 0 --> how to handle? Perhaps try one more time, then fail?
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
            self.emit_error(
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

    def __prepare_task_settings(self):
        """Prepare settings to use for the acquisition sections."""
        min_exposure, max_exposure = self._limits['exposure']
        min_gain, max_gain = self._limits['gain']
        min_dark, _ = self._limits['dark']
        _, (max_roi_w, max_roi_h), *_ = self.camera.get_roi_size_limits()
        full_roi = f"(0, 0, {max_roi_w}, {max_roi_h})"

        base_settings = {  # Common to both sections
            'camera_settings': {'roi': full_roi, 'binning': '1'},
            'measurement_settings': {'n_frames': '1'}
            }
        min_set, run_set = deepcopy(base_settings), base_settings

        min_set['camera_settings']['black_level'] = str(min_dark)
        min_set['measurement_settings']['gain'] = str(min_gain)
        min_set['measurement_settings']['exposure'] = str(min_exposure)

        run_set['camera_settings']['black_level'] = str(self.dark_level)
        min_set['measurement_settings']['gain'] = str(min(20, max_gain))
        run_set['measurement_settings']['exposure'] = str(
            min(500, max_exposure)
            )

        return {_DarkLevelOperation.ENSURE_CORRECT_MINIMUM: min_set,
                _DarkLevelOperation.ACQUIRE_IMAGES: run_set}

    def __report_progress(self):
        """Report info about the current state."""
        section = self.__current_section
        n_steps = section.n_steps

        steps_done = 0
        if section is _DarkLevelOperation.ACQUIRE_IMAGES:
            steps_done = self._frames_done
        elif section is _DarkLevelOperation.VERIFY_DARK_LEVEL:
            steps_done = self._frames_done - _N_DARK_LEVEL
            n_steps += steps_done

        _settings = self._task_settings
        exposure = _settings.getfloat("measurement_settings", "exposure")
        section.set_format(exposure=exposure, dark=self.dark_level)
        self.progress_occurred.emit(self.progress_name,
                                    *section, steps_done, n_steps)

    def __update_limits_and_continue(self):
        """Store new intensity limits, then go to next section."""
        # This is called after the first frame arrives in section
        # _DarkLevelOperation.ENSURE_CORRECT_MINIMUM.
        # Now the camera has updated its has_zero_minimum,
        # and thus returns the correct stuff in .intensity_limits
        px_min, px_max = self.camera.intensity_limits
        self._limits['intensity'] = (px_min, px_max, px_max - px_min)

        self.__current_section = self.__current_section.next_()
        self._task_settings.read_dict(self.__section_settings)
        self.__report_progress()

        self.update_device_settings()
        self.start_camera()
