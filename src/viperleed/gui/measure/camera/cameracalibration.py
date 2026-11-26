"""Module cameracalibration of viperleed.gui.measure.camera.

Defines abstract CalibrationTask subclasses that are specific to
cameras.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-10-12'
__license__ = 'GPLv3+'

from abc import abstractmethod

import numpy as np
from PyQt5 import QtCore as qtc

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.calibrationtask import CalibrationTask


_INVOKE = qtc.QMetaObject.invokeMethod


class CameraCalibrationErrors(base.ViPErLEEDErrorEnum):
    """Class for camera-calibration errors."""

    DARK_FRAME_TOO_BRIGHT = (
        520,
        "Dark frame has too much intensity. Camera "
        "optics is not appropriately covered."
        )


class CameraCalibrationTask(CalibrationTask):
    """Calibration task for cameras."""

    def __init__(self, device, *args, timeout=None, **kwargs):
        """Initialize this task.

        Parameters
        ----------
        device : CameraABC
            The device to which this calibration task belongs.
            This is used as the parent of this task.
        *args : object
            Unused positional arguments that subclasses may need.
        timeout : int, or None, optional
            The maximum duration for this task, in milliseconds.
            If negative or None, this task will not time out.
            Default is None.
        triggered : bool, optional
            Whether this calibration task acquires frames by
            triggering the camera. Default is True.
        saves_images : bool, optional
            Whether this calibration task should save the images
            it acquires. Default is False.
        **kwargs : object, optional
            Other unused keyword arguments that subclasses may need.
        """
        super().__init__(device, *args, timeout=None, **kwargs)

        self.__saves_images = kwargs.get('saves_images', False)

        if kwargs.get('triggered', True):
            self._to_be_connected.extend((
                (self.camera.image_processed, self._check_and_store_frame),
                (self.camera.started, self._trigger_next_frame),
                ))
            mode = 'triggered'
        else:
            self._to_be_connected.append(
                (self.camera.frame_ready, self._check_and_store_frame)
                )
            mode = 'live'

        if self.saves_images:
            self._to_be_connected.append(
                (self.camera.image_saved, self.__on_image_saved)
                )

        self._task_settings.set('camera_settings', 'mode', mode)
        self._original['process'] = self.camera.process_info.copy()
        self._info.setText(
            f"Frames will be acquired from camera {self.camera.name}."
            )
        self._frames_done = 0       # Counter; usable in subclasses

        px_min, px_max = self.camera.intensity_limits
        self._limits = {
            'gain': self.camera.get_gain_limits(),
            'exposure': self.camera.get_exposure_limits(),
            'intensity': (px_min, px_max, px_max - px_min),
            }

    @property
    def camera(self):
        """Return the camera to which this task belongs."""
        return self.device

    @property
    def camera_is_triggered(self):
        """Whether this task uses the camera in triggered mode."""
        mode = self._task_settings['camera_settings']['mode']
        return mode == 'triggered'

    @property
    def can_trigger_camera(self):
        """Return whether it is possible to trigger the camera."""
        if not self.camera_is_triggered:
            return False
        cam = self.camera
        return not cam.busy and cam.is_running and cam.mode == 'triggered'

    @property
    def frame_info(self):
        """Return height, width and dtype for the current camera."""
        width, height, n_bytes, *_ = self.camera.image_info
        dtype = np.uint8 if n_bytes == 1 else np.uint16
        return height, width, dtype

    @property
    def missing_frames(self):
        """Return whether frames are missing to complete acquisition.

        Subclasses can override this property, typically making
        use of self._frames_done. The base-class implementation
        always returns False.

        Returns
        -------
        missing_frames : bool
            True if more frames should be acquired.
        """
        return False

    @property
    def saves_images(self):
        """Return whether this task also saves images to disk."""
        return self.__saves_images

    @qtc.pyqtSlot(np.ndarray)
    @abstractmethod
    def _check_and_store_frame(self, frame):
        """Check that frame is acceptable and store it if it is.

        This slot must be overridden in subclasses. subclasses should
        call ._trigger_next_frame at the end of the overridden method
        if .camera_is_triggered. This will proceed to acquisition of
        the next image. This method is also a common place to increment
        self._frames_done to reflect the fact that more frames have
        been acquired.

        This is the slot connected to the camera image_processed signal
        when .camera_is_triggered, otherwise to the camera frame_ready
        signal. It is called every time a new image was acquired.

        Parameters
        ----------
        frame : numpy.ndarray
            The frame just acquired.

        Returns
        -------
        None.
        """

    def mark_as_done(self, is_done):
        """Mark this task as to be done or not after clearing _frames_done."""
        self._frames_done = 0
        super().mark_as_done(is_done)

    @qtc.pyqtSlot(tuple)
    def _on_device_error(self, error):
        """Respond to a device error.

        Parameters
        ----------
        error : tuple or ViPErLEEDErrorEnum
            Error information. If a tuple, the first element
            is the error code, the second the error message.

        Returns
        -------
        silenced : bool
            Whether the error was silenced. Silenced errors
            do not cause emission of self.error_occurred.
            The current implementation silences errors that
            are related to bad-pixels if this task is to be
            performed before a bad-pixels-finding routine.

        Emits
        -----
        self.error_occurred(error)
            Unless it returns True.
        """
        # If this task is performed before a bad-pixels-finding
        # routine, we should not cause abortion if there's no
        # bad-pixels information.
        if (self.camera.is_bad_pixels_error(error)
                and self in self.camera.calibration_tasks['bad_pixels']):
            return True
        return super()._on_device_error(error)

    def restore_device(self):
        """Restore settings and other device attributes.

        This method is called each time the task is aborted
        or has finished. The original device settings can be
        accessed via self.original_settings.

        Subclasses can extend this method, and may modify the
        _original dictionary before calling super().restore_device()

        Returns
        -------
        None.
        """
        super().restore_device()
        self.camera.process_info.restore_from(self._original['process'])

    @qtc.pyqtSlot()
    @abstractmethod
    def start(self):
        """Start this task.

        This method should be extended in subclasses, i.e.,
        subclasses should check that super().start() returns
        True. Only in this case they can start the task.

        Returns
        -------
        ok_to_start : bool
            Whether starting the device is possible.
        """
        if not super().start():
            return False
        fname = ''
        if self.saves_images:
            fname = base.as_valid_filename(self.name) + '_{__count__}.tif'
        self.camera.process_info.filename = fname
        return True

    def start_camera(self):
        """Start camera in the correct thread."""
        self._connect_device_signals()
        _INVOKE(self.camera, 'start')

    def trigger_camera(self):
        """Trigger camera in the correct thread."""
        _INVOKE(self.camera, 'trigger_now')

    @qtc.pyqtSlot()
    @abstractmethod
    def _trigger_next_frame(self, *_):
        """Trigger acquisition of a new frame if necessary.

        This slot must be overridden in subclasses. This is
        the slot connected with the camera .started signal if
        self.camera_is_triggered. Otherwise this method should
        not be used. It is safest to check .can_trigger_camera
        and return if this is False-y. Subclasses can also
        reset self._frames_done in here, if necessary (e.g.,
        when a new 'acquisition section' has to be started).

        Returns
        -------
        None.
        """
        if not self.can_trigger_camera:
            return

    @qtc.pyqtSlot(str)
    def __on_image_saved(self, *_):
        """Increase filename count after an image was saved."""
        self.camera.process_info.count += 1
