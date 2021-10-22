"""Module cameraabc of viperleed.?????.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-12
Author: Michele Riva
Author: Florian Doerr

"""

from abc import abstractmethod

import ast
import numpy as np
from PyQt5 import QtCore as qtc

from viperleed.guilib.measure.hardwarebase import (
    ViPErLEEDErrorEnum, config_has_sections_and_options, emit_error, QMetaABC
    )
from viperleed.guilib.measure.camera.imageprocess import (ImageProcessor,
                                                          ImageProcessInfo)


class CameraErrors(ViPErLEEDErrorEnum):
    """Data class for base camera errors."""

    INVALID_SETTINGS = (200,
                        "Invalid camera settings: Required "
                        "settings {!r} missing or values "
                        "inappropriate. Check configuration file.")
    MISSING_SETTINGS = (201,
                        "Camera cannot operate without settings. Load "
                        "an appropriate settings file before proceeding.")
    CAMERA_NOT_FOUND = (202, "Could not find camera {}.")

    SETTINGS_MISMATCH = (203,
                         "Different {} settings found in camera and "
                         "configuration file: camera={}, settings={}.")
    INVALID_SETTING_WITH_FALLBACK = (204,
                                     "Invalid camera settings value {} for "
                                     "setting {!r}. Using {} instead. Consider"
                                     " fixing your configuration file.")
    UNSUPPORTED_OPERATION = (205, "Cannot {} in {} mode. Switch mode before.")
    BINNING_ROI_MISMATCH = (206,
                            "Region of interest size ({} x {}) is incompatible"
                            " with binning factor {}.Reducing region of "
                            "interest to ({} x {}). A few pixels on the "
                            "lower-right corner may be removed.")
    UNSUPPORTED_WHILE_BUSY = (207, "Cannot {} while camera is busy.")


class CameraABC(qtc.QObject, metaclass=QMetaABC):
    """Abstract base class for ViPErLEED cameras."""

    camera_busy = qtc.pyqtSignal(bool)
    error_occurred = qtc.pyqtSignal(tuple)

    # frame_ready should be emitted in any callback function that the
    # camera may call. If the camera does not allow to call a callback
    # function, the reimplementation should take care of emitting a
    # frame_ready signal as soon as a new frame arrives at the PC.
    # Whenever this signal is emitted, the new frame will be processed
    # (in a separate thread).
    frame_ready = qtc.pyqtSignal(np.ndarray)

    # __process_frame is a private signal used to pass frames to
    # the processing thread
    __process_frame = qtc.pyqtSignal(np.ndarray)

    _mandatory_settings = [
            ('camera_settings', 'class_name'),
            ('camera_settings', 'device_name'),
            ('measurement_settings', 'camera_exposure'),
            ]

    # _abstract is a list of features for which setter and getter
    # methods must always be reimplemented in concrete classes
    _abstract = ('exposure', 'gain', 'mode')

    # hardware_supported_features contains a list of feature names (as
    # strings) that the camera supports at the hardware level.
    # A concrete subclass should define a read-only @property with each
    # of these names, as well as get_<name> and set_<name> methods that
    # retrieve and set, respectively, the value from/to the camera. The
    # property should read the setting from self.settings. Properties,
    # getters and setters are already implemented (some as abstract
    # methods) for the common features 'binning', 'n_frames', and 'roi'
    # The setter methods will be called after set_exposure, set_gain,
    # and set_mode, and in the same order they are listed in this list.
    hardware_supported_features = []

    def __init__(self, driver, *args, settings=None, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        driver : object
            A reference to the library/class that handles
            the communication with the camera firmware.
        *args : object
            Other unused positional arguments
        settings : dict or ConfigParser, optional
            Settings for the camera. Default is None.
            The following sections/options are mandatory:
            'camera_settings'/'class_name'
                Name of the concrete class that implements
                the abstract methods of the base ABC.
            'camera_settings'/'device_name'
                Name of the device as it appears in the device
                list (e.g., "DMK 33GX265").
            'measurement_settings'/'camera_exposure'
                Exposure time in milliseconds to be used for
                acquiring each frame.
        **kwargs : object
            Other unused keyword arguments
        """
        self.driver = driver
        self.__busy = False
        self.__settings = None
        self.set_settings(settings)

        self.process_info = ImageProcessInfo()
        self.process_info.camera = self
        self.__process_thread = qtc.QThread()

        # Number of frames accumulated for averaging
        self.__n_frames_done = 0

        self.frame_ready.connect(self.__on_frame_ready)

    @property
    def binning(self):
        """Return the binning factor from the settings.

        Returns
        -------
        binning_factor : int
            Positive integer used to reduce the size of frames.
            Binning consists in averaging pixel intensities over
            a (binning_factor x binning_factor) mesh.
        """
        try:
            binning_factor = self.settings.getint('camera_settings',
                                                  'binning', fallback=1)
        except ValueError:  # Cannot be read as int
            binning_factor = -1

        min_bin, max_bin = self.get_binning_limits()
        if binning_factor < min_bin or binning_factor > max_bin:
            emit_error(self, CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                       binning_factor, 'camera_settings/binning', 1)
            binning_factor = 1
            self.settings.set('camera_settings', 'binning', '1')
        return binning_factor

    @property
    def busy(self):
        """Return whether the camera is busy.

        Returns
        -------
        busy : bool
            True if the camera is busy. A camera is to be considered
            busy after a measurement started and until the last frame
            requested arrived back to the PC.
        """
        return self.__busy

    @busy.setter
    def busy(self, is_busy):
        """Set the busy state of the camera.

        If the busy state of the camera changes, the camera_busy
        signal will be emitted, carrying the current busy state.

        Emits
        -----
        camera_busy(self.busy)
            If the busy state changes.
        """
        was_busy = self.busy
        is_busy = bool(is_busy)
        if was_busy is not is_busy:
            self.__busy = is_busy
            self.camera_busy.emit(self.busy)

    @property
    def exposure(self):
        """Return the exposure time from the settings.

        Returns
        -------
        exposure_time : float
            Exposure time in milliseconds
        """
        try:
            exposure_time = self.settings.getfloat('measurement_settings',
                                                   'camera_exposure')
        except ValueError:  # Cannot be read as float
            exposure_time = -1

        min_exposure, max_exposure = self.get_exposure_limits()
        if exposure_time < min_exposure or exposure_time > max_exposure:
            emit_error(self, CameraErrors.INVALID_SETTINGS,
                       'measurement_settings/camera_exposure')
        return exposure_time

    @property
    def gain(self):
        """Return the gain from the settings.

        Returns
        -------
        gain : float
            Gain in decibel. If there is no gain setting,
            this method returns 0.
        """
        try:
            gain = self.settings.getfloat('measurement_settings',
                                          'camera_gain', fallback=0)
        except ValueError:  # Cannot be read as float
            gain = -1

        min_gain, max_gain = self.get_gain_limits()
        if gain < min_gain or gain > max_gain:
            emit_error(self, CameraErrors.INVALID_SETTINGS,
                       'measurement_settings/camera_gain')
        return gain

    @property
    def mode(self):
        """Return the camera mode from the settings.

        Returns
        -------
        mode : {'live', 'triggered'}
            The mode from the current settings
        """
        mode = self.settings.get('camera_settings', 'mode',
                                 fallback='triggered')

        if mode not in ('live', 'triggered'):
            emit_error(self, CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                       mode, 'camera_settings/mode', 'triggered')
            mode = 'triggered'
            self.settings.set('camera_settings', 'mode', mode)
        return mode

    @property
    def n_frames(self):
        """Return the number of frames used for averaging."""
        try:
            n_frames = self.settings.getint('measurement_settings',
                                            'n_frames', fallback=1)
        except ValueError:  # Cannot be read as int
            n_frames = -1

        min_n, max_n = self.get_n_frames_limits()
        if n_frames < min_n or n_frames > max_n:
            emit_error(self, CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                       n_frames, 'measurement_settings/n_frames', 1)
            n_frames = 1
            self.settings.set('measurement_settings', 'n_frames', '1')

        return n_frames

    @property
    def name(self):
        """Return the device name for this camera.

        Can only be set by changing
        self.settings['camera_settings']['device_name']
        """
        return self.settings.get('camera_settings', 'device_name')

    @property
    def roi(self):
        """Return the settings of the roi.

        Returns
        -------
        tuple
            Returns an empty tuple if the value in the settings is
            'None' or if any error occurred while trying to interpret
            the setting.

            Otherwise the following elements are returned:
            roi_x, roi_y : int
                Coordinates of the top-left pixel. Zero is the
                topmost/leftmost pixel.
            roi_width, roi_height : int
                Width and height of the region of interest in pixels
        """
        try:
            roi = ast.literal_eval(
                self.settings.get('camera_settings', 'roi', fallback='None')
                )
        except (SyntaxError, ValueError):  # Invalid ROI string
            roi = tuple()
        if roi:
            try:
                roi = tuple(int(r) for r in roi)
            except (TypeError, ValueError):  # Not iterable or not ints
                roi = tuple()

        if not self.__is_valid_roi(roi):
            emit_error(self, CameraErrors.INVALID_SETTINGS,
                       'camera_settings/roi')
            self.settings.set('camera_settings', 'roi', 'None')
            return tuple()

        # See if the roi width and height fit with the binning factor
        *roi_offsets, roi_w, roi_h = roi
        if roi_w % self.binning or roi_h % self.binning:
            new_roi_w = (roi_w // self.binning)*self.binning
            new_roi_h = (roi_h // self.binning)*self.binning
            emit_error(self, CameraErrors.BINNING_ROI_MISMATCH,
                       roi_w, roi_h, self.binning, new_roi_w, new_roi_h)
            roi = (*roi_offsets, new_roi_w, new_roi_h)
            self.settings.set('camera_settings', 'roi', str(roi))
        return roi

    # .settings getter
    def __get_settings(self):
        """Return the settings used for the camera."""
        return self.__settings

    # .settings setter
    def set_settings(self, new_settings):
        """Set new settings for the camera.

        The new settings are accepted only if valid. Otherwise, the
        previous settings are kept. If valid, the settings are also
        loaded to the camera hardware.

        Parameters
        ----------
        new_settings : dict or ConfigParser
            The new settings. The following sections/options
            are mandatory:
            'camera_settings'/'class_name'
                Name of the concrete class that implements
                the abstract methods of CameraABC.
            'camera_settings'/'device_name'
                Name of the device as it appears in the device
                list (e.g., "DMK 33GX265").
            'measurement_settings'/'camera_exposure'
                Exposure time in milliseconds to be used for
                acquiring each frame.
        """
        if new_settings is None:
            self.error_occurred.emit(CameraErrors.MISSING_SETTINGS)
            return

        # Checking of non-mandatory data is done in property getters
        new_settings, invalid = config_has_sections_and_options(
            self, new_settings,
            self._mandatory_settings
            )
        for setting in invalid:
            emit_error(self, CameraErrors.INVALID_SETTINGS, setting)
        if invalid:
            return

        self.__settings = new_settings
        self.stop()
        self.connect()  # Also loads self.settings to camera

    settings = property(__get_settings, set_settings)

    @property
    @abstractmethod
    def supports_trigger_burst(self):
        """Return whether the camera allows triggering multiple frames.
        
        This property should be reimplemented in concrete subclasses.
        Should the camera support trigger burst, get_n_frames_limits
        should also be reimplemented.
        
        Returns
        -------
        supported : bool
            True if the camera internally supports triggering multiple
            frames, i.e., only one 'trigger_now()' call is necessary to
            deliver all the frames needed.
        """
        return False

    def check_loaded_settings(self):
        """Check that camera and configuration settings are the same.

        Returns
        -------
        setting_match : bool
            True if the settings in the camera and in self.settings
            are the same.

        Emits
        -----
        error_occurred(CameraErrors.SETTINGS_MISMATCH)
            If any of the settings in the camera does not match
            the one in self.settings.
        """
        to_check = {k: (getattr(self, k), getattr(self, f"get_{k}"))
                    for k in (*self._abstract,
                              *self.hardware_supported_features)}

        error_code, error_txt = CameraErrors.SETTINGS_MISMATCH
        error_msg = []
        for setting, (in_config, getter) in to_check.items():
            in_device = getter()
            if in_config != in_device:
                error_msg.append(
                    error_txt.format(setting, in_device, in_config)
                    )

        if error_msg:
            self.error_occurred.emit((error_code, '\n'.join(error_msg)))
            return False
        return True

    @abstractmethod
    def close(self):
        """Close the camera device.

        This method is guaranteed to be called exactly once every
        time .disconnect() runs.

        Returns
        -------
        None
        """
        return

    def connect(self):
        """Connect to the camera."""
        if not self.open():
            emit_error(self, CameraErrors.CAMERA_NOT_FOUND, self.name)
            return
        self.load_camera_settings()

    def disconnect(self):
        """Disconnect the device."""
        self.stop()
        self.close()

    @abstractmethod
    def image_info(self):
        """Return information about the last image.

        Returns
        -------
        width, height : int
            Width and height of the image in pixels
        n_bytes : int
            Number of bytes per pixel and per color
        n_colors : int
            Number of color channels
        """
        return

    @abstractmethod
    def list_devices(self):
        """Return a list of available device names.

        Returns
        -------
        devices : list
            Each element is a string representing
            the name of a camera device.
        """
        return

    def load_camera_settings(self):
        """Load settings from self.settings into the camera.

        This method should only be reimplemented if the order in which
        setter-method calls are performed is inappropriate. The default
        order of the calls is:
            set_exposure(),  set_gain(),  set_mode()
            <setters associated to self.hardware_supported_features>
            [set_binning()], [set_n_frames()], [set_roi()]
        The last 3 are called only if they are not hardware-supported
        Should the method be reimplemented, the reimplementation should
        call self.check_loaded_settings() at the end, to verify that
        the settings have been correctly loaded.

        Returns
        -------
        None.
        """
        setters = dict.fromkeys((*self._abstract,
                                 *self.hardware_supported_features,
                                 'binning', 'n_frames', 'roi'))
        for key in setters:
            setter = getattr(self, f"set_{key}")
            setter()
        self.check_loaded_settings()

    @abstractmethod
    def open(self):
        """Open the camera device.

        This method is guaranteed to be called exactly once every time
        .connect() is called. If the device is already open, the method
        should return True, and typically do nothing.

        The reimplementation can use self.name to access the device
        name for this camera.

        Returns
        -------
        successful : bool
            True if the device was opened successfully.
        """
        return False

    @abstractmethod
    def get_binning(self):
        """Return the binning factor internally used for averaging.

        This method should be reimplemented to query the camera for the
        current number of pixels used internally for binning. Binning
        is done on a (binning_factor x binning_factor) mesh. Make sure
        to self.open() before querying the camera.

        This method must return 0 if the camera does not support
        binning internally, i.e., at the firmware level.

        Returns
        -------
        binning_factor : int
            Linear number of pixels used for binning.
        """
        return 0

    def set_binning(self, no_binning=True):
        """Define binning factor for images.

        Do not reimplement this method, unless the camera internally
        supports binning.  Make sure to self.open() before setting.

        Parameters
        ----------
        no_binning : bool
            True if binning should be deactivated, independently of
            the settings. If False, self.binning must be used as the
            binning factor.

        Returns
        -------
        None.
        """
        if self.get_binning():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports binning, "
                "but self.set_binning() was not reimplemented."
                )
        binning_factor = 1 if no_binning else self.binning
        self.process_info.binning = binning_factor

    def get_binning_limits(self):
        """Return the minimum and maximum value for binning.

        This method should be reimplemented only if the camera
        internally supports binning. Make sure to self.open()
        before querying the camera.

        The base implementation uses information on the largest
        available region of interest (i.e., the sensor size) to
        determine the largest binning factor possible (i.e., the
        one giving a 1-pixel-wide (or -high) binned image).

        Returns
        -------
        min_binning, max_binning : int
            The smallest and largest supported binning factors
        """
        if self.get_binning():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports binning, "
                "but self.get_binning_limits() was not reimplemented."
                )
        min_binning = 1
        _, (max_roi_w, max_roi_h) = self.get_roi_size_limits()
        max_binning = min(max_roi_w, max_roi_h)
        return min_binning, max_binning

    @abstractmethod
    def get_exposure(self):
        """Get the exposure time from the camera device.

        Make sure to self.open() before querying the camera.

        Returns
        -------
        exposure_time : float
            Exposure time in milliseconds.
        """
        return

    @abstractmethod
    def set_exposure(self):
        """Set the exposure time for one frame.

        This method must be reimplemented to call the internal driver
        function that sets the appropriate exposure time in the camera
        device. Make sure to self.open() before setting.

        The reimplementation must use the exposure time (in
        milliseconds) returned by self.exposure (i.e., the one
        in self.settings).

        Returns
        -------
        None
        """
        return

    @abstractmethod
    def get_exposure_limits(self):
        """Return the the minimum and maximum exposure times supported.

        Make sure to self.open() before querying the camera.

        Returns
        -------
        min_exposure, max_exposure : float
            Shortest and longest exposure times in milliseconds
        """
        return

    @abstractmethod
    def get_gain(self):
        """Get the gain from the camera device.

        Make sure to self.open() before querying the camera.

        Returns
        -------
        gain : float
            Gain in decibel.
        """
        return

    @abstractmethod
    def set_gain(self):
        """Set the gain of the camera in decibel.

        This method must be reimplemented to call the internal driver
        function that sets the appropriate gain in the camera device.
        Make sure to self.open() before setting.

        The reimplementation must use the gain factor (in decibel)
        returned by self.gain (i.e., the one in self.settings).
        """
        return

    @abstractmethod
    def get_gain_limits(self):
        """Return the the minimum and maximum gains supported.

        Make sure to self.open() before querying the camera.

        Returns
        -------
        min_gain, max_gain : float
            Smallest and largest gain factors in decibel
        """
        return

    @abstractmethod
    def get_mode(self):
        """Return the mode set in the camera.

        This method should be reimplemented to query the camera for the
        current mode, and return the appropriate string. Make sure to
        self.open() before querying the camera.

        Should the camera not support natively a live vs. triggered
        mode, this function can return self.mode.

        Returns
        -------
        mode : {'live', 'triggered'}
            The mode the camera is operating in. Continuous mode
            corresponds to receiving a 'stream' of images at
            well-defined time intervals (i.e., synchronous), while
            'triggered' is asynchronous: the camera returns a frame
            only when asked by self.trigger_now().
        """
        return

    @abstractmethod
    def set_mode(self):
        """Set the camera mode from self.settings.

        This method should be reimplemented to instruct the camera
        to switch to a desired mode, which must be retrieved with
        self.mode from the settings. Make sure to self.open() before
        setting.

        Should the camera not support natively a live vs. triggered
        mode, live mode can be simulated with a QTimer.
        """
        # TODO: we may want to not make this an abstract method,
        # and have get_mode return None or '' if the camera does
        # not support live-view. Then the base implementation would
        # prepare the system for a live-view simulation using QTimer
        return

    @abstractmethod
    def get_n_frames(self):
        """Return the number of frames internally used for averaging.

        This method should be reimplemented to query the camera
        for the current number of frames internally used for
        averaging. Make sure to self.open() before querying.

        Returns
        -------
        n_frames : int
            Number of frames used for averaging. Returns 0 if the
            camera does not support internally frame averaging.
        """
        return 0

    def set_n_frames(self):
        """Set the number of frames internally used for averaging.

        This method should be reimplemented to instruct the camera
        to use self.n_frames frames for averaging internally, before
        returning the image. Make sure to self.open() before setting.

        Should the camera not natively support frame averaging, this
        method should not be reimplemented, as frame averaging will
        be done via software during image processing.
        """
        if self.get_n_frames():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports frame "
                "averaging, but self.set_n_frames() was not reimplemented"
                )
        self.process_info.n_frames = self.n_frames

    def get_n_frames_limits(self):
        """Return the minimum and maximum number of frames supported.

        This method should be reimplemented only if the camera
        internally supports frame averaging. Make sure to .open()
        before querying the camera.

        Returns
        -------
        min_n_frames, max_n_frames : int
            Minimum and maximum number of frames usable for averaging
        """
        if self.get_n_frames():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports frame averaging,"
                " but self.get_n_frames_limits() was not reimplemented"
                )
        if self.supports_trigger_burst:
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports trigger burst, "
                "but self.get_n_frames_limits() was not reimplemented"
                )
        return 1, np.inf

    @abstractmethod
    def get_roi(self):
        """Return the region of interest set in the camera.

        This method should be reimplemented in concrete subclasses to
        query the camera for the current settings of the region of
        interest. Make sure to self.open() before querying the camera

        If the camera does not support setting a region of
        interest this method should return an empty tuple.

        Returns
        -------
        tuple
            The settings of the region of interest. If the camera
            supports hardware ROI, this method returns:
            roi_x, roi_y : int
                Coordinates of the top-left pixel. Zero is the
                topmost/leftmost pixel.
            roi_width, roi_height : int
                Width and height of the region of interest in pixels
        """
        return tuple()

    def set_roi(self, no_roi=True):
        """Set up region of interest in the camera.

        This method should be reimplemented only if the camera natively
        supports setting a region of interest at the hardware level.
        The reimplementation must retrieve the ROI information using
        self.roi. Make sure to .open() before setting.

        If the camera does not support setting internally a ROI, the
        the ROI setting is used to crop the frames during processing.

        Parameters
        ----------
        no_roi : bool, optional
            If True, set the ROI to the full size of the sensor rather
            than using the value from the settings. Default is True.

        Returns
        -------
        None.
        """
        if self.get_roi():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports setting a "
                "region of interest, but self.set_roi() was not reimplemented"
                )
        roi = tuple() if no_roi else self.roi
        self.process_info.roi = roi

    @abstractmethod
    def get_roi_size_limits(self):
        """Return minimum, maximum and granularity of the ROI.

        This method must be reimplemented in concrete subclasses.
        Typically it would return (1, 1) for the minimum ROI,
        (sensor_width, sensor_height) for the maximum ROI, and
        (1, 1) for the minimum increments. Make sure to .open()
        before querying the camera.

        Returns
        -------
        roi_min : tuple
            Two elements, both integers, corresponding to the
            minimum width and minimum height
        roi_max : tuple
            Two elements, both integers, corresponding to the
            maximum width and maximum height
        roi_increments : tuple
            Two elements, both integers, corresponding to the
            minimum allowed increments for width and height of
            the region of interest
        """
        return tuple(), tuple(), tuple()

    @abstractmethod
    def reset(self):
        """Reset the camera to default mode.

        After calling this method, the camera settings in .settings
        and those effectively 'loaded' in the camera may differ.

        Call self.load_camera_settings() to reload the settings.

        Returns
        -------
        None
        """
        return

    @abstractmethod
    def start(self):
        """Start the camera.

        The camera is started in the mode returned by self.mode.

        This method should be extended in concrete subclasses, i.e.,
        super().start() must be called in the reimplementation.

        Returns
        -------
        None.
        """
        if self.mode == 'triggered':
            self.__process_thread.start()

    @abstractmethod
    def stop(self):
        """Stop the camera.

        This method should be extended in concrete subclasses, i.e.,
        super().stop() must be called in the reimplementation.

        Returns
        -------
        None.
        """
        if self.__process_thread.isRunning():
            self.__process_thread.quit()

    @abstractmethod
    def trigger_now(self):
        """Start acquiring one (or more) frames now.

        This method should be reimplemented in concrete subclasses
        to signal the camera that acquisition of one (or multiple)
        frames should start now.

        The reimplementation must call super().trigger_now() before,
        as this will emit an error_occurred if the camera mode is
        inappropriate (i.e., the camera is not in triggered mode).
        It will set the camera to busy and set n_frames_done to 0
        as well.

        Returns
        -------
        None

        Emits
        -----
        error_occurred(CameraErrors.UNSUPPORTED_OPERATION)
        """
        if self.mode != 'triggered':
            emit_error(self, CameraErrors.UNSUPPORTED_OPERATION,
                       'trigger', 'live')
        self.busy = True
        self.n_frames_done = 0

    def __is_valid_roi(self, roi):
        """Check that ROI is OK and fits with limits.

        Parameters
        ----------
        roi : tuple
            The region of interest to be checked

        Returns
        -------
        roi_OK : bool
            True if roi is a valid ROI and if it fits in the limits.
        """
        if not roi:
            return False
        if len(roi) != 4:
            return False
        roi_x, roi_y, roi_w, roi_h = roi
        (min_w, min_h), (max_w, max_h), (d_w, d_h) = self.get_roi_size_limits()
        if (roi_x < 0 or roi_y < 0 or roi_w < min_w or roi_h < min_h
            or roi_x + roi_w > max_w or roi_y + roi_h > max_h
                or roi_w % d_w or roi_h % d_h):
            return False
        return True

    def __on_image_saved(self):
        """React to an image being saved."""
        # Disconnect the last processor to prevent
        # duplicate processing of the next frame
        processor = qtc.QObject().sender()
        if processor:
            try:
                self.__process_frame.disconnect(processor.process_frame)
            except TypeError:
                # Not connected
                pass

    def __on_frame_ready(self, image):
        """React to receiving a new frame from the camera.

        This is the slot connected to the frame_ready signal.

        Parameters
        ----------
        image : numpy.ndarrray
            The new frame as an array

        Returns
        -------
        None
        """
        if self.mode != 'triggered':
            # In live mode, the frames will be handled by the GUI
            return

        try:
            self.process_info.bad_pixels = self.settings.get(
                'camera_settings', 'bad_pixels', fallback=tuple()
                )
        except (TypeError, ValueError):
            emit_error(self, CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                       'camera_settings/bad_pixels', 'no bad pixels')
        # TODO: file name. Probably needs a temp path + a file name
        #       both should be taken care of by the measurement.
        # TODO: do I need to keep a reference here??

        processor = ImageProcessor()
        processor.image_saved.connect(self.__on_image_saved)
        processor.image_saved.connect(processor.deleteLater)
        self.__process_thread.finished.connec(processor.deleteLater)
        self.__process_frame.connect(processor.process_frame)
        processor.prepare_to_process(self.process_info.copy(), image)
        processor.moveToThread(self.__process_thread)                     # TODO: may need to be moved before connecting
        # processor.moveToThread(self.__process_thread.currentThread())   # TODO: may be the way to go as the thread is running already

        self.__process_frame.emit(image)
        self.n_frames_done += 1

        if self.n_frames_done < self.n_frames:
            # More frames to be acquired
            if self.supports_trigger_burst:
                return
            self.trigger_now()
        else:
            # All frames are done
            self.busy = False

    @abstractmethod
    def set_callback(self, on_frame_ready):
        """Pass a frame-ready callback to the camera driver.

        If the camera does not support having a callback function,
        a similar behavior can be obtained using an appropriate
        pyqtSignal, emitted as soon as a frame has been acquired.

        Parameters
        ----------
        on_frame_ready : callable
            The function that will be called by the camera each time
            a new frame arrives. The callback should only care of
            converting the data from the camera into a numpy.array
            of appropriate shape and data type, then emitting a
            frame_ready signal carrying the array. It must be able
            to take a reference to self as part of its arguments:
            It may do so either taking self directly or taking
            self.process_info (.camera is a reference to self).
            It can then access methods of the driver via self.driver.

        Returns
        -------
        None
        """
        return
