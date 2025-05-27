"""Module abc of viperleed.gui.measure.camera.

Defines the abstract functionality for all cameras that should
be supported by ViPErLEED. This includes error codes/messages
relative to the cameras (CameraErrors class) and the abstract
base-class CameraABC. All classes handling cameras should be
concrete subclasses of CameraABC.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-12'
__license__ = 'GPLv3+'

from abc import abstractmethod

import numpy as np
from PyQt5 import QtCore as qtc

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.camera.imageprocess import ImageProcessInfo
from viperleed.gui.measure.camera.imageprocess import ImageProcessor
from viperleed.gui.measure.camera import badpixels
from viperleed.gui.measure.classes.settings import NoSettingsError
from viperleed.gui.measure.classes.settings import NotASequenceError
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings


# pylint: disable=too-many-lines,too-many-public-methods
# Disabled too-many-lines because many are documentation
# of the abstract methods, and too-many-public-methods
# because they are only meant to be reimplemented in
# concrete subclasses and used internally, not called
# from user code.


class CameraErrors(base.ViPErLEEDErrorEnum):
    """Data class for base camera errors."""

    INVALID_SETTINGS = (200,
                        "Invalid camera settings: Required settings "
                        "{!r} missing or values inappropriate. "
                        "Check configuration file.\n{}")
    MISSING_SETTINGS = (201,
                        "Camera cannot operate without settings. Load "
                        "an appropriate settings file before proceeding.")
    CAMERA_NOT_FOUND = (202,
                        "Could not find camera {}.\nMake sure the camera is "
                        "connected and has power. Try (re-)plugging it, and "
                        "give the camera enough time to boot up.")
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
                            " with binning factor {}. Reducing region of "
                            "interest to ({} x {}). A few pixels on the "
                            "lower-right corner may be removed.")
    UNSUPPORTED_WHILE_BUSY = (207, "Cannot {} while camera is busy.")
    TIMEOUT = (208,   # Only in triggered mode
               "No frames returned by camera {} in the last {} seconds. "
               "Check that the camera is plugged in and powered. If it "
               "is, try rebooting the camera.")


class CameraABC(qtc.QObject, metaclass=base.QMetaABC):
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

    # image_processed is emitted when operating in triggered mode
    # after all post-processing steps have been performed.
    image_processed = qtc.pyqtSignal(np.ndarray)

    # image_saved is emitted after an image has been saved to
    # disk, carrying the path to the file saved.
    image_saved = qtc.pyqtSignal(str)

    # started/stopped are emitted whenever the camera is
    # started/stopped. When started, self.mode can be used
    # to deduce the camera mode.
    started = qtc.pyqtSignal()
    stopped = qtc.pyqtSignal()

    # preparing is emitted every time a camera begins (True)
    # or finishes (False) any pre-acquisition tasks.
    preparing = qtc.pyqtSignal(bool)

    # __process_frame is a private signal used to pass frames to
    # the processing thread
    __process_frame = qtc.pyqtSignal(np.ndarray)

    _mandatory_settings = [
            ('camera_settings', 'class_name'),
            ('camera_settings', 'device_name'),
            ('measurement_settings', 'exposure'),
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

    def __init__(self, driver, *__args,
                 settings=None, parent=None, **__kwargs):
        """Initialize instance.

        Parameters
        ----------
        driver : object
            An instance of the library/class that handles
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
        parent : QObject
            The parent QObject of this Camera.
        **kwargs : object
            Other unused keyword arguments.
        """
        super().__init__(parent=parent)
        self.driver = driver
        self.__busy = False
        self.__connected = False
        self.__settings = ViPErLEEDSettings()
        self.__bad_pixels = None
        self.__timeout = qtc.QTimer(parent=self)
        self.__timeout.setSingleShot(True)
        self.__timeout.timeout.connect(self.__on_timed_out)

        # Store some property values to limit the number
        # of potentially slow checks. Values are set in
        # the property getters where appropriate
        self.__properties = {'binning': 1, 'exposure': 1.0, 'gain': 0.0,
                             'n_frames': 1, 'roi': tuple()}

        # Keep track of some errors that we want to report
        # only once. This 'list' is cleared every time a
        # camera is disconnected.
        self.__reported_errors = set()

        self.process_info = ImageProcessInfo()
        self.process_info.camera = self
        self.__image_processors = []
        self.__process_thread = qtc.QThread()
        self.__retry_stop_timer = qtc.QTimer(parent=self)
        self.__retry_stop_timer.setSingleShot(True)
        self.__retry_stop_timer.timeout.connect(self.__stop_now_or_retry)

        self.__init_errors = []  # Report these with a little delay
        self.__init_err_timer = qtc.QTimer(parent=self)
        self.__init_err_timer.setSingleShot(True)

        self.error_occurred.connect(self.__on_init_errors)
        self.__init_err_timer.timeout.connect(self.__report_init_errors)

        # Number of frames accumulated for averaging
        self.n_frames_done = 0

        try:
            self.set_settings(settings)
        except self.exceptions:
            # Start a short QTimer to report errors that occurred here
            # AFTER the whole __init__ is done, i.e., when we suppose
            # that the error_occurred signal is already connected
            # externally.
            pass
        if self.__init_errors:
            self.__init_err_timer.start(20)
        self.error_occurred.disconnect(self.__on_init_errors)

    def __deepcopy__(self, memo):
        """Return self rather than a deep copy of self."""
        return self

    @property
    def bad_pixels(self):
        """Return a BadPixels object for this camera."""
        # If __bad_pixels is False-y (i.e., info never read),
        # attempt reading the bad pixels info from file,
        # based on the path in settings.
        if not self.__bad_pixels:
            if self.__bad_pixels is None:
                self.__bad_pixels = badpixels.BadPixels(self)
            self.update_bad_pixels()
        return self.__bad_pixels

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
        except (TypeError, ValueError):  # Cannot be read as int
            binning_factor = -1

        _stored = self.__properties['binning']
        if _stored == binning_factor:
            return binning_factor

        min_bin, max_bin = self.get_binning_limits()
        if binning_factor < min_bin or binning_factor > max_bin:
            base.emit_error(self, CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                            f"{binning_factor} "
                            f"[out of range ({min_bin}, {max_bin})]",
                            'camera_settings/binning', 1)
            binning_factor = 1
            self.settings.set('camera_settings', 'binning', '1')
        self.__properties['binning'] = binning_factor
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
    def connected(self):
        """Return whether the camera is connected.

        Notice that a camera may be connected [i.e., .connect_()
        was called] but not running [i.e., .start() was not].

        Returns
        -------
        connected : bool
            True if the camera is connected.
        """
        return self.__connected

    @property
    @abstractmethod
    def exceptions(self):
        """Return a tuple of camera exceptions.

        Returns
        -------
        exceptions : tuple
            Each element is a Exception subclass of exceptions
            that the camera may raise in case internal driver
            errors occur.
        """
        return tuple()

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
                                                   'exposure')
        except (TypeError, ValueError):  # Cannot be read as float
            # pylint: disable=redefined-variable-type
            # Seems a bug: getfloat always returns a float
            exposure_time = -1

        _stored = self.__properties['exposure']
        if _stored == exposure_time:
            return exposure_time

        min_exp, max_exp = self.get_exposure_limits()
        if exposure_time < min_exp or exposure_time > max_exp:
            base.emit_error(self, CameraErrors.INVALID_SETTINGS,
                            'measurement_settings/exposure',
                            f'Info: out of range ({min_exp}, {max_exp})')
            exposure_time = min_exp
        self.__properties['exposure'] = exposure_time
        return exposure_time

    @property
    def frame_interval(self):
        """Return the time interval (msec) between frames."""
        return max(self.exposure, 1000/self.get_frame_rate())

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
                                          'gain', fallback=0)
        except (TypeError, ValueError):  # Cannot be read as float
            # pylint: disable=redefined-variable-type
            # Seems a bug: getfloat always returns a float
            gain = -1

        _stored = self.__properties['gain']
        if _stored == gain:
            return gain

        min_gain, max_gain = self.get_gain_limits()
        if gain < min_gain or gain > max_gain:
            base.emit_error(self, CameraErrors.INVALID_SETTINGS,
                            'measurement_settings/gain',
                            f'Info: out of range ({min_gain}, {max_gain})')
            gain = min_gain
            self.settings.set('measurement_settings', 'gain', str(gain))
        self.__properties['gain'] = gain
        return gain

    @property
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

    @property
    @abstractmethod
    def intensity_limits(self):
        """Return the minimum and maximum value for a pixel.

        Returns
        -------
        pixel_min : int
            Minimum intensity for a pixel
        pixel_max : int
            Maximum intensity for a pixel
        """
        return

    @property
    @abstractmethod
    def is_running(self):
        """Return whether the camera is currently running."""
        return False

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
            base.emit_error(self, CameraErrors.INVALID_SETTING_WITH_FALLBACK,
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
        except (TypeError, ValueError):  # Cannot be read as int
            # pylint: disable=redefined-variable-type
            # Seems a bug: getint always returns an integer
            n_frames = -1

        _stored = self.__properties['n_frames']
        if _stored == n_frames:
            return n_frames

        min_n, max_n = self.get_n_frames_limits()
        if n_frames < min_n or n_frames > max_n:
            base.emit_error(self, CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                            f"{n_frames} [out of range ({min_n}, {max_n})]",
                            'measurement_settings/n_frames', 1)
            n_frames = 1
            self.settings.set('measurement_settings', 'n_frames', '1')
        self.__properties['n_frames'] = n_frames
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
            roi = self.settings.getsequence('camera_settings', 'roi',
                                            fallback=tuple())
        except NotASequenceError:  # Invalid ROI string
            roi = tuple()

        try:
            roi = tuple(int(r) for r in roi)
        except (TypeError, ValueError):  # Not iterable or not ints
            roi = tuple()

        _stored = self.__properties['roi']
        if _stored and _stored == roi:
            return roi

        _, (max_w, max_h), _, (roi_dx, roi_dy) = self.get_roi_size_limits()

        if not self.__is_valid_roi(roi):
            base.emit_error(self, CameraErrors.INVALID_SETTINGS,
                            'camera_settings/roi', 'Info: ROI is invalid.')
            self.settings.set('camera_settings', 'roi', 'None')
            return 0, 0, max_w, max_h

        # See if the ROI position offsets fit the
        # increments, and adjust if necessary
        *roi_offsets, roi_w, roi_h = roi
        roi_x, roi_y = roi_offsets
        if roi_x % roi_dx or roi_y % roi_dy:
            roi_offsets = ((roi_x // roi_dx) * roi_dx,
                           (roi_y // roi_dy) * roi_dy)
            new_roi = (*roi_offsets, roi_w, roi_h)
            if not self.__is_valid_roi(roi):
                base.emit_error(self, CameraErrors.INVALID_SETTINGS,
                                'camera_settings/roi',
                                f'Info: ROI {new_roi} is invalid after '
                                'adjusting top-left corner position')
                self.settings.set('camera_settings', 'roi', 'None')
                return 0, 0, max_w, max_h
            self.settings.set('camera_settings', 'roi', str(new_roi))

        # See if the ROI width and height fit with the binning factor
        if roi_w % self.binning or roi_h % self.binning:
            new_roi_w = (roi_w // self.binning) * self.binning
            new_roi_h = (roi_h // self.binning) * self.binning
            error = (CameraErrors.BINNING_ROI_MISMATCH,
                     roi_w, roi_h, self.binning, new_roi_w, new_roi_h)
            if not self.__already_reported(error):
                # base.emit_error(self, *error)                                 # TODO: make this a non-critical warning!
                self.__reported_errors.add(error)
            # Update the ROI in the settings only if the
            # new ROI is a valid one for the camera. The
            # image processor will anyway crop off the
            # extra pixels on the lower-right corner of
            # the image to apply binning.
            new_roi = (*roi_offsets, new_roi_w, new_roi_h)
            if self.__is_valid_roi(new_roi):
                self.settings.set('camera_settings', 'roi', str(new_roi))
                roi = new_roi
        self.__properties['roi'] = roi
        return roi

    @property
    def sensor_size(self):
        """Return the full width and height of the sensor.

        The base implementation relies on the results returned
        by .get_roi_size_limits(). This property should be
        reimplemented if the sensor size is larger than the
        largest ROI applicable.

        Returns
        -------
        width, height : int
            Width and height of the sensor in pixels.
        """
        return self.get_roi_size_limits()[1]

    @property
    def settings(self):
        """Return the settings used for the camera."""
        return self.__settings

    @settings.setter
    def settings(self, new_settings):
        """Set new settings for the camera."""
        self.set_settings(new_settings)

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
        try:
            new_settings = ViPErLEEDSettings.from_settings(new_settings)
        except (ValueError, NoSettingsError):
            base.emit_error(self, CameraErrors.MISSING_SETTINGS)
            return

        # Checking of non-mandatory data is done in property getters
        invalid = new_settings.has_settings(*self._mandatory_settings)
        if invalid:
            base.emit_error(self, CameraErrors.INVALID_SETTINGS,
                            ', '.join(invalid), '')
            return

        # Keep track of the old ROI: will recalculate
        # bad pixel coordinates if ROI changes
        old_roi = tuple()
        if self.settings:
            old_roi = self.roi

        self.__settings = new_settings
        if self.is_running:
            self.stop()
        self.close()
        self.connect_()  # Also loads self.settings to camera

        if self.roi != old_roi and self.bad_pixels:
            self.bad_pixels.apply_roi()

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
        time .disconnect_() runs.

        Returns
        -------
        None
        """
        return

    def connect_(self):
        """Connect to the camera."""
        if not self.open():
            base.emit_error(self, CameraErrors.CAMERA_NOT_FOUND, self.name)
            self.__connected = False
            return
        self.load_camera_settings()
        self.__connected = True

    def disconnect_(self):
        """Disconnect the device."""
        self.__reported_errors = set()
        self.stop()
        self.close()
        self.__connected = False

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
        The last 3 are called only if they are not hardware-supported.
        Should the method be reimplemented, the reimplementation should
        call self.check_loaded_settings() at the end, to verify that
        the settings have been correctly loaded.

        Returns
        -------
        None.
        """
        self.busy = True
        self.preparing.emit(self.busy)
        extra = (k for k in ('binning', 'n_frames', 'roi')
                 if k not in self.hardware_supported_features)
        setters = dict.fromkeys((*self._abstract,
                                 *self.hardware_supported_features,
                                 *extra))
        for key in setters:
            setter = getattr(self, f"set_{key}")
            setter()
        self.check_loaded_settings()
        self.busy = False
        self.preparing.emit(self.busy)

    @abstractmethod
    def open(self):
        """Open the camera device.

        This method is guaranteed to be called exactly once every
        time .connect_() is called. If the device is already open,
        the method should return True, and typically do nothing.

        The reimplementation can use self.name to access the device
        name for this camera.

        Returns
        -------
        successful : bool
            True if the device was opened successfully.
        """
        return False

    def update_bad_pixels(self):
        """Update the list of bad pixels from the settings.

        This function should be explicitly called after the
        path containing bad pixels files is changed in the
        settings. Calling this function may take a little
        time, as bad pixel coordinates are read from scratch
        and recalculated.

        Emits
        -----
        error_occurred(CameraErrors.INVALID_SETTINGS)
            If settings file does not have a 'bad_pixels_path' option
            in section 'camera_settings'.
        error_occurred(CameraErrors.INVALID_SETTINGS)
            If the path specified does not contain a valid bad-pixels
            file for this camera.
        """
        bad_pix_path = self.settings.get("camera_settings", "bad_pixels_path",
                                         fallback='')
        if not bad_pix_path:
            base.emit_error(self, CameraErrors.INVALID_SETTINGS,
                            'camera_settings/bad_pixels_path',
                            'Info: No bad_pixel_path found.')
            return
        try:
            self.__bad_pixels.read(bad_pix_path)
        except (FileNotFoundError, ValueError) as err:
            base.emit_error(self, CameraErrors.INVALID_SETTINGS,
                            'camera_settings/bad_pixels_path', f'Info: {err}')
            return
        self.__bad_pixels.apply_roi()

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

    def set_binning(self, no_binning=False):
        """Define binning factor for images.

        Do not reimplement this method, unless the camera internally
        supports binning.  Make sure to self.open() before setting.

        Parameters
        ----------
        no_binning : bool, optional
            True if binning should be deactivated, independently of
            the settings. If False, self.binning must be used as the
            binning factor. Default is False.

        Raises
        ------
        NotImplementedError
            If this method is not reimplemented for a camera
            that natively supports binning.
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

        Raises
        ------
        NotImplementedError
            If this method is not reimplemented for a camera
            that natively supports binning.
        """
        if self.get_binning():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports binning, "
                "but self.get_binning_limits() was not reimplemented."
                )
        min_binning = 1
        _, (max_roi_w, max_roi_h), *_ = self.get_roi_size_limits()
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
    def get_frame_rate(self):
        """Return the number of frames delivered per second.

        This property should be reimplemented in concrete subclasses.

        Notice that this quantity is unrelated to 1000/self.exposure.
        If self.exposure/1000 is larger than 1/frame_rate, frames are
        returned at self.exposure time intervals. Otherwise at intervals
        of 1/frame_rate duration. You can use self.frame_interval to
        get the correct time interval between frames.

        Returns
        -------
        frame_rate : float
            Number of frames delivered per second.
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

        Raises
        ------
        NotImplementedError
            If this method is not reimplemented for cameras that
            natively support averaging over multiple frames, or
            for those that natively support a trigger-burst mode
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

    def set_roi(self, no_roi=False):
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
            than using the value from the settings. Default is False.

        Raises
        -------
        NotImplementedError
            If this method is not reimplemented when the camera
            natively supports applying a region of interest.
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
        roi_offset_increments : tuple
            Two elements, both integers, corresponding to the
            minimum allowed increments for the horizontal and
            vertical position of the roi.
        """
        return tuple(), tuple(), tuple(), tuple()

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

    @abstractmethod
    def start(self, *_):
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
            self.n_frames_done = 0

        # Connect the frame_ready signal if not connected already. Should
        # be done here rather than in __init__ as the camera may need to
        # be internally started (via .driver) to perform some of the
        # other pre-starting operations, and may deliver frames in the
        # meantime. The signal is also disconnected in .stop() for the
        # same reason.
        base.safe_connect(self.frame_ready, self.__on_frame_ready,
                          type=qtc.Qt.UniqueConnection)

        # Warn if exposure settings are somewhat stupid
        frame_interval = self.frame_interval
        if self.exposure < frame_interval:
            print(                                                              # TODO: should become a non-fatal warning
                f"WARNING: Exposure ({self.exposure} ms) of camera "
                f"{self.name} is shorter than the time it takes to "
                f"deliver frames ({frame_interval:.2f} ms). Increase "
                "the exposure time to avoid wasting time."
                )
        self.started.emit()

    @abstractmethod
    def stop(self):
        """Stop the camera.

        The camera may not be immediately stopped if there is any
        image processing currently running. In this case, stopping
        will be attempted multiple times with a 100 ms interval.
        Once the camera has effectively stopped, the .stopped()
        signal is emitted.

        This method should be extended in concrete subclasses, i.e.,
        super().stop() must be called in the reimplementation.

        Returns
        -------
        None.
        """
        self.__stop_now_or_retry()

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
            base.emit_error(self, CameraErrors.UNSUPPORTED_OPERATION,
                            'trigger', 'live')
            return
        self.busy = True

        if self.supports_trigger_burst:
            self.process_info.clear_times()

        # Start the timeout timer: will fire if it takes
        # more than 5 s longer than the expected time to
        # receive all frames needed for averaging.
        # The correct time interval would be
        #         exposure + 1/frame_rate
        #        + (n_frames - 1) * self.frame_interval
        # with .frame_interval == max(exposure, 1/frame_rate).
        # However, we don't need to be super precise here,
        # and we can overestimate the interval to:
        frames_time = (self.n_frames + 1) * self.frame_interval
        self.__timeout.start(int(frames_time + 5000))

    def __already_reported(self, error_details):
        """Return whether an error was already reported.

        Parameters
        ----------
        error_details : tuple
            The error to be checked. Typically the first element
            is a ViPErLEEDErrorEnum. The following entries are
            parameters of the specific error. They should be
            given only in case a specific error with certain
            specific values of the parameters should be checked.
        """
        return error_details in self.__reported_errors

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
        ((min_w, min_h),
         (max_w, max_h),
         (d_w, d_h),
          _) = self.get_roi_size_limits()

        roi_outside = (roi_x < 0 or roi_y < 0
                       or roi_x + roi_w > max_w
                       or roi_y + roi_h > max_h)
        roi_small = roi_w < min_w or roi_h < min_h
        roi_not_multiple = roi_w % d_w or roi_h % d_h
        return not (roi_outside or roi_small or roi_not_multiple)

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

        if self.n_frames_done == 0:  # pylint: disable=compare-to-zero
            processor = ImageProcessor()
            processor.image_processed.connect(self.image_processed)
            processor.image_saved.connect(self.__on_image_saved)
            self.__process_thread.finished.connect(processor.deleteLater)
            self.__process_frame.connect(processor.process_frame)
            processor.prepare_to_process(self.process_info.copy(), image)
            processor.moveToThread(self.__process_thread)

            # Keep a reference to the processor to
            # prevent it from being garbage-collected
            self.__image_processors.append(processor)

        self.__process_frame.emit(image)
        self.n_frames_done += 1

        if self.n_frames_done < self.n_frames:
            # More frames to be acquired
            if self.supports_trigger_burst:
                return
            self.trigger_now()
        else:
            # All frames are done
            self.__timeout.stop()
            self.n_frames_done = 0
            self.busy = False

    @qtc.pyqtSlot(str)
    def __on_image_saved(self, fpath):
        """React to an image being saved."""
        processor = self.sender()
        if processor:
            # Disconnect the sender processor to prevent
            # duplicate processing of the next frame
            base.safe_disconnect(self.__process_frame, processor.process_frame)
            try:
                self.__image_processors.remove(processor)
            except ValueError:
                # Not in list
                pass

        # Emit image_saved only if a file was actually saved
        if fpath:
            self.image_saved.emit(fpath)

    def __on_init_errors(self, err):
        """Collect initialization errors to report later."""
        self.__init_errors.append(err)

    def __on_timed_out(self):
        """Report a timeout error while in triggered mode."""
        base.emit_error(self, CameraErrors.TIMEOUT, self.name, 5)

    def __report_init_errors(self):
        """Emit error_occurred for each initialization error."""
        for error in self.__init_errors:
            self.error_occurred.emit(error)
        self.__init_errors = []

    def __stop_now_or_retry(self):
        """Stop camera if image processing is over."""
        # Delay the stopping as long as there are image
        # processors that have not finished their tasks
        if (self.is_running
            and self.__process_thread.isRunning()
                and any(p.busy for p in self.__image_processors)):
            if not self.__retry_stop_timer.isActive():
                self.__retry_stop_timer.start(100)
            return

        # Now it is safe to clean up and stop
        self.__retry_stop_timer.stop()
        self.__timeout.stop()
        self.process_info.clear_times()
        self.n_frames_done = 0
        self.busy = False
        if self.__process_thread.isRunning():
            self.__process_thread.quit()
        base.safe_disconnect(self.frame_ready, self.__on_frame_ready)
        self.stopped.emit()
