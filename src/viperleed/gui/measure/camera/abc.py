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
from collections.abc import Sequence
from functools import wraps
import operator

import numpy as np
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.camera.imageprocess import ImageProcessInfo
from viperleed.gui.measure.camera.imageprocess import ImageProcessor
from viperleed.gui.measure.camera import badpixels
from viperleed.gui.measure.classes import settings as _m_settings
from viperleed.gui.measure.classes.abc import DeviceABC
from viperleed.gui.measure.classes.abc import DeviceABCErrors
from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.dialogs.settingsdialog import SettingsTag
from viperleed.gui.measure.widgets.roieditor import ROIEditor
from viperleed.gui.measure.widgets.pathselector import PathSelector
from viperleed.gui.measure.widgets.spinboxes import CoercingDoubleSpinBox
from viperleed.gui.measure.widgets.spinboxes import CoercingSpinBox


# pylint: disable=too-many-lines,too-many-public-methods
# Disabled too-many-lines because many are documentation
# of the abstract methods, and too-many-public-methods
# because they are only meant to be reimplemented in
# concrete subclasses and used internally, not called
# from user code.


def fallback_if_disconnected(func):
    """Return the parent class implementation if not connected."""
    @wraps(func)
    def _wrapper(self, *args, **kwargs):
        if not self.connected:
            parent_method = getattr(CameraABC, func.__name__)
            return parent_method(self, *args, **kwargs)
        return func(self, *args, **kwargs)
    return _wrapper


class CameraErrors(base.ViPErLEEDErrorEnum):
    """Data class for base camera errors."""

    SETTINGS_MISMATCH = (200,
                         'Different {} settings found in camera and '
                         'configuration file: camera={}, settings={}.')
    UNSUPPORTED_OPERATION = (201, 'Cannot {} in {} mode. Switch mode before.')
    BINNING_ROI_MISMATCH = (202,
                            'Region of interest size ({} x {}) is incompatible'
                            ' with binning factor {}. Reducing region of '
                            'interest to ({} x {}). A few pixels on the '
                            'lower-right corner may be removed.')
    UNSUPPORTED_WHILE_BUSY = (203, 'Cannot {} while camera is busy.')
    TIMEOUT = (204,   # Only in triggered mode
               'No frames returned by camera {} in the last {} seconds. '
               'Check that the camera is plugged in and powered. If it '
               'is, try rebooting the camera.')


class CameraABC(DeviceABC):
    """Abstract base class for ViPErLEED cameras."""

    # frame_ready should be emitted in any callback function that the
    # camera may call. If the camera does not allow to call a callback
    # function, the subclass should take care of emitting a frame_ready
    # signal as soon as a new frame arrives at the PC. Whenever this
    # signal is emitted, the new frame will be processed (in a separate
    # thread).
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

    _mandatory_settings = (
            ('camera_settings', 'class_name'),
            ('camera_settings', 'device_name'),
            ('measurement_settings', 'exposure'),
            )

    # _abstract is a list of features for which setter and getter
    # methods must always be overridden in concrete classes
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
        super().__init__(settings=settings, parent=parent)
        if not driver:
            raise ValueError(f'{driver} is not a valid camera driver.')
        self.driver = driver
        self.__connected = False
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
        self.__retry_stop_timer.timeout.connect(self.stop)

        # Number of frames accumulated for averaging
        self.n_frames_done = 0

        # Calibration tasks. See calibration_tasks property.
        self.__calibration_tasks = {'bad_pixels': [], 'starting': []}

        with self.errors_delayed():
            try:
                self.set_settings(self._settings_to_load)
            except self.exceptions:
                pass

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
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
                f'{binning_factor} '
                f'[out of range ({min_bin}, {max_bin})]',
                'camera_settings/binning', 1
                )
            binning_factor = 1
            self.settings.set('camera_settings', 'binning', '1')
        self.__properties['binning'] = binning_factor
        return binning_factor

    @property
    def calibration_tasks(self):
        """Return a dictionary of CameraCalibrationTask for self.

        If tasks are to be added, this should be done at runtime,
        i.e., by extending this property! This ensures that all
        the information is present in the camera by the time the
        task is created.

        Returns
        -------
        tasks : dict
            Keys can be used to discern when the task is to be
            performed. Values are lists of CameraCalibrationTask
            instances. Tasks are executed in order.
            Typical keys are:
                'bad_pixels'
                    Tasks that should be executed before running
                    the bad-pixels finding CameraCalibrationTask
                'starting'
                    Tasks to execute before completing a call to
                    .start(). The camera will not emit .started
                    until these tasks have been performed.
                    CURRENTLY UNUSED.
        """
        return self.__calibration_tasks

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

    @connected.setter                                                           # TODO: Check if this is used outside of CameraABC. Replace with private set method
    def connected(self, is_connected):
        """Set whether the camera is currently connected. For internal use."""
        was_connected = self.connected
        self.__connected = bool(is_connected)
        if was_connected != self.connected:
            self.connection_changed.emit(self.connected)

    @property
    @abstractmethod
    def exceptions(self):
        """Return a tuple of camera exceptions.

        Returns
        -------
        exceptions : tuple
            Each element is an Exception subclass of exceptions
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
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'measurement_settings/exposure',
                            f'\nInfo: out of range ({min_exp}, {max_exp})')
            exposure_time = min_exp
        self.__properties['exposure'] = exposure_time
        return exposure_time

    @property
    @abstractmethod
    def extra_delay(self):
        """Return the interval spent not measuring when triggered (msec).

        This quantity is used to judge how long the camera needs to
        perform one acquisition in triggered mode. The time it takes
        is typically
            self.exposure + 1000/self.get_frame_rate()   # First frame
            + (self.n_frames - 1) * self.frame_interval  # Other frames
            + self.extra_delay
        as implemented in self.time_to_image_ready

        Returns
        -------
        extra_delay : float
            Extra time in milliseconds required by the camera
            to complete a triggering cycle.
        """
        return 0.0

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
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'measurement_settings/gain',
                            f'\nInfo: out of range ({min_gain}, {max_gain})')
            gain = min_gain
            self.settings.set('measurement_settings', 'gain', str(gain))
        self.__properties['gain'] = gain
        return gain

    @property
    def has_valid_settings(self):
        """Return whether self.settings is valid for this device."""
        return bool(self.settings)

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
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
                mode, 'camera_settings/mode', 'triggered'
                )
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
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
                f'{n_frames} [out of range ({min_n}, {max_n})]',
                'measurement_settings/n_frames', 1
                )
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
    def name_clean(self):
        """Return a version of .name suitable for file names."""
        return base.as_valid_filename(self.name)

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
        # One special case: when the ROI setting is the empty
        # string, attempt to get ROI info from the camera
        if not self.settings.get('camera_settings', 'roi', fallback='_'):
            roi = self.get_roi()
            if roi:
                self.settings.set('camera_settings', 'roi', str(roi))
                self.settings.update_file()

        try:
            # pylint: disable=redefined-variable-type  # pylint bug
            roi = self.settings.getsequence('camera_settings', 'roi',
                                            fallback=tuple())
        except _m_settings.NotASequenceError:  # Invalid ROI string
            roi = tuple()

        try:
            roi = tuple(int(r) for r in roi)
        except (TypeError, ValueError):  # Not iterable or not ints
            roi = tuple()

        _stored = self.__properties['roi']
        if _stored and _stored == roi:
            return roi

        _, size_max, *_ = limits = self.get_roi_size_limits()
        full_roi = 0, 0, *size_max

        if not self.__is_valid_roi(roi, limits):
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'camera_settings/roi', '\nInfo: ROI is invalid. '
                            f'You can use the full sensor: roi = {full_roi}')
            self.settings.set('camera_settings', 'roi', 'None')
            return full_roi

        # Adjust the ROI position to fit increments, if necessary
        roi = self.__adjust_roi_position(roi, limits)
        if not roi:
            return full_roi

        # Adjust ROI width and height to fit the binning factor
        roi = self.__adjust_roi_size_to_bin(roi, limits)
        self.__properties['roi'] = roi
        return roi

    @property
    def sensor_size(self):
        """Return the full width and height of the sensor.

        The base implementation relies on the results returned
        by .get_roi_size_limits(). This property should be
        overridden if the sensor size is larger than the
        largest ROI applicable.

        Returns
        -------
        width, height : int
            Width and height of the sensor in pixels.
        """
        return self.get_roi_size_limits()[1]

    # .settings setter
    @qtc.pyqtSlot(object)
    def set_settings(self, new_settings):
        """Set new settings for the camera.

        The new settings are accepted only if valid. Otherwise, the
        previous settings are kept. If valid, the settings are also
        loaded to the camera hardware. It is important to ensure that
        the camera is currently not waiting to for completion of an
        image acquisition or image processing step, as this could
        potentially stall the camera.

        Returns
        -------
        settings_valid : bool
            True if the new settings given were accepted.

        Parameters
        ----------
        new_settings : dict or ConfigParser or str or Path or ViPErLEEDSettings
            Whatever can be used to create a ViPErLEEDSettings.
            The following sections/options are mandatory:
            'camera_settings'/'class_name'
                Name of the concrete class that implements
                the abstract methods of CameraABC.
            'camera_settings'/'device_name'
                Name of the device as it appears in the device
                list (e.g., "DMK 33GX265").
            'measurement_settings'/'exposure'
                Exposure time in milliseconds to be used for
                acquiring each frame.
        """
        # Will recalculate bad pixel coordinates if ROI changes.
        old_roi = tuple()
        if self.settings:
            try:
                old_roi = self.roi
            except self.exceptions:
                pass

        # Checking of non-mandatory data is done in property getters.
        # The base-class implementation takes care of _mandatory_settings.
        if not super().set_settings(new_settings):
            return False

        self.disconnect_()

        if self.uses_default_settings:
            # When loading from default settings, there's no point in           # TODO: should also clear bad pixels if present
            # trying to connect: the device name is certainly wrong.
            return True

        self.connect_()  # Also loads self.settings to camera.

        if self.roi != old_roi and self.bad_pixels:
            self.bad_pixels.apply_roi()
        return True

    @property
    @abstractmethod
    def supports_trigger_burst(self):
        """Return whether the camera allows triggering multiple frames.

        This property should be overridden in concrete subclasses.
        Should the camera support trigger burst, get_n_frames_limits
        should also be overridden.

        Returns
        -------
        supported : bool
            True if the camera internally supports triggering multiple
            frames, i.e., only one 'trigger_now()' call is necessary to
            deliver all the frames needed.
        """
        return False

    @property
    def time_to_image_ready(self):
        """Return the ideal interval (msec) to fully acquire an image.

        Returns
        -------
        time_to_image_ready : float
            The time interval between when .trigger_now() is called
            and when a full image is ready to be processed (i.e.,
            the time needed to acquire all .n_frames). This is also
            the time interval between when camera.busy == True (when
            triggering) and when camera.busy == False again. This
            attribute is an estimate.
        """
        exposure = self.exposure
        frame_rate = self.get_frame_rate()

        # Could use self.frame_interval, but we save one call
        # to .get_frame_rate() by calculating it explicitly
        interval = max(exposure, 1000/frame_rate)
        return (exposure + 1000/frame_rate         # first frame
                + interval * (self.n_frames - 1)   # all other frames
                + self.extra_delay)

    @staticmethod
    def is_bad_pixels_error(error_info):
        """Return whether error_info relates to 'bad pixels'."""
        error_code, error_msg = error_info
        try:
            error = QObjectSettingsErrors.from_code(error_code)
        except AttributeError:
            return False

        error_msg = error_msg.replace('_', ' ')
        if (error is QObjectSettingsErrors.INVALID_SETTINGS
                and 'bad pixel' in error_msg):
            return True
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
            if isinstance(in_config, str):
                equal = operator.eq
            elif isinstance(in_config, (Sequence, np.ndarray)):
                equal = np.allclose
            else:
                equal = np.isclose
            if not equal(in_config, in_device):
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
            self.emit_error(DeviceABCErrors.DEVICE_NOT_FOUND, self.name)
            self.connected = False
            return
        self.load_camera_settings()
        self.connected = True

    @qtc.pyqtSlot()
    def disconnect_(self):
        """Disconnect the device."""
        self.__reported_errors = set()
        self.stop()
        self.close()
        self.connected = False

    def get_settings_handler(self):
        """Return a SettingsHandler object for displaying settings.

        This method can be extended in subclasses, i.e., do
        handler = super().get_settings_handler(), and then add
        appropriate sections and/or options to it using the
        handler.add_section, and handler.add_option methods.

        Use the QNoDefaultPushButton from the widgets.buttons
        module in order to prevent any button from being set
        as the default button of the dialog.

        The base-class implementation returns a handler that
        already contains the following settings:
        - 'camera_settings'/'mode' (read-only, isolated option)
        - 'measurement_settings'/'exposure'
        - 'measurement_settings'/'gain'
        - 'measurement_settings'/'n_frames'
        - 'camera_settings'/'roi'
        - 'camera_settings'/'binning'
        - 'camera_settings'/'bad_pixel_path' (read-only)

        Returns
        -------
        handler : SettingsHandler
            The handler used in a SettingsDialog to display the
            settings of this controller to users.
        """
        handler = super().get_settings_handler()
        handler.add_option('camera_settings', 'mode',
                           handler_widget=qtw.QLabel(self.mode.capitalize()),
                           tags=SettingsTag.READ_ONLY | SettingsTag.REGULAR)
        handler.add_section('measurement_settings', display_name="Acquisition",
                            tags=SettingsTag.MEASUREMENT | SettingsTag.REGULAR)
        handler.add_section('camera_settings', display_name="Image Properties",
                            tags=SettingsTag.REGULAR)

        # pylint: disable=redefined-variable-type
        # Triggered for _widget. While this is true, it clear what
        # _widget is used for in each portion of filling the handler

        # Exposure time in ms
        _widget = CoercingDoubleSpinBox(soft_range=self.get_exposure_limits(),  # TODO: use prefix-adaptive widget
                                        decimals=2, step=5.0, suffix=' ms')
        _widget.setStepType(_widget.AdaptiveDecimalStepType)
        _widget.setAccelerated(True)
        _widget.setMinimum(0)  # Can't have negative exposure
        _tip = (
            "<nobr>Exposure time used for each frame. For LEED\u2011IV "
            "videos it is best</nobr> to choose the exposure time <b>as long "
            "as possible while not causing the camera to saturate</b>. You "
            "can find the best exposure time by (i) selecting the electron "
            "energy where you have the maximum intensity of spot(s), and "
            "(ii) change the exposure time to have close-to-saturated images"
            )
        handler.add_option('measurement_settings', 'exposure',
                           handler_widget=_widget, tooltip=_tip)

        # Gain in dB                                                            # TODO: add also a gain factor!
        _widget = CoercingDoubleSpinBox(soft_range=self.get_gain_limits(),
                                        decimals=1, suffix=" dB")
        _widget.setAccelerated(True)
        _tip = (
            "<nobr>ADC gain used for each frame. For LEED\u2011IV videos it "
            "is</nobr> best to choose the gain <b>as small as possible</b>, "
            "as larger gains <b>amplify noise</b>. Use a longer exposure time "
            "(without causing saturation) rather than a larger gain"
            )
        handler.add_option('measurement_settings', 'gain',
                           handler_widget=_widget, tooltip=_tip)

        # n_frames
        if self.get_n_frames():
            _range = self.get_n_frames_limits()
        else:
            _range = (1, float('inf'))
        _widget = CoercingSpinBox(soft_range=_range)
        _widget.setMinimum(0)
        _widget.setAccelerated(True)
        _tip = (
            "<nobr>Number of frames to be averaged for saving images.</nobr> "
            "This is <b>used only in triggered mode</b> and not when snapping "
            "an image. For LEED\u2011IV videos it is a good idea to choose "
            "the number of frames such that acquiring images takes roughly "
            "the same amount of time as measuring other quantities (e.g., I0)"
            )
        handler.add_option('measurement_settings', 'n_frames',
                           handler_widget=_widget, tooltip=_tip,
                           display_name="No. frames")

        if not self.connected:
            return handler
        # ROI
        _widget = ROIEditor()
        size_min, size_max, _d_size, _d_offset = self.get_roi_size_limits()

        # pylint: disable=no-value-for-parameter
        # Seems a bug related to calling a function with variadic
        # arguments when it has a signature without.
        _widget.set_increments(*_d_offset, *_d_size)
        # pylint: enable=no-value-for-parameter

        _widget.set_ranges(size_min, size_max)
        _widget.original_roi = self.roi
        _widget.set_ = _widget.set_from_string
        _widget.get_ = _widget.get_as_string
        _widget.notify_ = _widget.roi_changed
        _tip = (
            "<nobr>Region of interest used to crop camera images.</nobr> "
            "Since a LEED screen is circular, it is best to (i) place "
            "the camera so that the <b>LEED screen fills most of the field "
            "of view</b>, and (ii) crop the image to a <b>square-ish "
            "shape.</b> This reduces file sizes and computation times."
            )
        handler.add_option('camera_settings', 'roi', handler_widget=_widget,
                           display_name="ROI", label_alignment='bottom',
                           tooltip=_tip)

        # Binning factor                                                        # TODO: add a mention of the image size after processing. Also, add some warning if > 10
        _widget = CoercingSpinBox(soft_range=self.get_binning_limits())
        _widget.setMinimum(0)
        _tip = (
            "<nobr>Binning factor used to reduce the size of images. Binning "
            "consists</nobr> in averaging pixel intensities over a (bin "
            "\u00d7 bin) mesh. It is best if images in LEED\u2011IV videos "
            "are reduced to roughly 500\u2013600 pixels. This saves space and "
            "makes data extraction faster, without compromising the quality"
            )
        handler.add_option('camera_settings', 'binning',
                           handler_widget=_widget, tooltip=_tip)

        # Bad pixels path + note on using Tools...
        _tip = (
            "<nobr>Path to folder containing bad\u2011pixels files. This "
            "path</nobr> can be set while finding bad pixels <b>using the "
            "Tools\u2011>Find bad pixels... menu</b>. Bad pixels are "
            "particularly bright, dark, or noisy pixels that should be "
            "corrected to improve data quality"
            )
        if self.settings.misses_settings(('camera_settings',
                                          'bad_pixels_path')):
            # We may ignore errors related to the bad-pixels path
            self.settings['camera_settings']['bad_pixels_path'] = ''
        handler.add_option('camera_settings', 'bad_pixels_path',
                           handler_widget=PathSelector,
                           tags=SettingsTag.READ_ONLY, tooltip=_tip)
        return handler

    @abstractmethod
    def list_devices(self):
        """Return a list of available devices.

        This method must return a list of SettingsInfo instances. The
        SettingsInfo class is located in the classes.abc module. Each
        camera is represented by a single SettingsInfo instance. The
        SettingsInfo object must contain a .unique_name,
        .has_hardware_interface which is True if the device has a
        hardware interface present, and a dict holding .more information
        about the device. If there is no additional information about
        the camera, then this dict can be empty. The information
        contained within a SettingsInfo must be enough to determine a
        suitable settings file for the device from it. Subclasses should
        raise a DefaultSettingsError if they fail to create instances
        from the default settings.

        Returns
        -------
        devices : list
            Each element is a SettingsInfo instance containing the
            name of a camera and additional information as a dict.
        """
        return

    def load_camera_settings(self):
        """Load settings from self.settings into the camera.

        This method should only be overridden if the order in which
        setter-method calls are performed is inappropriate. The default
        order of the calls is:
            set_exposure(),  set_gain(),  set_mode()
            <setters associated to self.hardware_supported_features>
            [set_binning()], [set_n_frames()], [set_roi()]
        The last 3 are called only if they are not hardware-supported.
        Should the method be overridden, the new implementation should
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

        Subclasses can use self.name to access the device name for
        this camera.

        Returns
        -------
        successful : bool
            True if the device was opened successfully.
        """
        return False

    @qtc.pyqtSlot()
    def update_bad_pixels(self):
        """Update the list of bad pixels from the settings.

        This function should be explicitly called after the
        path containing bad pixels files is changed in the
        settings. Calling this function may take a little
        time, as bad pixel coordinates are read from scratch
        and recalculated.

        Emits
        -----
        error_occurred(QObjectSettingsErrors.INVALID_SETTINGS)
            If settings file does not have a 'bad_pixels_path' option
            in section 'camera_settings'.
        error_occurred(QObjectSettingsErrors.INVALID_SETTINGS)
            If the path specified does not contain a valid bad-pixels
            file for this camera.
        """
        bad_pix_path = self.settings.get("camera_settings", "bad_pixels_path",
                                         fallback='')
        if not bad_pix_path:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'camera_settings/bad_pixels_path',
                            '\nInfo: No bad_pixels_path found.')
            return
        try:
            self.__bad_pixels.read(bad_pix_path)
        except (FileNotFoundError, ValueError) as err:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'camera_settings/bad_pixels_path',
                            f'\nInfo: {err}')
            return
        self.__bad_pixels.apply_roi()

    @abstractmethod
    def get_binning(self):
        """Return the binning factor internally used for averaging.

        This method should be overridden to query the camera for the
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

        Do not override this method, unless the camera internally
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
            If this method is not overridden for a camera
            that natively supports binning.
        """
        if self.get_binning():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports binning, "
                "but self.set_binning() was not overridden."
                )
        binning_factor = 1 if no_binning else self.binning
        self.process_info.binning = binning_factor

    def get_binning_limits(self):
        """Return the minimum and maximum value for binning.

        This method should be overridden only if the camera
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
            If this method is not overridden for a camera
            that natively supports binning.
        """
        if self.get_binning():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports binning, "
                "but self.get_binning_limits() was not overridden."
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

        This method must be overridden to call the internal driver
        function that sets the appropriate exposure time in the camera
        device. Make sure to self.open() before setting.

        Subclasses must use the exposure time (in milliseconds) that
        self.exposure returns (i.e., the one in self.settings).

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
        return 0, np.inf

    @abstractmethod
    def get_frame_rate(self):
        """Return the number of frames delivered per second.

        This property should be overridden in concrete subclasses.

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
        return self.settings.getfloat('measurement_settings', 'gain')

    @abstractmethod
    def set_gain(self):
        """Set the gain of the camera in decibel.

        This method must be overridden to call the internal driver
        function that sets the appropriate gain in the camera device.
        Make sure to self.open() before setting.

        Subclasses must use the gain factor (in decibel) returned
        by self.gain (i.e., the one in self.settings).
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
        return -np.inf, np.inf

    @abstractmethod
    def get_mode(self):
        """Return the mode set in the camera.

        This method should be overridden to query the camera for the
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

        This method should be overridden to instruct the camera
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

        This method should be overridden to query the camera
        for the current number of frames internally used for
        averaging. Make sure to self.open() before querying.

        Returns
        -------
        n_frames : int
            Number of frames used for averaging. Returns 0 if the
            camera does not internally support frame averaging.
        """
        return 0

    def set_n_frames(self):
        """Set the number of frames internally used for averaging.

        This method should be overridden to instruct the camera to use
        self.n_frames frames for averaging internally, before returning
        the image. Make sure to self.open() before setting.

        Should the camera not natively support frame averaging, this
        method should not be overridden, as frame averaging will be
        done via software during image processing.
        """
        if self.get_n_frames():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports frame "
                "averaging, but self.set_n_frames() was not overridden"
                )
        self.process_info.n_frames = self.n_frames

    def get_n_frames_limits(self):
        """Return the minimum and maximum number of frames supported.

        This method should be overridden only if the camera
        internally supports frame averaging. Make sure to
        .open() before querying the camera.

        Returns
        -------
        min_n_frames, max_n_frames : int
            Minimum and maximum number of frames usable for averaging

        Raises
        ------
        NotImplementedError
            If this method is not overridden for cameras that
            natively support averaging over multiple frames, or
            for those that natively support a trigger-burst mode
        """
        if self.get_n_frames():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports frame averaging,"
                " but self.get_n_frames_limits() was not overridden"
                )
        if self.supports_trigger_burst:
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports trigger burst, "
                "but self.get_n_frames_limits() was not overridden"
                )
        return 1, np.inf

    @abstractmethod
    def get_roi(self):
        """Return the region of interest set in the camera.

        This method should be overridden in concrete subclasses to
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

        This method should be overridden only if the camera natively
        supports setting a region of interest at the hardware level.
        Subclasses must retrieve the ROI information using self.roi.
        Make sure to .open() before setting.

        If the camera does not internally support setting a ROI, then
        the ROI setting is used to crop the frames during processing.

        Parameters
        ----------
        no_roi : bool, optional
            If True, set the ROI to the full size of the sensor rather
            than using the value from the settings. Default is False.

        Raises
        -------
        NotImplementedError
            If this method is not overridden when the camera
            natively supports applying a region of interest.
        """
        if self.get_roi():
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports setting a "
                "region of interest, but self.set_roi() was not overridden"
                )
        roi = tuple() if no_roi else self.roi
        self.process_info.roi = roi

    @abstractmethod
    def get_roi_size_limits(self):
        """Return minimum, maximum and granularity of the ROI.

        This method must be overridden in concrete subclasses.
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
    @qtc.pyqtSlot()
    @qtc.pyqtSlot(object)
    def start(self, *_):
        """Start the camera.

        The camera is started in the mode returned by self.mode.

        This method should be extended in concrete subclasses, i.e.,
        super().start() must be called in the reimplementation. The
        reimplementation must self.started.emit() after the driver
        has been appropriately started.

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
        if self.exposure < frame_interval - 0.1:
            print(                                                              # TODO: should become a non-fatal warning
                f"WARNING: Exposure ({self.exposure} ms) of camera "
                f"{self.name} is shorter than the time it takes to "
                f"deliver frames ({frame_interval:.2f} ms). Increase "
                "the exposure time to avoid wasting time."
                )

    @abstractmethod
    @qtc.pyqtSlot()
    def stop(self):
        """Stop the camera.

        The camera may not be immediately stopped if there is any
        image processing currently running. In this case, stopping
        will be attempted multiple times with a 100 ms interval.
        Once the camera has effectively stopped, the .stopped()
        signal is emitted.

        This method should be extended in concrete subclasses. The
        pattern of the extended method must be of the form:
            if not super().stop():
                return False
            ... actually stop self.driver ...
            self.stopped.emit()
            return True

        Notice that the base implementation DOES NOT EMIT .stopped().

        Returns
        -------
        stopped : bool
            True if the call led to stopping the camera. False if
            the camera was not running, or if it cannot be stopped
            yet because some frames are missing.
        """
        if not self.is_running:
            return False

        # Delay the stopping as long as there are image
        # processors that have not finished their tasks,
        # i.e., we are still waiting for some frames, or
        # waiting for images to be processed and saved.
        if (self.__process_thread.isRunning()
                and any(p.busy for p in self.__image_processors)):
            if not self.__retry_stop_timer.isActive():
                self.__retry_stop_timer.start(100)
            return False

        # Now it is safe to clean up and stop
        self.__retry_stop_timer.stop()
        self.__timeout.stop()
        self.process_info.clear_times()
        self.n_frames_done = 0
        self.busy = False
        if self.__process_thread.isRunning():
            self.__process_thread.quit()
        base.safe_disconnect(self.frame_ready, self.__on_frame_ready)
        return True

    @abstractmethod
    @qtc.pyqtSlot()
    def trigger_now(self):
        """Start acquiring one (or more) frames now.

        This method should be extended in concrete subclasses to
        signal the camera that acquisition of one (or multiple)
        frames should start now.

        Subclasses should proceed only if super().trigger_now()
        is True. The base-class implementation will (i) emit an
        error_occurred if the camera mode is inappropriate (i.e.,
        the camera is not in triggered mode); (ii) mark the camera
        as busy, (iii) reset n_frames_done to zero, and (iv) start
        a timeout timer.

        Returns
        -------
        successfully_triggered : bool
            True if the camera was successfully triggered.

        Emits
        -----
        error_occurred(CameraErrors.UNSUPPORTED_OPERATION)
        """
        if self.mode != 'triggered':
            self.emit_error(CameraErrors.UNSUPPORTED_OPERATION,
                            'trigger', 'live')
            return False
        self.busy = True

        if self.supports_trigger_burst:
            self.process_info.clear_times()

        # Start the timeout timer: will fire if it takes
        # more than 5 s longer than the expected time to
        # receive all frames needed for averaging
        self.__timeout.start(int(self.time_to_image_ready + 5000))
        return True

    def __adjust_roi_position(self, roi, limits):
        """Check and adjust the position in roi fits the increments."""
        pos_increments = limits[-1]
        *roi_pos, roi_w, roi_h = roi
        if not any(p % i for p, i in zip(roi_pos, pos_increments)):
            # Offsets are OK
            return roi

        # Offsets are not OK. Let's try to fix them.
        roi_pos = ((p // i) * i for p, i in zip(roi_pos, pos_increments))
        roi = (*roi_pos, roi_w, roi_h)
        if self.__is_valid_roi(roi, limits):
            # Successfully adjusted
            self.settings.set('camera_settings', 'roi', str(roi))
            return roi

        self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                        'camera_settings/roi',
                        f'\nInfo: ROI {roi} is invalid after '
                        'adjusting top-left corner position')
        self.settings.set('camera_settings', 'roi', 'None')
        return None

    def __adjust_roi_size_to_bin(self, roi, limits):
        """Check and adjust roi width/height to fit bin size."""
        _bin = self.binning
        if _bin == 1:  # Nothing to check. Will always be OK
            return roi

        *roi_offsets, roi_w, roi_h = roi
        if not (roi_w % _bin or roi_h % _bin):
            # Fits binning
            return roi

        # Does not fit. Fix it, if possible, otherwise
        # let the image processor take care of this
        new_roi_w, new_roi_h = ((s // _bin) * _bin for s in (roi_w, roi_h))
        error = (CameraErrors.BINNING_ROI_MISMATCH,
                 roi_w, roi_h, self.binning, new_roi_w, new_roi_h)
        if not self.__already_reported(error):
            # self.emit_error(*error)                                           # TODO: make this a non-critical warning!
            self.__reported_errors.add(error)

        # Update the ROI in the settings only if the
        # new ROI is a valid one for the camera. The
        # image processor will anyway crop off the
        # extra pixels on the lower-right corner of
        # the image to apply binning.
        new_roi = (*roi_offsets, new_roi_w, new_roi_h)
        if self.__is_valid_roi(new_roi, limits):
            self.settings.set('camera_settings', 'roi', str(new_roi))
            roi = new_roi
        return roi

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

    def __is_valid_roi(self, roi, limits):
        """Check that ROI is OK and fits with limits.

        Parameters
        ----------
        roi : tuple
            The region of interest to be checked
        limits : tuple
            Limits for ROI, used for the check. Elements are
            2-tuples with the following meaning: size_min,
            size_max, size_delta, pos_delta. pos_delta is
            not used here.

        Returns
        -------
        roi_OK : bool
            True if roi is a valid ROI and if it fits limits.
        """
        if not roi:
            return False
        if len(roi) != 4:
            return False
        pos, size = roi[:2], roi[2:]
        size_min, size_max, size_delta, _ = limits

        roi_outside = any(p < 0 or p + s > s_max
                          for p, s, s_max in zip(pos, size, size_max))
        roi_small = any(s < s_min for s, s_min in zip(size, size_min))
        roi_not_multiple = any(
            (s - s_min) % d
            for s, s_min, d in zip(size, size_min, size_delta)
            )
        return not (roi_outside or roi_small or roi_not_multiple)

    @qtc.pyqtSlot()
    def __on_all_frames_done(self):
        """Disconnect the current processor to prevent duplicate processing."""
        processor = self.sender()
        if processor:
            # Disconnect the sender processor to prevent
            # double processing of the next frame (issue #122)
            base.safe_disconnect(self.__process_frame, processor.process_frame)

    @qtc.pyqtSlot(np.ndarray)
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

        if not self.busy:
            # Some frames arrived, but we haven't asked for them:
            # we should be busy if we are in triggered mode and
            # went through self.trigger_now()
            return

        if not self.n_frames_done:
            processor = ImageProcessor()
            processor.image_processed.connect(self.image_processed)
            processor.image_saved.connect(self.__on_image_saved)
            processor.all_frames_acquired.connect(self.__on_all_frames_done)
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
            try:
                self.__image_processors.remove(processor)
            except ValueError:
                # Not in list
                pass

        # Emit image_saved only if a file was actually saved
        if fpath:
            self.image_saved.emit(fpath)

    @qtc.pyqtSlot()
    def __on_timed_out(self):
        """Report a timeout error while in triggered mode."""
        self.emit_error(CameraErrors.TIMEOUT, self.name, 5)
