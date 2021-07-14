"""Module cameraabc of viperleed.?????.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-12
Author: Michele Riva
Author: Florian Doerr

"""

from abc import ABCMeta, abstractmethod
import sys
import ctypes

import numpy as np
from PyQt5 import QtCore as qtc

from viperleed.guilib.measure import hardwarebase

# check this out!
# def bin2d(a, K):
    # m_bins = a.shape[0]//K
    # n_bins = a.shape[1]//K
    # return a.reshape(m_bins, K, n_bins, K).sum(3).sum(1)

def bin_image(image, bin_factor):
    """
    Returns a binned copy of image given a binning factor. Notice that the
    bottom-right corner of the image will be cropped if the image dimensions
    are not an integer multiple of the binning factor

    Parameters
    ----------
    image: np.ndarray
    bin_factor: int

    Returns
    -------
    image: np.dnarray, dtype=float64
    """

    # LOOK AT https://stackoverflow.com/questions/28906298/python-code-to-quickly-reduce-the-resolution-of-an-image-using-numpy
    # for a possibly faster implementation. Over there, it may be possible to
    # also speed it up a bit by using np.sum *(1/bin_fat**2)

    # n, m = (image.shape//bin_factor)*bin_factor
    n = (image.shape[0]//bin_factor)*bin_factor
    m = (image.shape[1]//bin_factor)*bin_factor
    if (n, m) != image.shape:
        image = image[0:n, 0:m]

    image = image.reshape(n, m//bin_factor,
                          bin_factor).sum(axis=2, dtype=np.int32).T
    image = image.reshape(m//bin_factor, n//bin_factor,
                          bin_factor).sum(axis=2, dtype=np.int32).T
    return image*(1/bin_factor**2)


class CameraErrors(hardwarebase.ViPErLEEDErrorEnum):
    """Data class for basic serial errors not available in QSerialPort."""
    INVALID_SETTINGS = (200,
                        "Invalid camera settings: Required "
                        "settings {!r} missing or values "
                        "inappropriate. Check configuration file.")
    MISSING_SETTINGS = (201,
                        "Camera cannot operate without settings. "
                        "Load an appropriate settings file before "
                        "proceeding.")
    CAMERA_NOT_FOUND = (202, "Could not find camera {}.")

    SETTINGS_MISMATCH = (203,
                         "Different {} settings found in camera and "
                         "configuration file: camera={}, settings={}.")
    INVALID_SETTING_WITH_FALLBACK = (204,
                                     "Invalid camera settings value {} for "
                                     "setting {!r}. Using {} instead. Consider "
                                     "fixing your configuration file.")
    UNSUPPORTED_OPERATION = (205, "Cannot {} in {} mode. Switch mode before.")

class CallbackData(ctypes.Structure):
    """ Data passed to the callback function. """
    # Each subclass must define a _fields_ attribute. _fields_ must be a list of
    # 2-tuples, containing a field name and a field type. The field type must be
    # a ctypes type like c_int, or any other derived ctypes type: structure,
    # union, array, pointer.
    # WE'RE NOT DOING THIS, SO PROBABLY WE DO NOT REALLY NED THE CTYPES.STRUCTURE
    # AT ALL
    def __init__(self):
        self.image = 0
        self.busy = False
        self.directory = None
        self.dac_busy = False
        # keys 'binning' and 'roi' should be left to None unless the camera
        # does not support binning and ROI definitions
        # WHERE IS THE KEY roi????
        self.process_params = {'filename_prefix': None,
                               'filename_energy': None,
                               'directory': None,
                               'n_frames': 1,
                               'counter_frames': 0,
                               'binning': None,
                               'roi': None}  # added. Check if OK
        self.camera_lib = None  # Reference to the camera object


class CameraABC(metaclass=ABCMeta):
    """Abstract base class for ViPErLEED cameras."""

    error_occurred = qtc.pyqtSignal(tuple)

    _mandatory_settings = [  # TODO: is it better this way or the way we do it in serial and controller
            ('camera_settings', 'class_name'),
            ('camera_settings', 'device_name'),
            ('measurement_settings', 'camera_exposure'),
            # ('measurement_settings', 'n_frames'),  # default=1
            ]

    def __init__(self, *args, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        *args : object
            Unused positional parameters
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
        self.__settings = None
        self.set_settings(kwargs.get("settings", None))

        self.callback_data = CallbackData()

    # @staticmethod
    # def camera_from_ini(camera_class_name):
        # """
        # Given a camera class name as a string, returns an instance of the class
        # if implemented. Raises ValueError otherwise.
        # """
        # try:
            # camera_class = getattr(sys.modules[__name__], camera_class_name)
        # except AttributeError:
            # raise ValueError(f"camera class {camera_class_name} not found")
        # return camera_class()

    def __get_settings(self):
        """Return the settings used for the camera."""
        return self.__settings

    def set_settings(self, new_settings):
        """Set new settings for the camera.

        Parameters
        ----------
        new_settings : dict or ConfigParser
            The new settings. The following sections/options
            are mandatory:
            'camera_settings'/'class_name'
                Name of the concrete class that implements
                the abstract methods of the base ABC.
            'camera_settings'/'device_name'
                Name of the device as it appears in the device
                list (e.g., "DMK 33GX265").
            'measurement_settings'/'camera_exposure'
                Exposure time in milliseconds to be used for
                acquiring each frame.
        """
        new_settings, invalid = hardwarebase.config_has_sections_and_options(
            self,
            new_settings,
            self._mandatory_settings
            )
        if invalid:
            print("invalid in viperleed", invalid)
            if hasattr(invalid, '__len__'):
                for setting in invalid:
                    (error_code,
                     error_msg) = CameraErrors.INVALID_SETTINGS
                    error_msg = error_msg.format(setting)
                    self.error_occurred.emit((error_code, error_msg))
            else:
                self.error_occurred.emit(CameraErrors.MISSING_SETTINGS)
        self.__settings = new_settings

        # Store constant optional data with fallback values
        # also checking that, if given, they are acceptable
        # (1) number of frames for averaging
        if self.n_frames <= 0:
            error_code, error_msg = CameraErrors.INVALID_SETTING_WITH_FALLBACK
            error_msg = error_msg.format(self.n_frames,
                                         'measurement_settings/n_frames',
                                         1)
            self.error_occurred.emit((error_code, error_msg))
            self.settings.set('measurement_settings', 'n_frames', '1')

        # (2) measurement mode
        if self.mode not in ('live', 'triggered'):
            error_code, error_msg = CameraErrors.INVALID_SETTING_WITH_FALLBACK
            error_msg = error_msg.format(self.mode,
                                         'camera_settings/mode',
                                         'triggered')
            self.error_occurred.emit((error_code, error_msg))
            self.settings.set('camera_settings', 'mode', 'triggered')

        # (3) region of interest (None or a 4-tuple)
        if self.roi is not None and len(roi) != 4:
            error_code, error_msg = CameraErrors.INVALID_SETTING_WITH_FALLBACK
            error_msg = error_msg.format(self.roi,
                                         'camera_settings/roi',
                                         'triggered')
            self.error_occurred.emit((error_code, error_msg))
            self.settings.set('camera_settings', 'roi', 'None')

    settings = property(__get_settings, set_settings)

    @property
    def exposure(self):
        """Return the exposure time from the settings.

        Returns
        -------
        exposure_time : float
            Exposure time in milliseconds
        """
        return self.settings.getfloat('measurement_settings', 'camera_exposure')

    @property
    def gain(self):
        """Return the gain from the settings.

        Returns
        -------
        gain : float
            Gain in decibel. If there is no gain setting,
            this method returns 0.
        """
        return self.settings.getfloat('measurement_settings',
                                      'camera_gain',
                                      fallback=0)

    @property
    def mode(self)
        """Return the camera mode from the settings.

        Returns
        -------
        mode : {'live', 'triggered'}
            The mode from the current settings
        """
        return self.settings.get('camera_settings',
                                 'mode',
                                 fallback='triggered')

    @property
    def n_frames(self):
        """Return the number of frames used for averaging."""
        return self.settings.getint('measurement_settings',
                                    'n_frames',
                                    fallback=1)

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
        except (SyntaxError, ValueError):
            error_code, error_msg = CameraErrors.INVALID_SETTINGS
            error_msg = error_msg.format('camera_settings/roi')
            self.error_occurred.emit((error_code, error_msg))
            roi = None
        if roi is not None:
            try:
                roi = tuple(int(r) for r in roi)
            except (TypeError, ValueError):
                roi = None
        self.settings.set('camera_settings', 'roi', str(roi))
        return roi

    def connect(self):
        """Connect to the camera."""
        if not self.open():
            error_code, error_msg = CameraErrors.CAMERA_NOT_FOUND
            error_msg = error_msg.format(self.name)
            self.error_occurred.emit((error_code, error_msg))
            return

        self.load_camera_settings()
        self.check_loaded_settings()

    def load_camera_settings(self):
        """Load settings from self.settings into the camera."""
        self.set_exposure()
        self.set_gain()
        self.set_mode()

    def check_loaded_settings(self):
        """Check that camera and confguration settings are the same.

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
        to_check = {
            'exposure time': (self.exposure, self.get_exposure()),
            'gain': (self.gain, self.get_gain()),
            'mode': (self.mode, self.get_mode()),
            }

        error_code, error_txt = CameraErrors.SETTINGS_MISMATCH
        error_msg = []
        for setting, (config, device) in to_check.items():
            if config != device:
                error_msg.append(error_txt.format(setting, device, config))

        # The next settings may be natively unsupported by the
        # camera, and be done at the software level while
        # procesing images. Unsupported features return None
        may_be_unsupported = {
            'roi': (self.roi, self.get_roi()),
            'binning': (self.binning, self.get_binning()),
            }
        for setting, (config, device) in may_be_unsupported.items():
            if device is None:
                continue
            if device != config:
                error_msg.append(error_txt.format(setting, device, config))

        if error_msg:
            self.error_occurred.emit((error_code, '\n'.join(error_msg)))
            return False
        return True

    @abstractmethod
    def list_devices(self):
        """Return a list of avaliable device names.

        Returns
        -------
        devices : list
            Each element is a string representing
            the name of a camera device.
        """
        return

    @abstractmethod
    def open(self)
        """Open the camera device.

        This method is guaranteed to be called exactly once
        every time .connect() is called.

        The reimplementation can use self.name to access the
        device name for this camera.

        Returns
        -------
        successful : bool
            True if the device was opened successfully.
        """
        return False

    @abstractmethod
    def initialize(self, *args, **kwargs):
        """Initialize camera.

        TODO: see what we're using this for and add doc
        """
        return

    def is_busy(self):
        """TODO: doc."""
        return self.callback_data.busy

    @abstractmethod
    def get_exposure(self):
        """Get the exposure time from the camera device.

        Returns
        -------
        exposure_time : float
            Exposure time in milliseconds.
        """
        return

    @abstractmethod
    def set_exposure(self):
        """Set the exposure time for one frame.

        This method must be rimplemented to call the internal
        driver function that sets the appropriate exposure time
        in the camera device.

        The exposure time in the settings should be accessed as
        a float in milliseconds using self.exposure.
        """
        return

    @abstractmethod
    def get_gain(self):
        """Get the gain from the camera device.

        Returns
        -------
        gain : float
            Gain in decibel.
        """
        return

    @abstractmethod
    def set_gain(self, new_gain):
        """Set the gain of the camera in decibel.

        This method must be rimplemented to call the internal
        driver function that sets the appropriate gain in the
        camera device.

        The gain in the settings should be accessed as a float
        in decibel using self.gain.

        Parameters
        ----------
        new_gain : float
            Gain in decibel
        """
        return

    @abstractmethod
    def get_mode(self):
        """Return the mode set in the camera.

        This method should be reimplemented to query the camera
        for the current mode, and return the appropriate string.
        Should the camera not support natively a live vs. triggered
        mode, this function can return self.mode.

        Returns
        -------
        mode : {'live', 'triggered'}
            The mode the camera is operating in. Continuous mode
            corresponds to receiving a 'stream' of images at
            well-defined time intervals (i.e., synchronous), while
            'triggered' is asyncrhronous: the camera returns a frame
            only when asked by self.trigger_now().
        """
        return

    @abstractmethod
    def set_mode(self):
        """Set the camera mode from self.settings.

        This method should be reimplemented to instruct the camera
        to switch to a desired mode, which should be retrieved with
        self.mode from the settings.

        Should the camera not support natively a live vs. triggered
        mode, live mode can be simulated with a QTimer.
        """
        # TODO: we may want to not make this an abstractmethod,
        # and have get_mode return None if the camera does not
        # support live-view. Then the base implementation would
        # prepare the system for a live-view simulation using
        # QTimer
        return

    @abstractmethod
    def get_n_frames(self):
        """Return the number of frames internally used for averaging.

        This method should be reimplemented to query the camera
        for the current number of frames internally used for
        averaging.
            
        Returns
        -------
        n_frames : int or None
            Number of frames used for averaging. Returns None if the
            camera does not support internally frame averaging.
        """
        return
    
    def set_n_frames(self):
        """Set the number of frames internally used for averaging.
        
        This method should be reimplemented to instruct the camera
        to use self.n_frames frames for averging internally, before
        returning the image.
        
        Should the camera not natively support frame averaging, this
        method should not be reimplemented, as frame averaging will
        be done via software during image processing.
        """
        if self.get_n_frames() is not None:
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports frame "
                "averaging, but self.set_n_frames() was not reimplemented"
                )
        self.callback_data.process_params['n_frames'] = self.n_frames

    @abstractmethod
    def trigger_now(self):
        """Start acquiring one (or more) frames now.

        This method should be reimplemented in concrete subclasses
        to signal the camera that acquisition of one (or multiple)
        frames should start now.

        The reimplementation must call super().trigger_now() before,
        as this will emit an error_occurred if the camera mode is
        inappropriate (i.e., the camera is not in triggered mode).

        Returns
        -------
        None

        Emits
        -----
        error_occurred(CameraErrors.UNSUPPORTED_OPERATION)
        """
        if self.mode != 'triggered':
            error_code, error_msg = CameraErrors.UNSUPPORTED_OPERATION
            error_msg = error_msg.format('trigger', 'continuous')
            self.error_occurred.emit((error_code, error_msg))

    # TODO: keep going from here

    @abstractmethod
    def get_roi(self):
        """Return the region of interest set in the camera.

        This method should be reimplemented in concrete
        subclasses to query the camera for the current
        settings of the region of interest. If the camera
        does not support setting a region of interest this
        method should return None.

        Returns
        -------
        tuple or None
            The settings of the region of interest. If a tuple,
            the elements are:
            roi_x, roi_y : int
                Coordinates of the top-left pixel. Zero is the
                topmost/leftmost pixel.
            roi_width, roi_height : int
                Width and height of the region of interest in pixels
        """
        return

    def set_roi(self):
        """Set up region of interest in the camera from self.settings.

        If the camera does not support setting internally a ROI, the
        the ROI setting should be used in the callback function that
        processes the images to crop the images to the ROI.
        
        Returns
        -------
        None
        """
        if self.get_roi() is not None:
            raise NotImplementedError(
                f"{self.__class__.__name__} natively supports setting a "
                "region of interest, but self.set_roi() was not reimplemented"
                )
        self.callback_data.process_params['roi'] = self.roi

################################################################################
# TODO: have the callback rather emit a frame_ready signal that is always
# connected to the very same image processing function

    @abstractmethod
    def set_binning(self, bin_factor):
        """Define binning for images.

        If the camera does not support setting binning internally, the
        binning values set here should be used in the callback function
        that processes the images.

        Parameters
        ----------
        bin_factor : int
            Binning will be executed on a (bin_factor, bin_factor) mesh
        """
        return

    @abstractmethod
    def set_callback(self, on_frame_ready):
        """Define the frame-ready callback in the camera driver.

        If the camera does not support having a callback function,
        a similar behavior can be obtained using an appropriate
        pyqtSignal, emitted as soon as a frame has been acquired.

        Parameters
        ----------
        on_frame_ready : callable
            Function that will be called in triggered mode
            after each frame is completed.
        """
        return

    @abstractmethod
    def snap_image(self):
        """Return a snapshot while in live-view mode.

        This method can be called only when not runnning in triggered
        mode, i.e., after self.enable_trigger(False).

        Returns
        -------
        snapshot : numpy.ndarray
            A snapshot of the current frame
        """
        return

    @abstractmethod
    def set_imagepath(self):  # TODO: useless. takes no arguments. Better use a hard-coded (relative) temp directory
        """Set the directory where images are saved."""
        return

    @abstractmethod
    def start_camera(self, live_mode_active=False):  # TODO: rename to start only
        """Start the camera.

        The camera can be started in live-view mode or in triggered mode.

        Parameters:
        ----------
        live_mode_active : bool, optional
            If True, the camera is in live-view mode. Default is False.
        """
        return

    @abstractmethod
    def stop_camera(self):    # TODO: rename to stop only
        """
        Stops the camera and resets all the camera properties to their default
        values. If a property has automation, the automatic will be enabled. WHAT DOES THIS MEAN??
        """
        return

    # TODO: perhaps a def reset() to disentangle from the stop()?

# if __name__ == "__main__":
# #    ImagingSourceDMKCamera()
#     camera_instance = Camera.camera_from_ini("ImagingSourceDMKCamera")


