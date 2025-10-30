"""Module tisgrabber of viperleed.gui.measure.camera.drivers.imagingsource.

This is an extension of the original tisgrabber.py module available
in the original form at:
github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/
          Python/Open%20Camera%2C%20Grab%20Image%20to%20OpenCV

Edits include code (and code-style) improvements.

The original tisgrabber.py module was:
Created on Mon Nov 21 09:44:40 2016
@author: Daniel Vassmer, Stefan_Geissler
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-18'
__license__ = 'GPLv3+'

# pylint: disable=too-many-lines
# Does not really make sense to split this module apart, especially
# considering that the WindowsCamera class already exceeds the 1000
# lines limit.

from enum import IntEnum
from ctypes import CFUNCTYPE
from ctypes import POINTER
from ctypes import Structure as CStructure
from ctypes import cast as c_cast
from ctypes import c_char
from ctypes import c_char_p
from ctypes import c_float
from ctypes import c_int
from ctypes import c_long
from ctypes import c_ubyte
from ctypes import c_ulong
from ctypes import c_void_p
from ctypes import py_object
from ctypes import windll
import os
from pathlib import Path
import sys

from viperleed.gui.measure.camera.drivers.imagingsource.winerrors import (
    DLLReturns,
    ImagingSourceError,
    check_dll_return,
    )
from viperleed.gui.measure.camera.drivers.imagingsource.properties import (
    CameraProperty,
    PropertyDiscoveryCallbackType,
    SwitchProperty,
    VCDPropertyInterface,
    VideoProperty,
    store_property,
    )
from viperleed.gui.measure.camera.drivers.imagingsource.models import ISModels
from viperleed.gui.measure.classes.settings import SystemSettings


c_float_p = POINTER(c_float)
c_long_p = POINTER(c_long)
c_int_p = POINTER(c_int)


def get_dll_path():
    """Return the path to the dll files."""
    drivers_path = SystemSettings().get('PATHS', 'drivers', fallback=None)
    if not drivers_path:
        raise ImportError
    try:
        dll_path = next(Path(drivers_path).rglob('TIS_UDSHL11*.dll')).parent
    except StopIteration:
        raise ImportError from None
    return dll_path


def _to_bytes(string):
    """Return string as bytes.

    Parameters
    ----------
    string : str, bytes, bytearray or None
        The string to be encoded

    Returns
    -------
    encoded_string : bytes or None
        UTF-8-encoded version of string. Returns None
        if string None.

    Raises
    ------
    TypeError
        If string is neither str, bytes or bytearray
    """
    if string is None or isinstance(string, bytes):
        return string
    if isinstance(string, bytearray):
        return bytes(bytearray)
    if isinstance(string, str):
        return string.encode()
    raise TypeError(f"Invalid string type {type(string).__name__!r}. "
                    "Expected, 'str', 'bytes', or 'bytearray")


class SinkFormat(IntEnum):
    """Available video formats for Imaging Source cameras."""

    Y800 = 0      # Monochrome 1-byte;    top-left
    RGB24 = 1     # 3 bytes: B G R;       top-left
    RGB32 = 2     # 4 bytes: B G R alpha; top-left; alpha unused == 0
    UYVY = 3      # 4-bytes color;        top-left
    Y16 = 4       # Monochrome 2-bytes;   top-left
    NONE = 5      # Used only as return value to signal 'invalid'
    MEGA = 65536  # Guard for max enum value, i.e., 2**sizeof(int)
    # RGB64       # 8 bytes: 2B 2G 2R 2A; top-left; seems unsupported

    @classmethod
    def get(cls, sink_format):
        """Return an appropriate enum from value or name.

        Parameters
        ----------
        sink_format : str, int or SinkFormat
            Some reference to a sink format.

        Returns
        -------
        sink_format : SinkFormat
            The appropriate attribute

        Raises
        ------
        ValueError
            If sink_format is not a valid sink format
        ImagingSourceError
            If sink_format would yield .NONE or .MEGA
        """
        # pylint: disable=redefined-variable-type
        # Easier not to have different names.
        if isinstance(sink_format, str):
            sink_format = getattr(cls, sink_format, None)
            if sink_format is None:
                raise ValueError("Unknown sink/video format "
                                 f"{sink_format}")
        else:
            try:
                sink_format = cls(sink_format)
            except ValueError as err:
                raise ValueError("Unknown sink/video format "
                                 f"{sink_format}") from err
        if sink_format == SinkFormat.NONE:
            raise ImagingSourceError("Invalid sink/video "
                                     f"format {sink_format.name}",
                                     err_code=DLLReturns.INVALID_SINK_FORMAT)
        if sink_format == SinkFormat.MEGA:
            raise ImagingSourceError("Camera should never return "
                                     "SinkFormat.MEGA",
                                     err_code=DLLReturns.INVALID_SINK_FORMAT)
        return sink_format

    @property
    def display_name(self):
        """Return a descriptive name for this format."""
        if self.n_colors == 1:
            return f"{8*self.n_bytes:d}-bit monochrome"
        _map = {SinkFormat.RGB24: "RGB color (8-bits)",
                SinkFormat.RGB32: "RGB+alpha color (8-bits)",
                SinkFormat.UYVY: "UYVY color"}
        try:
            return _map[self]
        except KeyError:
            raise ImagingSourceError(
                f"Format {self} does not have a display_name.",
                err_code=DLLReturns.INVALID_SINK_FORMAT
                ) from None

    @property
    def n_colors(self):
        """Return the number of color channels for a format."""
        if self.name.startswith('Y'):
            return 1
        if self is SinkFormat.RGB24:
            return 3
        if self is SinkFormat.RGB32:
            return 4
        raise ImagingSourceError(f"Format {self} does not have a "
                                 "defined number of color channels.",
                                 err_code=DLLReturns.INVALID_SINK_FORMAT)

    @property
    def n_bytes(self):
        """Return the number of bytes per pixel and per color channel."""
        # Always 1 byte, except Y16, and, when will be supported, RGB64
        if self in (SinkFormat.Y16,):
            return 2
        return 1

    @property
    def green_channel(self):
        """Return the 0-based index of the green channel."""
        if self.name.startswith('RGB'):
            return 1  # B, G, R
        if self.name.startswith('Y'):
            return 0  # monochrome
        raise ValueError(f"No green channel for {self.name} format.")


class StreamMode(IntEnum):
    """Enum for setting the stream mode of the camera.

    Attributes
    ----------
    CONTINUOUS
        The camera will pass each frame acquired to the frame-ready
        callback function
    SNAP
        Images are passed to the frame-ready callback only when
        snap_live_image() is called.
    """

    CONTINUOUS = 0
    SNAP = 1


# pylint: disable=too-few-public-methods
# Just an interface class to which we can make C pointers to.
class GrabberHandle(CStructure):
    """C-interface class for Imaging Source camera handle."""

    _fields_ = [('unused', c_int)]

    def __repr__(self):
        """Return string representation of self."""
        return f"{self.__class__.__name__}({self.unused})"
# pylint: enable=too-few-public-methods


# ctypes.POINTER returns a class
GrabberHandlePtr = POINTER(GrabberHandle)
GrabberHandlePtr.__repr__ = lambda self: f"GrabberHandlePtr({self.contents})"


# Prototype of disconnected callbacks (best to use as a decorator):
DisconnectedCallbackType = CFUNCTYPE(
    c_void_p,          # return type
    GrabberHandlePtr,  # HGRABBER
    py_object,         # python camera object
    )

DisconnectedCallbackType.__ctypeswrapper__ = 'DisconnectedCallbackType'


# Prototype of frame-ready callbacks (best to use as a decorator):
FrameReadyCallbackType = CFUNCTYPE(
    c_void_p,          # return type
    GrabberHandlePtr,  # HGRABBER
    POINTER(c_ubyte),  # pointer to first byte of image data
    c_ulong,           # current frame number
    py_object          # other data (python) passed along
    )

FrameReadyCallbackType.__ctypeswrapper__ = 'FrameReadyCallbackType'


# pylint: disable=too-many-public-methods,too-many-instance-attributes
# Bugs? I count 17/20 (not 35/20) methods and 7/7 (not 8/7) attributes
class WindowsCamera:
    """Interface class for Imaging Source camera in Windows."""

    # To prevent creating bindings to the dll functions each time
    # methods are called, each method that calls a dll function
    # has ctypes definition of the function used right above it.
    # Note: ctypes assumes restype == c_int. Thus restype is
    # omitted below where this is the case.

    dll_path = get_dll_path()
    if sys.maxsize > 2**32:
        try:
            windll.LoadLibrary(str(dll_path / "TIS_UDSHL11_x64.dll"))
        except FileNotFoundError:
            raise ImportError
        _dll = windll.LoadLibrary(str(dll_path / "tisgrabber_x64.dll"))
    else:
        try:
            windll.LoadLibrary(str(dll_path / "TIS_UDSHL11.dll"))
        except FileNotFoundError:
            raise ImportError
        _dll = windll.LoadLibrary(str(dll_path / "tisgrabber.dll"))

    __initalized = False

    def __init__(self):
        """Initialize grabber instance."""
        # Call init_library with None straight away. It was done the
        # same way in the original tisgrabber.py. init_library will
        # do something only the first time an instance is created.
        self.__init_library()
        self.__handle = self.create_grabber()
        self.__vcd_properties = {}
        self._has_disconnected_callback = False
        self._has_frame_ready_callback = False
        # Store video format info, as this is static for an open
        # device, and takes a while to retrieve the first time.
        self.__video_fmt_info = {'min_w': None, 'min_h': None,
                                 'max_w': None, 'max_h': None,
                                 'd_w': None, 'd_h': None,
                                 'types': None}
        self.__model = ''

    @classmethod
    def __init_library(cls, license_key=None):
        """Initialize the IC Imaging Control DLL library.

        This function must be called at least once before any
        other functions of the DLL are called.

        Parameters
        ----------
        license_key : str, bytes or None, optional
            The license key. Pass None if using a 'trial' version.
            Default is None.

        Returns
        -------
        None.
        """
        if not cls.__initalized:
            _dll_init = cls._dll.IC_InitLibrary
            _dll_init.argtypes = (c_char_p,)
            _dll_init.errcheck = check_dll_return()

            _dll_init(_to_bytes(license_key))
            cls.__initalized = True

    _dll_create_grabber = _dll.IC_CreateGrabber
    _dll_create_grabber.restype = GrabberHandlePtr
    _dll_create_grabber.argtypes = None
    _dll_create_grabber.errcheck = check_dll_return('pointer')

    def create_grabber(self):
        """Return a new grabber handle.

        Unused grabber handles should be released by calling
        IC_ReleaseGrabber(). There is no need to release a handle
        when closing a device, as handles can be reused.

        Returns
        -------
        handle : GrabberHandlePtr
            New handle. Can be used to open devices.
        """
        return self._dll_create_grabber()

    @property
    def devices(self):
        """Return the unique names of available devices.

        Returns
        -------
        unique_names : list
            List of unique device names, each as a string.
        """
        # This is not great, as the loop is done in python and requires
        # multiple interactions with the camera, but IC_ListDevices does
        # not return unique names (i.e., no serial number info).
        names = (self.__unique_name_from_index(i)
                 for i in range(self.__n_devices))
        return [n for n in names if n]

    @property
    def dynamic_range(self):
        """Return the dynamic range in bits for the open camera."""
        return ISModels[self.__model].dynamic_range

    @property
    def exposure_range(self):
        """Return min and max exposure times in milliseconds."""
        min_exp, max_exp = self.get_vcd_property_range("Exposure", "Value",
                                                       method=int)
        return min_exp/1000, max_exp/1000

    @property
    def exposure(self):
        """Return the exposure time in milliseconds."""
        # Get the exposure in microseconds as an integer. The float
        # method would rather return seconds, and it is less accurate.
        int_exposure = self.get_vcd_property("Exposure", "Value", method=int)

        return int_exposure/1000  # milliseconds

    @exposure.setter
    def exposure(self, exposure_time):
        """Set the exposure time in milliseconds."""
        exp_min, exp_max = self.exposure_range
        if not exp_min <= exposure_time <= exp_max:
            raise ValueError("Cannot set exposure time to "
                             f"{exposure_time} ms: Value out of range.")
        exposure_time = round(exposure_time*1000)  # microseconds
        self.set_vcd_property("Exposure", "Value", exposure_time, method=int)

    @property
    def frame_rate_range(self):     # TODO: slow. Store property? When does it change? Probably only with video/sink format
        """Return the minimum and maximum available image readout rates."""
        # This is a bit of a workaround: IC_GetAvailableFrameRates
        # seems not to return the correct range. Will instead set
        # unreasonably small and unreasonably large frame rates,
        # and see what the camera actually sets internally.
        try:
            old_frame_rate = self.frame_rate
        except ImagingSourceError:
            old_frame_rate = 1024

        self.frame_rate = 0
        min_rate = self.frame_rate

        self.frame_rate = 1024
        max_rate = self.frame_rate

        self.frame_rate = old_frame_rate

        return min_rate, max_rate

    _dll_get_frame_rate =  _dll.IC_GetFrameRate
    _dll_get_frame_rate.restype = c_float
    _dll_get_frame_rate.argtypes = (GrabberHandlePtr,)
    _dll_get_frame_rate.errcheck = check_dll_return('>0')

    @property
    def frame_rate(self):
        """Return the rate at which images are read from the device.

        This corresponds to the rate at which images are transmitted
        back from the camera through the GigE communication line.
        However, mind that frames may get lost if images are sent
        back too quickly, typically due to non-ideal behavior of the
        Ethernet hardware.

        Returns
        -------
        frame_rate : float
            A positive number if image readout rate is supported,
            zero otherwise.
        """
        return self._dll_get_frame_rate(self.__handle)

    _dll_set_frame_rate = _dll.IC_SetFrameRate
    _dll_set_frame_rate.argtypes = GrabberHandlePtr, c_float
    _dll_set_frame_rate.errcheck = check_dll_return(
        exclude_errors=('NOT_IN_LIVE_MODE', 'NO_PROPERTYSET',
                        'DEFAULT_WINDOW_SIZE_SET', 'WRONG_XML_FORMAT',
                        'WRONG_INCOMPATIBLE_XML',
                        'NOT_ALL_PROPERTIES_RESTORED', 'DEVICE_NOT_FOUND')
        )

    @frame_rate.setter
    def frame_rate(self, rate):
        """Set a new rate for reading out images.

        This corresponds to the rate at which images are transmitted
        back from the camera through the GigE communication line.
        However, mind that frames may get lost if images are sent
        back too quickly, typically due to non-ideal behavior of the
        Ethernet hardware.

        Parameters
        ----------
        rate : float
            The new frame rate. If rate is larger than the maximum
            allowed rate (i.e., the ideal transmission rate for GigE
            at the current image size), the camera will rather use
            this maximum (theoretical) frame rate.
        """
        was_running = self.is_running
        if was_running:
            self.pause()
        self._dll_set_frame_rate(self.__handle, rate)
        if was_running:
            self.start()

    @property
    def gain_range(self):
        """Return min and max available gains in decibel."""
        min_gain, max_gain = self.get_vcd_property_range("Gain", "Value",
                                                         method=float)
        return min_gain, max_gain

    @property
    def gain(self):
        """Return the current gain used in decibel."""
        return self.get_vcd_property("Gain", "Value", method=float)

    @gain.setter
    def gain(self, new_gain):
        """Set a new gain in decibel."""
        self.set_vcd_property("Gain", "Value", new_gain, method=float)

    _dll_get_image_description = _dll.IC_GetImageDescription
    _dll_get_image_description.errcheck = check_dll_return()
    _dll_get_image_description.argtypes = (
        GrabberHandlePtr, c_long_p, c_long_p, c_int_p, c_int_p
        )

    @property
    def image_info(self):
        """Return properties of the current video format and sink type.

        Returns
        -------
        width, height : int
            Width and height of the image in pixels
        bytes_per_pixel : int
            Number of bytes for each pixel and color channel
        n_colors : int
            Number of color channels
        format : SinkFormat
            The current color format. Notice that this is the color
            format of the image stored in the internal memory, NOT
            the video format.

        Raises
        ---------
        ImagingSourceError
            If the format returned is unknown.
        """
        width, height = c_long(), c_long()
        bits_per_pixel, color_format_c = c_int(), c_int()

        try:
            self._dll_get_image_description(self.__handle, width, height,
                                            bits_per_pixel, color_format_c)
        except ImagingSourceError as err:
            # No info available. May happen when the sink format
            # is not up to date in the camera. Call the getter,
            # which updates it internally, then try again.
            if err.err_code is DLLReturns.ERROR:
                _ = self.sink_format
                self._dll_get_image_description(self.__handle, width, height,
                                                bits_per_pixel, color_format_c)
            else:
                raise

        try:
            color_format = SinkFormat.get(color_format_c.value)
        except ValueError as err:
            raise ImagingSourceError(
                err.args[0] + " received.",
                err_code=DLLReturns.INVALID_SINK_FORMAT
                ) from err

        n_colors = color_format.n_colors
        bytes_per_pixel = (bits_per_pixel.value / n_colors) // 8

        return (width.value, height.value, int(bytes_per_pixel),
                color_format.n_colors, color_format)

    _dll_is_running = _dll.IC_IsLive
    _dll_is_running.argtypes = (GrabberHandlePtr,)
    _dll_is_running.errcheck = check_dll_return(
        ">=0",
        exclude_errors=("WRONG_XML_FORMAT", "CAMERA_PROPERTY_NOT_AVAILABLE")
        )

    @property
    def is_running(self):
        """Return whether the camera is currently running."""
        return bool(self._dll_is_running(self.__handle))

    _dll_get_sink_format = _dll.IC_GetFormat
    _dll_get_sink_format.argtypes = (GrabberHandlePtr,)

    @property
    def sink_format(self):
        """Return the color format of the sink type currently set.

        The sink type describes the format of the buffer where images
        are stored internally in the camera. The method returns a valid
        value only after IC_PreprareLive() or start() are called.
        Notice that the sink type may differ from the video format.

        Returns
        -------
        color_format : SinkFormat
            The sink format currently in use.

        Raises
        ------
        ImagingSourceError
            If the value returned is not a valid SinkFormat
        """
        sink_format_c = self._dll_get_sink_format(self.__handle)
        try:
            sink_format = SinkFormat.get(sink_format_c)
        except ValueError as err:
            raise ImagingSourceError(
                err.args[0] + " received",
                err_code=DLLReturns.INVALID_SINK_FORMAT
                ) from err
        except ImagingSourceError:
            if sink_format_c == SinkFormat.NONE.value:
                # Camera was not yet started: sink format is unavailable.
                # Start/stop, then try again.
                self.start()
                self.stop()

                sink_format_c = self._dll_get_sink_format(self.__handle)
                try:
                    sink_format = SinkFormat.get(sink_format_c)
                except ValueError as err:
                    raise ImagingSourceError(
                        err.args[0] + " received",
                        err_code=DLLReturns.INVALID_SINK_FORMAT
                        ) from err
                except ImagingSourceError as err:
                    raise ImagingSourceError(
                        "Could not retrieve sink format.",
                        err_code=DLLReturns.INVALID_SINK_FORMAT
                        ) from err
            else:  # Mega
                raise
        return sink_format

    _dll_set_sink_format = _dll.IC_SetFormat
    _dll_set_sink_format.argtypes = GrabberHandlePtr, c_int
    _dll_set_sink_format.errcheck = check_dll_return()

    @sink_format.setter
    def sink_format(self, sink_format):
        """Set the sink type.

        This property must be set once before images can be snapped.
        The sink type describes the format of the buffer where
        snapped images are stored. Notice that the sink type may
        differ from the video format.

        Parameters
        ----------
        sink_format : str or int or SinkFormat
            One of the values in SinkFormat. Note that UYVY can
            only be used in conjunction with a UYVY video format
        """
        was_triggered = self.trigger_enabled
        was_running = self.is_running
        if was_running:
            self.stop()
        sink_format = SinkFormat.get(sink_format)
        self._dll_set_sink_format(self.__handle, sink_format.value)

        # Prepare camera to allow triggering. When the sink format
        # has been changed, one needs to snap an image. This is done
        # automatically when __set_frame_acquisition_mode is set to SNAP
        self.__set_frame_acquisition_mode(StreamMode.SNAP)

        # Reset stream mode to continuous as then it is just a
        # matter of having enable_trigger(True/False) for getting
        # single frames delivered upon request [start('triggered')
        # + send_software_trigger()] or a continuous stream of
        # images [start('continuous')].
        self.trigger_enabled = was_triggered
        self.__set_frame_acquisition_mode(StreamMode.CONTINUOUS)
        if was_running:
            self.start()

    @property
    def trigger_burst_count_range(self):
        """Return the max and min number of frames in a trigger burst."""
        return self.get_vcd_property_range("Trigger", "Burst Count")

    @property
    def trigger_burst_count(self):
        """Return the number of frames in a trigger burst."""
        return self.get_vcd_property("Trigger", "Burst Count")

    @trigger_burst_count.setter
    def trigger_burst_count(self, n_frames):
        """Set the number of frames in a trigger burst.

        Notice that when giving a trigger burst, the time to receive
        each frame is (a bit more than) the maximum between .exposure
        and 1/actual_frame_rate, where the actual frame rate may be
        much smaller than .frame_rate in case frames are lost in
        communication. This is commonly the case when large images
        with short exposure times are being acquired.

        Parameters
        ----------
        n_frames : int
            Number of frames that should be used in a trigger
            burst. Notice that if trigger_burst_count is
            changed after send_software_trigger() and before
            all frames are acquired, the updated number of
            frames will be honored. If n_frames is larger than
            the maximum allowed (typically 1000), the maximum
            will be set.
        """
        self.set_vcd_property("Trigger", "Burst Count", n_frames)

    @property
    def trigger_enabled(self):
        """Return whether triggering is currently enabled."""
        return self.get_vcd_property('Trigger', 'Enable') is SwitchProperty.ON

    @trigger_enabled.setter
    def trigger_enabled(self, enabled):
        """Enable or disable capability to respond to a trigger signal.

        Trigger signals may either be via hardware or via software.
        Notice that disabling the trigger while the camera was
        already .start()ed makes it return a stream of frames.

        Parameters
        ----------
        enabled : bool
            If True, the device will afterwards be responsive to a
            (hardware or software) trigger signal.
        """
        self.set_vcd_property('Trigger', 'Enable', bool(enabled))

    @property
    def vcd_properties(self):
        """Return a dict of VCD properties/elements/interfaces.

        Returns
        -------
        vcd_properties : dict
            The dict has the following structure:
            {prop_name_1:
                {el_name_1_1 : (VCDPropertyInterface, vcd_properties, ...),
                 ...},
             prop_name_2:
                {...},
             ...}
        """
        return self.__vcd_properties

    @property
    def video_format_shape_range(self):
        """Return min_width, min_height, max_width, max_height."""
        if not self.__video_fmt_info['min_w']:
            widths = []
            heights = []
            types = set()
            for vid_fmt in self.__default_video_formats:
                # formats are of the form "<type> (wxh)"
                fmt_type, vid_fmt = vid_fmt.split()
                vid_fmt = vid_fmt.replace('(','').replace(')','')
                width, height = vid_fmt.split('x')
                widths.append(int(width))
                heights.append(int(height))
                types.add(fmt_type)
            self.__video_fmt_info['min_w'] = min(widths)
            self.__video_fmt_info['min_h'] = min(heights)
            self.__video_fmt_info['max_w'] = max(widths)
            self.__video_fmt_info['max_h'] = max(heights)
            self.__video_fmt_info['types'] = sorted(
                SinkFormat.get(t) for t in types
                )
        return (self.__video_fmt_info['min_w'], self.__video_fmt_info['min_h'],
                self.__video_fmt_info['max_w'], self.__video_fmt_info['max_h'])

    @property
    def video_format_shape_increments(self):
        """Return the minimum increments allowed for width and height."""
        if not self.__video_fmt_info['d_w']:
            min_w, min_h, max_w, max_h = self.video_format_shape_range
            old_w, old_h, *_, old_format = self.image_info
            fmt = "{} ({}x{})"
            incr_w = incr_h = -1

            for width in range(min_w + 1, max_w):
                # setting an unavailable format raises errors
                try:
                    self._dll_set_video_format(
                        self.__handle,
                        _to_bytes(fmt.format('Y16', width, min_h))
                        )
                except ImagingSourceError:
                    pass
                else:
                    incr_w = width - min_w
                    break
            for height in range(min_h + 1, max_h):
                # setting an unavailable format raises errors
                try:
                    self._dll_set_video_format(
                        self.__handle,
                        _to_bytes(fmt.format('Y16', min_w, height))
                        )
                except ImagingSourceError:
                    pass
                else:
                    incr_h = height - min_h
                    break
            self.set_video_format(fmt.format(old_format.name, old_w, old_h))
            self.__video_fmt_info['d_w'] = incr_w
            self.__video_fmt_info['d_h'] = incr_h
        return (self.__video_fmt_info['d_w'], self.__video_fmt_info['d_h'])

    _dll_video_format_width = _dll.IC_GetVideoFormatWidth
    _dll_video_format_height = _dll.IC_GetVideoFormatHeight
    _dll_video_format_width.argtypes = (GrabberHandlePtr,)
    _dll_video_format_height.argtypes = (GrabberHandlePtr,)
    _dll_video_format_width.errcheck = check_dll_return('>0')
    _dll_video_format_height.errcheck = check_dll_return('>0')

    @property
    def video_format_shape(self):
        """Return (width, height) (pixels) of the current video format."""
        width = self._dll_video_format_width(self.__handle)
        height = self._dll_video_format_height(self.__handle)
        return width, height

    @property
    def video_formats_available(self):
        """Return a list of available SinkFormat for this camera."""
        if self.__video_fmt_info['types'] is None:
            # No video_format_shape_range was never called.
            _ = self.video_format_shape_range
        return  self.__video_fmt_info['types']

    def abort_trigger_burst(self):
        """Abort frame delivery."""
        if not self.is_running or not self.trigger_enabled:
            return
        self.pause()
        self._dll_start_live(self.__handle, False)

    _dll_close_device = _dll.IC_CloseVideoCaptureDevice
    _dll_close_device.restype = None  # Returns void
    _dll_close_device.argtypes = (GrabberHandlePtr,)

    def close(self):
        """Close the currently open video-capture device."""
        self._dll_close_device(self.__handle)
        self.__vcd_properties = {}
        # Would be nicer to have it set by the device lost callback,
        # as this does not prevent one from closing/reopening and
        # trying to set a new callback. It is unclear whether the
        # problem with setting the callback twice is only there
        # when running from the terminal.
        self._has_disconnected_callback = False
        self._has_frame_ready_callback = False

        self.__video_fmt_info = {'min_w': None, 'min_h': None,
                                 'max_w': None, 'max_h': None,
                                 'd_w': None, 'd_h': None}

    _dll_enable_auto_camera_prop = _dll.IC_EnableAutoCameraProperty
    _dll_enable_auto_camera_prop.argtypes = GrabberHandlePtr, c_int, c_int
    _dll_enable_auto_camera_prop.errcheck = check_dll_return()
    _dll_get_auto_camera_prop = _dll.IC_GetAutoCameraProperty
    _dll_get_auto_camera_prop.argtypes = GrabberHandlePtr, c_int, c_int_p
    _dll_get_auto_camera_prop.errcheck = check_dll_return(
        valid_returns=("CAMERA_PROPERTY_NOT_AVAILABLE",)
        )
    _dll_enable_auto_video_prop = _dll.IC_EnableAutoVideoProperty
    _dll_enable_auto_video_prop.argtypes = GrabberHandlePtr, c_int, c_int
    _dll_enable_auto_video_prop.errcheck = check_dll_return()

    def enable_auto_properties(self, enabled):
        """Enable/disable automatic properties.

        Properties affected:
        Exposure, Gain, all available CameraProperty,
        Auto-centering of 'partial scan' (i.e., ROI)

        Parameters
        ----------
        enabled : bool
            Enable if True, disable otherwise.

        Raises
        -------
        ImagingSourceError
            If the method fails to enable/disable some
            of the available automatic properties.
        """
        enabled = bool(enabled)
        check_enabled = c_int()
        failed = []

        # Only GAIN can be set to 'Auto' in the Video properties
        self._dll_enable_auto_video_prop(self.__handle,
                                         VideoProperty.GAIN.value,
                                         enabled)
        for cam_prop in CameraProperty:
            if cam_prop is CameraProperty.MEGA:
                continue
            available = (self._dll_get_auto_camera_prop(self.__handle,
                                                        cam_prop.value,
                                                        check_enabled)
                         != DLLReturns.CAMERA_PROPERTY_NOT_AVAILABLE.value)
            if not available:
                continue
            try:
                self._dll_enable_auto_camera_prop(self.__handle,
                                                  cam_prop.value,
                                                  enabled)
            except ImagingSourceError:
                failed.append(cam_prop)

        for vcd_prop in ("Gain", "Exposure"):
            self.set_vcd_property(vcd_prop, 'Auto', enabled)
        self.set_vcd_property("Partial scan", "Auto-center", enabled)

        if failed:
            raise ImagingSourceError(
                f"Failed to {'en' if enabled else 'dis'}able: {failed}",
                err_code=DLLReturns.FAILED_TO_ENABLE
                )

    # TODO: check which errors
    _dll_open_by_unique_name = _dll.IC_OpenDevByUniqueName
    _dll_open_by_unique_name.argtypes = GrabberHandlePtr, c_char_p
    _dll_open_by_unique_name.errcheck = check_dll_return()

    def open(self, unique_name):
        """Open a video capture by using its unique name.

        Use IC_GetUniqueName() to retrieve the unique name of an
        open device. __unique_name_from_index(index) can be used
        before the device is open to retrieve the unique name.

        The available VCD properties of the device are also stored,
        and can be retrieved with .vcd_properties. The newly open
        device is set to deliver frames continuously to the frame-
        ready callback.

        Parameters
        ----------
        unique_name : bytes
            Unique name of the device to be opened

        Raises
        ------
        RuntimeError
            If the camera model contained in unique_name is
            unknown. Requires adding the model (and its info)
            to the models.py module.
        """
        self._dll_open_by_unique_name(self.__handle,
                                      _to_bytes(unique_name))
        self.__model = '_'.join(unique_name.split()[:2]).replace('.', '_')
        try:
            ISModels[self.__model]
        except AttributeError:
            self.close()
            raise RuntimeError(
                "Could not find Imaging Source camera model "
                f"{self.__model} in the known model list. "
                "Please update models.py accordingly."
                ) from None

        # Do some of the slow stuff already here, such that
        # following calls will be fast.
        self.__find_vcd_properties()
        _ = self.video_format_shape_increments

        # Set stream mode always to continuous as then it is just a
        # matter of having enable_trigger(True/False) for getting
        # single frames delivered upon request [start('triggered')
        # + send_software_trigger()] or a continuous stream of
        # images [start('continuous')].
        self.__set_frame_acquisition_mode(StreamMode.CONTINUOUS)

    _dll_pause_live = _dll.IC_SuspendLive
    _dll_pause_live.errcheck = check_dll_return()
    _dll_pause_live.argtypes = (GrabberHandlePtr,)

    def pause(self):
        """Pause camera.

        A puased camera can be .start()ed a bit more quickly than
        a .stop()ped one. However, one should not attempt to edit
        .video_format or .sink_format.

        Returns
        -------
        None.
        """
        self._dll_pause_live(self.__handle)

    _dll_set_callback = _dll.IC_SetCallbacks
    _dll_set_callback.errcheck = check_dll_return()
    # The py_objects at the end are the same that will later be passed
    # on to the callbacks as the last argument.
    _dll_set_callback.argtypes = (
        GrabberHandlePtr,
        FrameReadyCallbackType, py_object, # FRAME_READY_CALLBACK
        DisconnectedCallbackType, py_object, # DEVICE_LOST_CALLBACK
        )

    def set_callbacks(self, on_frame_ready, on_disconnected,
                      py_obj_frame_ready, camera):
        """Define callback functions.

        Parameters
        ----------
        on_frame_ready : callable
            Should be decorated with @FrameReadyCallbackType to ensure
            proper typing. This is the function that will be called
            whenever a new frame has been acquired (in 'continuous'
            mode) or whenever an image is 'snapped'. The callable
            should return None, and have signature:
            on_frame_ready(__handle, img_buffer, frame_no, py_obj_frame_ready)
            where:
                __handle : int
                    Pointer to grabber handle. Should not be used.
                img_buffer : ctypes.POINTER(c_ubyte)
                    Pointer to first byte of bottom line of image data
                frame_no : int
                    Frame number. Unclear what it corresponds to.
                py_obj_frame_ready : object
                    Same object as below
        on_disconnected : callable
            Should be decorated with @DisconnectedCallbackType to ensure
            proper typing. This is the function that will be called
            whenever the camera is disconnected. The callable
            should return None, and have signature:
            on_disconnected_(__handle, camera)
            where:
                __handle : int
                    Pointer to grabber handle. Should not be used.
                camera : ImagingSourceCamera
                    The camera python object.
        py_obj_frame_ready : object
            The python object that is passed as the last argument to
            on_frame_ready. If the callback is to modify this, it is
            necessary to pass a dictionary-like object (or a class
            instance with attributes).
        camera : ImagingSourceCamera
            The python camera object that is passed as the last argument
            to on_disconnected.

        Raises
        -------
        TypeError
            If one of the callbacks is not callable.
        ValueError
            If one of the callbacks was not properly decorated.
        ImagingSourceError
            If this method is called more than once
            before rebooting a camera.
        """
        # Make sure the callbacks have been wrapped correctly.
        allowed = ('FrameReadyCallbackType', 'DisconnectedCallbackType')
        for cb in (on_frame_ready, on_disconnected):
            if not hasattr(cb, '__call__'):
                raise TypeError(f'Callback {cb.__name__} '
                                'is not a callable.')
            if (not hasattr(cb, '__ctypeswrapper__')
                    or cb.__ctypeswrapper__ not in allowed):
                raise ValueError(f'Callback {cb.__name__} was not '
                                 'decorated with an allowed callback '
                                 'type. This is needed to ensure '
                                 'appropriate type checking/conversions')
        if (not self._has_disconnected_callback
                and not self._has_frame_ready_callback):
            self._dll_set_callback(self.__handle, on_frame_ready,
                                   py_obj_frame_ready, on_disconnected, camera)
            self._has_disconnected_callback = True
            self._has_frame_ready_callback = True
        else:
            raise ImagingSourceError(
                'Cannot set callbacks twice due to some bug in the '
                'underlying DLL. Power-down the camera, and retry.',
                err_code=DLLReturns.REQUIRES_REBOOT
                )

    _dll_set_frame_ready_callback = _dll.IC_SetFrameReadyCallback
    _dll_set_frame_ready_callback.errcheck = check_dll_return()
    # The py_object at the end is the same that will later be passed
    # on to the frame-ready callback as the last argument.
    _dll_set_frame_ready_callback.argtypes = (
        GrabberHandlePtr, FrameReadyCallbackType, py_object
        )

    def set_frame_ready_callback(self, on_frame_ready, py_obj_for_callback):
        """Define the frame-ready callback function.

        Parameters
        ----------
        on_frame_ready : callable
            Should be decorated with @FrameReadyCallbackType to ensure
            proper typing. This is the function that will be called
            whenever a new frame has been acquired (in 'continuous'
            mode) or whenever an image is 'snapped'. The callable
            should return None, and have signature:
            on_frame_ready(__handle, img_buffer, frame_no, py_obj_for_callback)
            where:
                __handle : int
                    Pointer to grabber handle. Should not be used.
                img_buffer : ctypes.POINTER(c_ubyte)
                    Pointer to first byte of bottom line of image data
                frame_no : int
                    Frame number. Unclear what it corresponds to.
                py_obj_for_callback : object
                    Same object as below
        py_obj_for_callback : object
            The python object that is passed as the last argument to
            on_frame_ready. If the callback is to modify this, it is
            necessary to pass a dictionary-like object (or a class
            instance with attributes).

        Raises
        -------
        TypeError
            If on_frame_ready is not callable
        ValueError
            If on_frame_ready was not decorated with
            @FrameReadyCallbackType
        ImagingSourceError
            If this method is called more than once
            before rebooting a camera.
        """
        raise NotImplementedError('This function should no longer be used.'
                                  'Use set_callbacks instead.')
        # Make sure the frame-ready callback has been wrapped correctly
        if not hasattr(on_frame_ready, '__call__'):
            raise TypeError(f"Frame-ready callback {on_frame_ready.__name__} "
                            "is not a callable")
        if (not hasattr(on_frame_ready, '__ctypeswrapper__')
                or on_frame_ready.__ctypeswrapper__ != 'FrameReadyCallbackType'):
            raise ValueError("Frame-ready callback was not decorated with "
                             "@FrameReadyCallbackType. This is needed to "
                             "ensure appropriate type checking/conversions")
        if not self._has_disconnected_callback:
            self._dll_set_frame_ready_callback(self.__handle, on_frame_ready,
                                               py_obj_for_callback)
            self._has_disconnected_callback = True
        else:
            raise ImagingSourceError(
                "Cannot set twice a callback due to some bug in the "
                "underlying DLL. Power-down the camera, and retry.",
                err_code=DLLReturns.REQUIRES_REBOOT
                )

    _dll_start_live = _dll.IC_StartLive
    _dll_start_live.argtypes = GrabberHandlePtr, c_int
    _dll_start_live.errcheck = check_dll_return()

    def start(self, mode=None, show_video=False):
        """Start camera.

        After starting, one has to wait a little time before frames
        can be retrieved. A similar wait is necessary after triggering
        and before calling stop(), otherwise frames are not delivered.

        When .trigger_enabled == False, it takes 1/frame_rate
        to deliver a frame; when .trigger_enabled == True, one needs
        trigger_delay + exposure + 1/frame_rate_max before a
        frame is delivered. A good minimum time for triggering
        is (exposure + 1/frame_rate_max)*1.05 [less than this
        will miss some frames].

        Parameters
        ----------
        mode : None or str, optional
            Which mode should the camera be started in. 'triggered'
            automatically enables triggering, 'continuous' (or any
            other value) will stream frames every 1/frame_rate seconds.
            If None, the mode will be derived from whether triggering
            is enabled or not. Default is None.
        show_video : bool, optional
            Do not show a video if False, show one if True. Frames
            will be delivered in any case (i.e., callbacks can be
            used). The value of this argument is used only when
            mode == 'continuous'; no video is ever shown in
            'triggered' mode. Default is False.
        """
        if mode is None:
            mode = 'triggered' if self.trigger_enabled else 'continuous'
        if mode == 'triggered':
            show_video = False
            self.trigger_enabled = True
        else:
            self.trigger_enabled = False
        if self.is_running:
            self.stop()             # could we use pause()?
        self._dll_start_live(self.__handle, bool(show_video))

    _dll_stop_live = _dll.IC_StopLive
    _dll_stop_live.argtypes = (GrabberHandlePtr,)

    def stop(self):
        """Stop the camera."""
        self._dll_stop_live(self.__handle)

    # #####################    VCD PROPERTIES    ######################
    #
    # The way properties are 'set', 'gotten' or 'range-gotten' is
    # always derived from their native 'interface', unless there
    # are multiple options, in which case the <method> kwarg is
    # used.

    __vcd_prop_checker = check_dll_return(
        include_errors=(
            'NO_HANDLE', 'NO_DEVICE', 'PROPERTY_ITEM_NOT_AVAILABLE',
            'PROPERTY_ELEMENT_NOT_AVAILABLE',
            'PROPERTY_ELEMENT_WRONG_INTERFACE'
            )
        )
    # Value as float: getter, range getter, setter.
    _dll_get_property_value_float = _dll.IC_GetPropertyAbsoluteValue
    _dll_get_property_value_float.errcheck = __vcd_prop_checker
    _dll_get_property_value_float.argtypes = (
        GrabberHandlePtr, c_char_p,c_char_p, c_float_p
        )
    _dll_get_property_value_range_float = _dll.IC_GetPropertyAbsoluteValueRange
    _dll_get_property_value_range_float.errcheck = __vcd_prop_checker
    _dll_get_property_value_range_float.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_float_p, c_float_p
        )
    _dll_set_property_value_float = _dll.IC_SetPropertyAbsoluteValue
    _dll_set_property_value_float.errcheck = __vcd_prop_checker
    _dll_set_property_value_float.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_float
        )

    # Value as int: getter, range getter, setter.
    _dll_get_property_value_int = _dll.IC_GetPropertyValue
    _dll_get_property_value_int.errcheck = __vcd_prop_checker
    _dll_get_property_value_int.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_long_p
        )
    _dll_get_property_value_range_int = _dll.IC_GetPropertyValueRange
    _dll_get_property_value_range_int.errcheck = __vcd_prop_checker
    _dll_get_property_value_range_int.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_long_p, c_long_p
        )
    _dll_set_property_value_int = _dll.IC_SetPropertyValue
    _dll_set_property_value_int.errcheck = __vcd_prop_checker
    _dll_set_property_value_int.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_int
        )

    # On/Off (== bool): getter, setter. [NO RANGE OBVIOUSLY]
    _dll_get_on_off_property = _dll.IC_GetPropertySwitch
    _dll_get_on_off_property.errcheck = __vcd_prop_checker
    _dll_get_on_off_property.argtypes= (
        GrabberHandlePtr, c_char_p, c_char_p, c_long_p
        )
    _dll_set_on_off_property = _dll.IC_SetPropertySwitch
    _dll_set_on_off_property.errcheck = __vcd_prop_checker
    _dll_set_on_off_property.argtypes= (
        GrabberHandlePtr, c_char_p, c_char_p, c_int
        )

    def get_vcd_property(self, property_name, element_name=None, method=None):
        """Return the current value of a property/element pair.

        Parameters
        ----------
        property_name : str or bytes
            The name of the property whose value is to be returned
        element_name : str, bytes or None, optional
            The name of the setting to be accessed (examples: 'Value',
            'Auto'). This argument is mandatory if the property has
            multiple elements. Otherwise, if not given, the only
            element present will be returned.
        method :  str, VCDPropertyInterface
                 or VCDPropertyInterface.value, optional
            Which method should be used to access the property/element.
            This is honored only in case the property/element has
            multiple access methods, in which case it is mandatory.

        Returns
        -------
        property_value : int, float or SwitchProperty
            The value of the property/element. The return type
            depends on the choice of method.

        Raises
        ------
        ImagingSourceError
            If errors occur in retrieving properties.
        NotImplementedError
            If the property required is only accessible via the
            MAPSTRINGS interface method. Currently unsupported.
        ValueError
            If method is not one of the interface methods known.
        """
        (prop_name,
         elem_name,
         method) = self.__get_property_element_interface(property_name,
                                                         element_name,
                                                         method)
        if method is VCDPropertyInterface.ABSOLUTEVALUE:
            getter = self._dll_get_property_value_float
            property_value = c_float()
        elif method is VCDPropertyInterface.RANGE:
            getter = self._dll_get_property_value_int
            property_value = c_long()
        elif method is VCDPropertyInterface.SWITCH:
            getter = self._dll_get_on_off_property
            property_value = c_long()
        elif method is VCDPropertyInterface.MAPSTRINGS:
            raise NotImplementedError("Map strings not yet implemented")
        else:
            raise ValueError(
                f"Unexpected method {method} for getting VCD "
                f"property/element {property_name}/{element_name}"
                )
        getter(self.__handle, prop_name, elem_name, property_value)

        if method is not VCDPropertyInterface.SWITCH:
            return property_value.value
        return SwitchProperty.get(property_value.value)

    def set_vcd_property(self, property_name, element_name,
                         value, method=None):
        """Set the value of a VCD property/element pair.

        Parameters
        ----------
        property_name : str or bytes
            The name of the property whose value is to be changed
        element_name : str or bytes
            The name of the setting to be accessed (e.g. 'Value',
            'Auto').
        value : object
            The value to be set for the property/element. For on/off
            settings (method == VCDPropertyInterface.SWITCH) accepts
            strings 'on'/'off' (not case sensitive), SwitchProperty,
            or any object with a truth value. Otherwise it should be
            a number, or a string that evaluates to a number.
        method : str, VCDPropertyInterface
                 or VCDPropertyInterface.value, optional
            Which method should be used to access the property/element.
            This is honored only in case the property/element has
            multiple access methods, in which case it is mandatory.

        Raises
        ------
        ImagingSourceError
            If errors occur in setting properties.
        NotImplementedError
            If the property required is only accessible via the
            MAPSTRINGS interface method. Currently unsupported.
        ValueError
            If method is not one of the interface methods known.
        """
        (prop_name,
         elem_name,
         method) = self.__get_property_element_interface(property_name,
                                                         element_name,
                                                         method)

        if method is VCDPropertyInterface.ABSOLUTEVALUE:
            setter = self._dll_set_property_value_float
            value = float(value)
        elif method is VCDPropertyInterface.RANGE:
            setter = self._dll_set_property_value_int
            value = int(value)
        elif method is VCDPropertyInterface.SWITCH:
            setter = self._dll_set_on_off_property
            value = SwitchProperty.get(value).value
        elif method is VCDPropertyInterface.MAPSTRINGS:
            raise NotImplementedError("Map strings not yet implemented")
        else:
            raise ValueError(
                f"Unexpected method {method} for setting VCD "
                f"property/element {property_name}/{element_name}"
                )
        setter(self.__handle, prop_name, elem_name, value)

    def get_vcd_property_range(self, property_name, element_name=None,
                               method=None):
        """Return minimum and maximum values for a property/element pair.

        Parameters
        ----------
        property_name : str or bytes
            The name of the property whose range should be returned
        element_name : str, bytes or None, optional
            The name of the setting to be accessed (examples: 'Value',
            'Auto'). This argument is mandatory if the property has
            multiple elements. Otherwise, if not given, the only
            element present will be returned.
        method : str, VCDPropertyInterface
                 or VCDPropertyInterface.value, optional
            Which method should be used to access the property/element.
            This is honored only in case the property/element has
            multiple access methods, in which case it is mandatory.

        Returns
        -------
        property_min, property_max : int or float
            The minimum and maximum values of the property/element.
            The return type depends on the method chosen: int for
            VCDPropertyInterface.RANGE, float for .ABSOLUTEVALUE.

        Raises
        ------
        ImagingSourceError
            If errors occur in retrieving property ranges.
        NotImplementedError
            If the property required is only accessible via the
            MAPSTRINGS interface method. Currently unsupported.
        ValueError
            If method is not one of the interface methods known,
            or if this method is called on a SWITCH property.
        """
        (prop_name,
         elem_name,
         method) = self.__get_property_element_interface(property_name,
                                                         element_name,
                                                         method)

        if method is VCDPropertyInterface.ABSOLUTEVALUE:
            range_getter = self._dll_get_property_value_range_float
            prop_min, prop_max = c_float(), c_float()
        elif method is VCDPropertyInterface.RANGE:
            range_getter = self._dll_get_property_value_range_int
            prop_min, prop_max = c_long(), c_long()
        elif method is VCDPropertyInterface.MAPSTRINGS:
            raise NotImplementedError("Map strings not yet implemented")
        else:
            raise ValueError(
                f"Unexpected method {method} for retrieving the range "
                f"of VCD property/element {property_name}/{element_name}"
                )
        range_getter(self.__handle, prop_name, elem_name, prop_min, prop_max)
        return prop_min.value, prop_max.value

    _dll_reset = _dll.IC_ResetProperties
    _dll_reset.argtypes = (GrabberHandlePtr,)
    _dll_reset.errcheck = check_dll_return()

    # TODO: funnily it always returns 0 (== ERROR). Should be
    # checked whether all the properties are actually reset.
    def reset(self):
        """Reset the camera to default properties."""
        self._dll_reset(self.__handle)

    # TODO: does is_trigger_available() [REMOVED] only check for hardware
    # trigger, or also software?

    def send_software_trigger(self):
        """Trigger the camera now via software."""
        self.trigger_vcd_property_once("Trigger", "Software Trigger")

    _dll_set_video_format = _dll.IC_SetVideoFormat
    _dll_set_video_format.argtypes = GrabberHandlePtr, c_char_p
    _dll_set_video_format.errcheck = check_dll_return()

    def set_video_format(self, video_format):
        """Set a video format for the current video-capture device.

        The video format must be supported by the current device.

        Parameters
        ----------
        video_format : str or bytes or bytearray
            The desired video format. Should be in the form
            '<format> (<width>x<height>)', where <format> is
            one of the names in SinkFormat, and <width> and
            <height> are the number of pixels requested along
            the two directions.

        Returns
        -------
        ret_val : {DLLReturns.SUCCESS.value, DLLReturns.ERROR.value}
            ERROR if something went wrong.
        """
        # NOTE:
        # There seems to be a C++ class method get_active_video_format
        # for CaptureDevice, which actually ends up calling the
        # implementation from the specific camera type (e.g., in
        # <arv.h>). However, there seems to be no way to call it via
        # the DLL, i.e., I cannot implement a video_format getter.

        self._dll_set_video_format(self.__handle,
                                   _to_bytes(video_format))
        # The .image_info is not updated until the sink format
        # is also fully updated. We can safely set video and sink
        # color formats equal here.
        self.sink_format = SinkFormat.get(video_format.split()[0])

        # Then update the info in the camera itself by getting back
        # the sink format just set.
        _ = self.sink_format

    def start_delivering_frames(self, readout_rate=1024):
        """Start delivering frames in continuous mode.

        Parameters
        ----------
        readout_rate : int, optional
            .frame_rate to be used for delivering frames.
            Default is setting the maximum frame_rate. Notice
            that, if the exposure time is larger than 1/rate, frames
            will be delivered at time intervals slightly larger than
            exposure. The extra time scales with the size of frames
            in pixels. Moreover, the first frame typically takes longer
            to be delivered (depending on when this method is called
            with respect to the internal clock of the camera, as one
            needs to wait for the first internal 'clock' pulse).
        """
        if self.is_running:
            self.stop()
        self.frame_rate = readout_rate
        self.start('continuous')

    def stop_delivering_frames(self):
        """Stop delivering frames.

        The camera is still open and active, and will accept trigger
        pulses if needed.
        """
        self.trigger_enabled = True

    _dll_trigger_property_once = _dll.IC_PropertyOnePush
    _dll_trigger_property_once.argtypes = GrabberHandlePtr, c_char_p, c_char_p
    _dll_trigger_property_once.errcheck = __vcd_prop_checker

    def trigger_vcd_property_once(self, property_name, element_name):
        """Trigger a property/element automatically once now.

        Parameters
        ----------
        property_name : str or bytes
            The name of the property to be triggered (e.g., 'Trigger')
        element_name : str or bytes or None, optional
            The name of the setting to be accessed (e.g., 'Software
            Trigger').
        """
        self._dll_trigger_property_once(
            self.__handle, _to_bytes(property_name), _to_bytes(element_name)
            )

    _dll_list_video_formats = _dll.IC_ListVideoFormats
    _dll_list_video_formats.argtypes = GrabberHandlePtr, c_char_p, c_int
    _dll_list_video_formats.errcheck = check_dll_return(">=0")

    @property
    def __default_video_formats(self):
        """Return a list of the default video formats available.

        This property is especially useful to determine the minimum
        and maximum width and height that a format supports.

        Returns
        -------
        video_formats : list
            Each element is a string of the form
            '<formatID> (<width>x<height>)'
        """
        # There is probably a better way to do this than making a 2D
        # char array and then casting to char_p (that's what the function
        # expects according to the .h signature).
        n_max_formats = 80
        n_max_chars = 40
        formats_arr = ((c_char * n_max_chars) * n_max_formats)()
        self._dll_list_video_formats(
            self.__handle, c_cast(formats_arr, c_char_p), n_max_chars
            )
        formats = [c_cast(f, c_char_p).value.decode()
                   for f in formats_arr
                   if f[0] != b'\x00']
        return formats

    _dll_n_devices = _dll.IC_GetDeviceCount
    _dll_n_devices.argtypes = None
    _dll_n_devices.errcheck = check_dll_return('>=0')

    @property
    def __n_devices(self):
        """Get the number of the currently available devices.

        This function creates an internal array of all connected video-
        capture devices. With each call to this function, this array
        is rebuilt. The unique name can be retrieved from the internal
        array using the method __unique_name_from_index(). This can
        then be used to open devices.

        Returns
        -------
        n_devices : int
            The number of devices found. DLLReturns.NO_HANDLE.value is
            returned in case of error.
        """
        return self._dll_n_devices()

    # ###############   DISCOVERY OF VCD PROPERTIES   #################
    # Three functions, one for getting the name of properties, one
    # for getting the names of elements of an (already discovered)
    # property, one for getting available interfaces for get/setting
    # the property/element pair.

    _dll_discover_vcd_properties = _dll.IC_enumProperties
    _dll_discover_vcd_properties.errcheck = check_dll_return()
    _dll_discover_vcd_properties.argtypes = (
        GrabberHandlePtr, PropertyDiscoveryCallbackType, py_object
        )

    _dll_discover_vcd_elements = _dll.IC_enumPropertyElements
    _dll_discover_vcd_elements.errcheck = check_dll_return()
    _dll_discover_vcd_elements.argtypes = (
        GrabberHandlePtr, c_char_p, PropertyDiscoveryCallbackType, py_object
        )

    _dll_discover_vcd_interfaces = _dll.IC_enumPropertyElementInterfaces
    _dll_discover_vcd_interfaces.errcheck = check_dll_return()
    _dll_discover_vcd_interfaces.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, PropertyDiscoveryCallbackType,
        py_object
        )

    def __find_vcd_properties(self):
        """Return all available VCD property/element/interface.

        The properties are also stored in the private __vcd_properties
        attribute.

        Returns
        -------
        properties : dict
            keys are property names, values are dict. Each
            key of each sub-dict is an element name, values
            are lists of VCDPropertyInterface attributes.
        """
        self.__vcd_properties = {}
        self._dll_discover_vcd_properties(self.__handle, store_property,
                                          self.__vcd_properties)
        for prop_name, elements in self.__vcd_properties.items():
            self._dll_discover_vcd_elements(
                self.__handle, _to_bytes(prop_name),
                store_property, elements
                )
            for elem_name, interfaces in elements.items():
                self._dll_discover_vcd_interfaces(
                    self.__handle, _to_bytes(prop_name), _to_bytes(elem_name),
                    store_property, interfaces
                    )
                # Replace interface names with VCDPropertyInterface
                elements[elem_name] = [VCDPropertyInterface.get(i)
                                       for i in interfaces]
        return self.__vcd_properties

    def __get_property_element_interface(self, prop_name, elem_name, method):
        """Return a tuple of VCD property/element/interface.

        The return values are suitable to be used in setter, getter,
        and range-getter DLL functions for VCD properties.

        Parameters
        ----------
        prop_name : str or bytes
            Name of the property
        elem_name : str or bytes or None
            Element to be accessed. Can be None only if the property
            has only one element.
        method : str or VCDPropertyInterface or VCDPropertyInterface.value
            Honored only if the property/element has multiple access methods.

        Returns
        -------
        prop_name : bytes
            Name of the property
        elem_name : bytes
            Name of the element
        method : VCDPropertyInterface
            Interface.

        Raises
        ------
        ImagingSourceError
            If this method is called but the device is not open.
        ValueError
            If the element/interface cannot be picked univocally.
        """
        if not self.vcd_properties:
            # Device is not open
            raise ImagingSourceError(
                "Cannot access properties when camera is not open",
                err_code=DLLReturns.NO_DEVICE
                )

        # Get the property
        vcd_prop = self.vcd_properties.get(prop_name, None)
        if vcd_prop is None:
            raise ValueError(f"VCD property {prop_name} not available "
                             "or misspelled. Use self.vcd_properties to "
                             "see the available VCD properties")
        # Pick an element
        if len(vcd_prop) > 1 and not elem_name:
            raise ValueError(f"Property {prop_name} has the following "
                             f"elements: {list(vcd_prop.keys())}. "
                             "Pick one.")
        if len(vcd_prop) == 1:
            elem_name = tuple(vcd_prop)[0]
        vcd_elem = vcd_prop.get(elem_name, None)
        if vcd_elem is None:
            raise ValueError(f"VCD element {elem_name} not available for "
                             f"property {prop_name} (or misspelled). Use "
                             "self.vcd_properties to see available properties")

        # Pick an interface (== method)
        if len(vcd_elem) > 1 and not method:
            raise ValueError(
                f"Property/Element {prop_name}/{elem_name} has the "
                f"following methods: {[i.name for i in vcd_elem]}. Pick one."
                )
        if len(vcd_elem) == 1:
            method = vcd_elem[0]

        method = VCDPropertyInterface.get(method)

        return _to_bytes(prop_name), _to_bytes(elem_name), method

    _dll_set_frame_acquisition_mode = _dll.IC_SetContinuousMode
    _dll_set_frame_acquisition_mode.argtypes = GrabberHandlePtr, c_int
    _dll_set_frame_acquisition_mode.errcheck = check_dll_return()

    def __set_frame_acquisition_mode(self, mode):
        """Set the mode for returning images.

        For correct operation, 'triggered' mode requires mode to be
        StreamMode.CONTINUOUS, .trigger_enabled == True, then start()
        after setting the mode and enabling the trigger [alternatively,
        start('triggered')]. For 'live' mode, both CONTINUOUS and
        SNAP work fine, but .trigger_enabled must be False [this is
        automatic if start('continuous') is called].

        Parameters
        ----------
        mode : StreamMode or str
            If CONTINUOUS, each frame acquired by the camera
            will be passed to the frame-ready callback function
            automatically. If SNAP, images are returned only
            when IC_SnapImage() is called.

        Raises
        ------
        ValueError
            If mode is not a valid StreamMode
        """
        try:
            mode = StreamMode(mode)
        except ValueError as err:
            raise ValueError(f"Invalid stream mode {mode}") from err

        was_running = self.is_running
        if was_running:
            self.stop()

        self._dll_set_frame_acquisition_mode(self.__handle, mode.value)

        # Snap once (if needed) to update the internal
        # settings of the camera
        if mode is StreamMode.SNAP:
            # The next line also sets self.trigger_enabled to False
            self.start('continuous')
            self._dll.IC_SnapImage(self.__handle, 2000)
            self.stop()

        # Restart in the previous mode
        if was_running:
            self.start()

    _dll_unique_name_from_index = _dll.IC_GetUniqueNamefromList
    _dll_unique_name_from_index.restype = c_char_p
    _dll_unique_name_from_index.argtypes = (c_int,)
    _dll_unique_name_from_index.errcheck = check_dll_return('pointer')

    def __unique_name_from_index(self, index):
        """Get unique name of a device given its index.

        The unique device name consist of the device name, any given
        name (in square brackets) and its serial number. It allows to
        differentiate between more devices of the same type connected
        to the computer. The unique device name is passed to the method
        open().

        Parameters
        ----------
        index : int
            The index of the device whose name is to be returned. It
            must be between 0 and self.__n_devices.

        Returns
        -------
        device_name : str
            The unique name of the device. An empty string in case of
            failure.
        """
        try:
            name = self._dll_unique_name_from_index(index)
        except ImagingSourceError as err:
            # Probably some fuckup with the index.
            print(f"IS driver warning: {err}. index={index} self.__n_devices={self.__n_devices}")    # TODO: remove this
            name = None
        if name is not None:
            return name.decode()
        return ''
