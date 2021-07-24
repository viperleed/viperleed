"""Module grabber of camera.drivers.imagingsource.

========================================
   ViPErLEED Graphical User Interface
========================================

This is an extension of the original tisgrabber.py module available
in the original form at:
github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/
          Python/Open%20Camera%2C%20Grab%20Image%20to%20OpenCV

Edits include code (and code-style) improvements.

Author: Michele Riva
Date edits started: 2021-07-18

The original tisgrabber.py module was:
Created on Mon Nov 21 09:44:40 2016
@author: Daniel Vassmer, Stefan_Geissler
"""
from enum import Enum
from ctypes import (Structure as CStructure, WinDLL, CFUNCTYPE, cast as c_cast,
                    sizeof as c_sizeof, c_uint8, c_int, c_long, c_ulong,
                    c_ubyte, c_float, py_object, POINTER, c_char_p, c_void_p)
import inspect
import os
from pathlib import Path
import sys
import numpy as np

from viperleed.guilib.measure.camera.drivers.imagingsource.winerrors import (
    DLLReturns, check_dll_return, ImagingSourceError
    )

dll_path = Path(__file__).parent.resolve()
print(f"{dll_path=}")
try:
    os.add_dll_directory(dll_path)
except AttributeError:
    pass
sys.path.append(dll_path)


c_float_p = POINTER(c_float)

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


class CameraProperty(Enum):
    """Class holding possible camera properties.

    Does not include VCD properties, which are accessed by name.
    """

    PAN = 0
    TILT = 1
    ROLL = 2
    ZOOM = 3
    EXPOSURE = 4
    IRIS = 5
    FOCUS = 6
    MEGA = 65536  # Guard for max enum value, i.e., 2**sizeof(int)

    @classmethod
    def get(cls, cam_property):
        """Return an appropriate enum from value or name.

        Parameters
        ----------
        cam_property : str, int, or CameraProperty
            Some reference to a camera property. If a string, lookup
            is done by name (not case sensitive), otherwise by value.

        Returns
        -------
        property : CameraProperty
            The attribute corresponding to cam_property

        Raises
        ------
        ValueError
            If no property could be found
        """
        if isinstance(cam_property, str):
            try:
                cam_property = getattr(cls, cam_property.upper())
            except AttributeError as err:
                raise ValueError(f"Unknown camera property {cam_property}. "
                                 "Perhaps you meant video property?") from err
        else:
            try:
                cam_property = cls(cam_property)
            except ValueError as err:
                raise ValueError(f"Unknown camera property {cam_property}. "
                                 "Perhaps you meant video property?") from err
        return cam_property


class SinkFormat(Enum):
    Y800 = 0      # monochrome 1-byte;    top-left
    RGB24 = 1     # 3 bytes: B G R;       top-left
    RGB32 = 2     # 4 bytes: B G R alpha; top-left; alpha unused == 0
    UYVY = 3      # 4-bytes color;        top-left; see fourcc.org; bit complicated (U&V shared by neighbor pixels, y for both)
    Y16 = 4       # gray 2-bytes;         top-left
    NONE = 5      # used only as return value
    MEGA = 65536  # Guard for max enum value, i.e., 2**sizeof(int)
    # RGB64       # 8 bytes: 2B 2G 2R 2A; top-left; seems unsupported

    @classmethod
    def get(cls, sink_format):
        """Return an appropriate enum from value or name.

        Parameters
        ----------
        sink_format : str, int or SinkFormat
            Some reference to a sink format.

        Raises
        ------
        ValueError
            If sink_format is not a valid sink format
        ImagingSourceError
            If sink_format would yield .NONE
        RuntimeError
            If sink_format would yield .MEGA
        """
        if isinstance(sink_format, str):
            try:
                sink_format = getattr(cls, sink_format)
            except AttributeError as err:
                raise ValueError("Unknown sink/video format "
                                 f"{sink_format}") from err
        else:
            try:
                sink_format = cls(sink_format)
            except ValueError as err:
                raise ValueError("Unknown sink/video format "
                                 f"{sink_format}") from err
        if sink_format == SinkFormat.NONE:
            raise ImagingSourceError("Invalid sink/video "
                                     f"format {sink_format.name}")
        elif sink_format == SinkFormat.MEGA:
            raise ImagingSourceError("Camera should never return "
                                     "SinkFormat.MEGA")
        return sink_format


class StreamMode(Enum):
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


class SwitchProperty(Enum):
    """On/Off VCD property enum."""

    ON = 1
    OFF = 0

    @classmethod
    def get(cls, on_off):
        """Return an appropriate enum from value or name.

        Parameters
        ----------
        on_off : str, SwitchProperty or object
            If a str, should be 'on'/'off' (not case sensitive).
            If an object, its truth value will be used.

        Returns
        -------
        switch_property : SwitchProperty
            The attribute corresponding to on_off.

        Raises
        ------
        ValueError
            If a string is passed that is neither 'on' nor 'off'
            (not case sensitive)
        """
        if isinstance(on_off, str):
            try:
                switch_property = getattr(cls, on_off.upper())
            except AttributeError as err:
                raise ValueError(
                    f"Invalid string {on_off} for SwitchProperty. "
                    "Should be 'on' or 'off' (not case sensitive)"
                    )
        else:
            switch_property = cls(bool(on_off))
        return switch_property


class GrabberHandle(CStructure):
    _fields_ = [('unused', c_int)]

    def __repr__(self):
        """Return string representation of self."""
        return f"{self.__class__.__name__}({self.unused})"


# ctypes.POINTER returns a class
GrabberHandlePtr = POINTER(GrabberHandle)
GrabberHandlePtr.__repr__ = lambda self: f"GrabberHandlePtr({self.contents})"


# Prototype of frame-ready callbacks:
FrameReadyCallbackType = CFUNCTYPE(
    c_void_p,          # return type
    GrabberHandlePtr,  # HGRABBER; funnily was c_int!
    POINTER(c_ubyte),  # pointer to first byte of image data
    c_ulong,           # current frame number
    py_object          # other data (python) passed along
    )

# A little trick to make functions decorated with
# @FrameReadyCallbackType inspect-able with inspect.signature
# NB: other inspect features may not work, as there is no
# explicit support for CFUNCTYPE() objects.
def __example_cb(__handle, img_buffer, frame_no, py_obj_for_callback):
    """A no-op frame-ready callback with appropriate signature."""
    return
FrameReadyCallbackType.__signature__ = (
    inspect.Signature.from_callable(__example_cb)
    )


class VCDPropertyInterface(Enum):
    """List of VCD property interfaces with their methods."""

    RANGE = int
    SWITCH = bool
    BUTTON = 'one_push'
    MAPSTRINGS = str
    ABSOLUTEVALUE = float
    UNKNOWN = None

    @classmethod
    def get(cls, name):
        """Return attribute by name."""
        if isinstance(name, cls):
            return name
        return getattr(cls, name.upper(), cls.UNKNOWN)



# Prototype of frame-ready callbacks:
PropertyDiscoveryCallbackType = CFUNCTYPE(
    c_int,             # return type
    c_char_p,           # name of property/element/interface
    py_object          # data (python) useful to store the name
    )

@PropertyDiscoveryCallbackType
def store_property(name, storage):
    """Store a property in a python dict.

    Parameters
    ----------
    name : bytes
        Name of property/element/interface
    storage : dict
        Property names will be stored keys, values are dictionaries

    Returns
    -------
    finished : int
        Whether property search is over or not. This function
        always returns 0, such that all properties can be discovered.
    """
    storage[name.decode()] = {}
    return 0 # DLLReturns.SUCCESS.value


# TODO: look at the documentation at https://github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/Python/Open%20Camera%2C%20Grab%20Image%20to%20OpenCV
# TODO: also at (ip-config) https://github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/c%23/Open%20GigE%20Camera%20by%20IP%20Address

# TODO: try removing some of the many functions, especially those we will never use

class  WindowsCamera: #GrabberDLL:
    # """Interface class that defines python equivalents to the DLL calls."""
    """Interface class for Imaging Source camera in Windows."""

    # To prevent creating bindings to the dll functions each time
    # methods are called, each method that calls a dll function
    # has ctypes definition of the function used right above it.
    # Note: ctypes assumes restype == c_int. Thus restype is
    # omitted below where this is the case.

    if sys.maxsize > 2**32:
        _dll = WinDLL("tisgrabber_x64.dll")
    else:
        _dll = WinDLL("tisgrabber.dll")

    __initalized = False

    def __init__(self):
        """Initialize grabber instance."""
        # Call init_library with None straight away. It was done the
        # same way in the original tisgrabber.py. init_library will
        # do something only the first time an instance is created.
        self.init_library()
        self.__handle = self.create_grabber()

    @classmethod
    def init_library(cls, license_key=None):
        """Initialize the IC Imaging Control dll library.

        This function must be called only once before any other
        functions of the dll are called.

        Parameters
        ----------
        license_key : str or None, optional
            The license key. Pass None if using a 'trial' version.
            Default is None.

        Returns
        -------
        ret_val : int
            SUCCESS on success, ERROR if license key is wrong or for
            other generic errors.
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

        Unused grabber handles should be release by calling
        IC_ReleaseGrabber(). There is no need to release a handle
        when closing a device, as handles can be reused.

        Returns
        -------
        handle : GrabberHandlePtr
            New handle. Can be used to open devices.
        """
        return self._dll_create_grabber()

    # TODO: an alternative to do the loop in C would be to use IC_ListDevices
    # but it is unclear whether it would return unique names or just device
    # names. Perhaps worth implementing and testing. Then we would not need
    # n_devices and unique_name_from_index.
    @property
    def devices(self):
        """Return the unique names of available devices.

        Returns
        -------
        unique_names : list
            List of unique device names, each as a string.
        """
        names = (self.unique_name_from_index(i) for i in range(self.n_devices))
        return [n for n in names if n]

    _dll_get_frame_rate =  _dll.IC_GetFrameRate
    _dll_get_frame_rate.restype = c_float
    _dll_get_frame_rate.argtypes = (GrabberHandlePtr,)
    _dll_get_frame_rate.errcheck = check_dll_return('>0')

    @property
    def frame_rate(self):
        """Return the frame rate of the current device.

        The device should be open and in live mode.

        Returns
        -------
        frame_rate : float
            A positive frame rate if frame rate is supported,
            zero otherwise.
        """
        # TODO: check it is OK. The original implementation was first
        # defining frame_rate_c = ctypes.c_long(), then overwriting
        # it with the return value
        return self._dll_get_frame_rate(self.__handle)

    _dll_set_frame_rate = _dll.IC_SetFrameRate
    _dll_set_frame_rate.argtypes = GrabberHandlePtr, c_float
    _dll_set_frame_rate.errcheck = check_dll_return(
        # NOT_AVAILABLE would be an error, but NOT_IN_LIVEMODE
        # seems more reasonable to return. Perhaps better check
        # live mode explicitly and leave NOT_AVAILABLE.
        # TODO: cross check exclude, as it seems not working properly.
        #       got a WRONG_INCOMPATIBLE_XML error which should have been
        #       NOT_IN_LIVEMODE
        exclude_errors=('NOT_AVAILABLE', 'NO_PROPERTYSET',
                        'DEFAULT_WINDOW_SIZE_SET', 'WRONG_XML_FORMAT',
                        'WRONG_INCOMPATIBLE_XML',
                        'NOT_ALL_PROPERTIES_RESTORED', 'DEVICE_NOT_FOUND')
        )

    @frame_rate.setter
    def frame_rate(self, rate):
        """Set a new frame rate.

        A new frame rate can not be set while in live mode. Call
        stop_live(), then try again.

        Parameters
        ----------
        frame_rate : float
            The new frame rate.
        """
        # TODO: if live, stop, then set, then restart
        self._dll_set_frame_rate(self.__handle, rate)

    # TODO: unclear what 'color format' means. Is it == video_format??
    #       is it == sink_type ??
    # TODO: eventually it would be nice to align the returns to the
    #       signature I want in CameraABC (i.e., return the no. of
    #       color channels rather than the format)
    _dll_get_image_description = _dll.IC_GetImageDescription
    _dll_get_image_description.argtypes = (
        GrabberHandlePtr, POINTER(c_long), POINTER(c_long),
        POINTER(c_int), POINTER(c_int)
        )
    _dll_get_image_description.errcheck = check_dll_return()

    @property
    def image_info(self):
        """Return properties of the current video format and sink type.

        Returns
        -------
        width, height : int
            Width and height of the image in pixels
        bytes_per_pixel : int
            Number of bytes for each pixel, accounting for color
            channels
        format : SinkFormat
            The current color format

        Raises
        ---------
        ImagingSourceError
            If the format returned is unknown.
        """
        # TODO: this returns an ERROR if the sink format was changed
        # but no image was ever snapped since then. The way to go:
        # check if it errors out with return == 0. If this is the case,
        # .start_live(), .stop_live(), then call again. The second time
        # it should not fail if the camera is fine.
        width, height = c_long(), c_long()
        bits_per_pixel, color_format_c = c_int(), c_int()

        self._dll_get_image_description(self.__handle, width, height,
                                        bits_per_pixel, color_format_c)
        bytes_per_pixel = bits_per_pixel.value // 8
        try:
            color_format = SinkFormat.get(color_format_c.value)
        except ValueError as err:
            raise ImagingSourceError(err.args[0] + " received.") from err

        return width.value, height.value, bytes_per_pixel, color_format

    _dll_n_devices = _dll.IC_GetDeviceCount
    _dll_n_devices.argtypes = None
    _dll_n_devices.errcheck = check_dll_return('>=0')

    @property
    def n_devices(self):    # TODO: do we need this to be public??
        """Get the number of the currently available devices.

        This function creates an internal array of all connected video-
        capture devices. With each call to this function, this array is
        rebuilt. The name and the unique name can be retrieved from the
        internal array using the function IC_GetDevice() and the method
        unique_name_from_index(). Both are useful to for retrieving
        device names to open devices.

        Returns
        -------
        n_devices : int
            The number of devices found. DLLReturns.NO_HANDLE.value is
            returned in case of error.
        """
        return self._dll_n_devices()

    _dll_get_sink_format = _dll.IC_GetFormat
    _dll_get_sink_format.argtypes = (GrabberHandlePtr,)

    @property
    def sink_format(self):
        """Return the color format of the sink type currently set.

        The sink type describes the format of the buffer where snapped
        images are stored. The method returns a valid value only after
        IC_PreprareLive() or start_live() are called. Notice that the
        sink type may differ from the video format.

        Returns
        -------
        color_format : int
            One of the values in SinkFormat. Returns !!UNCLEAR WHAT!!
            if not in live mode or if an error occurred.
        """
        # TODO: needs live mode? NO, requires the camera to have acquired
        # at least one frame. Before that, it errors out returning NONE.
        # TODO: old implementation was returning RGB24 in case
        #       of 'unknown' format
        sink_format_c = self._dll_get_sink_format(self.__handle)
        try:
            sink_format = SinkFormat.get(sink_format_c)
        except ValueError as err:
            raise ImagingSourceError(err.args[0] + " received")
        except ImagingSourceError as err:
            raise ImagingSourceError("Camera should be in live mode. Call "
                                     ".start_live(), then retry.") from err

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
        sink_type : str, int, or SinkFormat
            One of the values in SinkFormat. Note that UYVY can
            only be used in conjunction with a UYVY video format
        """
        # TODO: whenever changed, also .start_live() .stop_live() to
        # make sure that the getter returns data
        sink_format = SinkFormat.get(sink_format)
        self._dll_set_sink_format(self.__handle, sink_format.value)

    _dll_is_device_valid = _dll.IC_IsDevValid
    _dll_is_device_valid.argtypes = (GrabberHandlePtr,)
    _dll_is_device_valid.errcheck = check_dll_return()

    @property
    def valid(self):  # TODO: maybe I don't need this public
        """Return whether there is a valid open device."""
        # TODO: after checking possible returns, get rid of
        # the error checker if the only thing that can happen
        # is a generic ERROR. The purpose of the error checker
        # is just printing out the return values.
        try:
            ret_val = self._dll_is_device_valid(self.__handle)
        except ImagingSourceError:
            # No device open
            pass
        return  ret_val == DLLReturns.SUCCESS.value

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
        return width.value, height.value

    # TODO: perhaps would be nice to run the loop at the c level. This
    # can be done with IC_ListVideoFormats. Then one would not need
    # __n_video_formats() nor __video_format_by_index(idx)
    @property
    def video_formats(self):
        """Return a list of the video formats available for a device."""
        n_formats = self.__n_video_formats()
        formats = (self.__video_format_by_index(i) for i in range(n_formats))
        return [f for f in formats if f]

    _dll_set_video_format = _dll.IC_SetVideoFormat
    _dll_set_video_format.argtypes = GrabberHandlePtr, c_char_p
    _dll_set_video_format.errcheck = check_dll_return()

    # TODO: funnily enough there is no way to get the video format
    # back as a string. The only way to fetch info about the
    # video format is via get_image_description()
    def set_video_format(self, video_format):
        """Set a video format for the current video-capture device.

        The video format must be supported by the current device.

        Parameters
        ----------
        video_format : str, bytes, or bytearray
            The desired video format.

        Returns
        -------
        ret_val : {DLLReturns.SUCCESS.value, DLLReturns.ERROR.value}
            ERROR if something went wrong.
        """
        self._dll_set_video_format(self.__handle,
                                   _to_bytes(video_format))

    # TODO: check which errors
    _dll_open_by_unique_name = _dll.IC_OpenDevByUniqueName
    _dll_open_by_unique_name.argtypes = GrabberHandlePtr, c_char_p
    _dll_open_by_unique_name.errcheck = check_dll_return()

    def open_by_unique_name(self, unique_name):
        """Open a video capture by using its unique name.

        Use IC_GetUniqueName() to retrieve the unique name of an
        open device. unique_name_from_index(index) can be used
        before the device is open to retrieve the unique name.

        Parameters
        ----------
        unique_name : bytes
            Unique name of the device to be opened
        """
        self._dll_open_by_unique_name(self.__handle,
                                      _to_bytes(unique_name))

    # TODO: may it return other errors, like NO_HANDLE and NO_DEVICE?
    # TODO: include_errors / exclude_errors doe not work properly.
    #       maybe I have it the wrong way around?
    _dll_snap_image = _dll.IC_SnapImage
    _dll_snap_image.argtypes = GrabberHandlePtr, c_int
    _dll_snap_image.errcheck = check_dll_return(
        include_errors=('ERROR', 'NOT_IN_LIVEMODE')
        )

    # Would be nicer for it to actually return a np.array already
    # I could then call __get_image_ptr in here rather than having
    # it as a separate method
    def snap_live_image(self, timeout=2000):
        """Snap an image.

        The video capture device must be set to live mode and a sink
        type has to be set before this call. The format of the snapped
        images depends on the selected sink type.

        Parameters
        ----------
        timeout : int, optional
            Time interval in milliseconds after which the device will
            time out. A value of -1 corresponds to no time-out. Default
            is 2000, i.e., 2 sec.

        Returns
        -------
        None
        """
        # C-Returns
        # ---------
        # ret_val : int
            # SUCCESS if an image has been snapped, NOT_IN_LIVEMODE if
            # the device is not in live mode, i.e., start_live() was
            # not called, ERROR if something else went wrong.
        self._dll_snap_image(self.__handle, timeout)

    _dll_start_live = _dll.IC_StartLive
    _dll_start_live.argtypes = GrabberHandlePtr, c_int
    _dll_start_live.errcheck = check_dll_return()

    def start_live(self, show_video=False):
        """Start camera in live mode.

        Parameters
        ----------
        show_video : bool, optional
            Do not show a video if False, show one if True. Frames
            will be delivered in any case (i.e., callbacks can be
            used). Default is False.

        Returns
        -------
        ret_val : int
            SUCCESS on success, ERROR if something went wrong.
        """
        self._dll_start_live(self.__handle, bool(show_video))

    # TODO: unclear which errors it may return
    _dll_stop_live = _dll.IC_StopLive
    _dll_stop_live.argtypes = (GrabberHandlePtr,)
    _dll_stop_live.errcheck = check_dll_return()

    def stop_live(self):
        """Stop live mode."""
        self._dll_stop_live(self.__handle)

    _dll_show_device_selection_dialog = _dll.IC_ShowDeviceSelectionDialog
    _dll_show_device_selection_dialog.restype = GrabberHandlePtr
    _dll_show_device_selection_dialog.argtypes = (GrabberHandlePtr,)

    # TODO: implemented just for testing. will be later removed.
    # TODO: looks like it should never error out.
    def temp_selection_dialog(self):
        """Show a device-selection dialog.

        This dialog allows to select the video capture device, the
        video norm, video format, input channel and frame rate.

        Returns
        -------
        None.
        """
        self.__handle = self._dll_show_device_selection_dialog(self.__handle)

    # TODO: restype did not match signature. Also, the old python
    # implementation assigned the return to __handle. Should be
    # checked.
    # TODO: implemented just for testing. will be later removed.
    # Dialog looks nice. Perhaps we should mock it in Qt and make
    # it look more or less like this for all cameras.
    _dll_show_property_dialog = _dll.IC_ShowPropertyDialog
    _dll_show_property_dialog.restype = c_int  # was GrabberHandlePtr
    _dll_show_property_dialog.argtypes = (GrabberHandlePtr,)
    _dll_show_property_dialog.errcheck = check_dll_return()

    def temp_property_dialog(self):
        """Show the VCDProperty dialog.

        Requires a video-capture device to be open.

        Returns
        -------
        ret_val : int
            SUCCESS on success, NO_HANDLE if grabber_handle is invalid,
            NO_DEVICE if no device was open, ERROR otherwise.
        """
        self._dll_show_property_dialog(self.__handle)

    _dll_unique_name_from_index = _dll.IC_GetUniqueNamefromList
    _dll_unique_name_from_index.restype = c_char_p
    _dll_unique_name_from_index.argtypes = (c_int,)
    _dll_unique_name_from_index.errcheck = check_dll_return('pointer')

    def unique_name_from_index(self, index):
        """Get unique name of a device given its index.

        The unique device name consist of the device name and its
        serial number. It allows to differentiate between more devices
        of the same type connected to the computer. The unique device
        name is passed to the method open_by_unique_name().

        Parameters
        ----------
        index : int
            The index of the device whose name is to be returned. It
            must be between 0 and self.n_devices.

        Returns
        -------
        device_name : str
            The unique name of the device. An empty string in case of
            failure.
        """
        name = self._dll_unique_name_from_index(index)
        if name is not None:
            return name.decode()
        return ''

    _dll_video_format_by_index = _dll.IC_GetVideoFormat
    _dll_video_format_by_index.restype = c_char_p
    _dll_video_format_by_index.argtypes = GrabberHandlePtr, c_int
    _dll_video_format_by_index.errcheck = check_dll_return('pointer')

    def __video_format_by_index(self, index):
        """Return the video format specified by index as bytes.

        Parameters
        ----------
        index : int
            Index of the video format to be used. Must be between 0 and
            n_video_formats().

        Returns
        -------
        video_format : str
            The name of the specified video format. Returns an empty
            string in case of error.
        """
        video_format = self._dll_video_format_by_index(self.__handle, index)
        if video_format:
            return video_format
        return ''

    _dll_n_video_formats = _dll.IC_GetVideoFormatCount
    _dll_n_video_formats.argtypes = (GrabberHandlePtr,)
    _dll_n_video_formats.errcheck = check_dll_return('>=0')

    def __n_video_formats(self):
        """Return the number of video formats supported.

        Requires a video-capture device to be open.

        Returns
        -------
        n_formats : int
            The number of video formats supported by the open device
            if successful. Otherwise NO_DEVICE (no device open) or
            NO_HANDLE (invalid grabber_handle).
        """
        return self._dll_n_video_formats(self.__handle)

    _dll_get_image_ptr = _dll.IC_GetImagePtr
    _dll_get_image_ptr.restype = c_void_p  # void to also allow None
    _dll_get_image_ptr.argtypes = (GrabberHandlePtr,)
    _dll_get_image_ptr.errcheck = check_dll_return('pointer')

    def __get_image_ptr(self):
        """Get a pointer to the pixel data of the last snapped image.

        Returns
        -------
        pointer : ctypes.c_void_p
            A pointer to the the first byte of the bottommost line
            of the image (images are saved from bottom to top!).
            Returns None in case of error.
        """
        return self._dll_get_image_ptr(self.__handle)


    _dll_is_cam_property_available = _dll.IC_IsCameraPropertyAvailable
    _dll_is_cam_property_available.argtypes = GrabberHandlePtr, c_int
    _dll_is_cam_property_available.errcheck = (
        check_dll_return(exclude_errors=('ERROR',))  # == not available
        )

    def check_available_cam_properties(self):
        """Return a dict of available camera properties.

        The dictionary is also stored for check purposes.

        Returns
        -------
        properties : dict
            Availability of camera properties. Each entry is
            property_name: str
            available: bool
        """
        self.__available_cam_properties = {}
        for cam_property in CameraProperty:
            try:
                available = self._dll_is_cam_property_available(
                    self.__handle, cam_property.value
                    )
            except ImagingSourceError:
                # Check I get ERROR, i.e., not available
                continue
            self.__available_cam_properties[cam_property.name] = (
                available == DLLReturns.SUCCESS.value
                )
        return self.__available_cam_properties

    # TODO: unclear what it returns exactly.
    # TODO: check also whether setting exposure this way is the
    #       same as using set_property_value_int and absolute value
    #       If it is the same, probably better keep only one to
    #       limit overcrowding of class.
    _dll_get_camera_property = _dll.IC_GetCameraProperty
    _dll_get_camera_property.argtypes = (GrabberHandlePtr, c_int,
                                          POINTER(c_long))
    _dll_get_camera_property.errcheck = check_dll_return()

    def get_camera_property(self, cam_property):
        """Retrieve the value of a camera property.

        Parameters
        ----------
        cam_property : str, int, or CameraProperty
            The property whose value is to be read

        Returns
        -------
        property_value : int
            Value of the requested property

        Raises
        ------
        ValueError
            If cam_property is not a valid camera property
        ImagingSourceError
            if cam_property is not available
        """
        # C-Returns
        # ---------
        # ret_val : int
        #     SUCCESS or ERROR. Probably also NOT_AVAILABLE, and perhaps
        #     PROPERTY_ITEM_NOT_AVAILABLE, PROPERTY_ELEMENT_NOT_AVAILABLE,
        #     PROPERTY_ELEMENT_WRONG_INTERFACE
        cam_property = CameraProperty.get(cam_property)
        if not self.__available_cam_properties:
            self.check_available_cam_properties()
        if not self.__available_cam_properties.get(cam_property.name):
            raise ImagingSourceError(
                f"Property {cam_property.name} not available for this device"
                )
        prop = c_long()
        self._dll_get_camera_property(self.__handle, cam_property.value, prop)
        return prop.value

    _dll_set_camera_property = _dll.IC_SetCameraProperty
    _dll_set_camera_property.argtypes = GrabberHandlePtr, c_int, c_long
    _dll_set_camera_property.errcheck = check_dll_return()

    def set_camera_property(self, cam_property, value):
        """Set a camera property.

        Parameters
        ----------
        property : str, int, or CameraProperty
            The property to be set
        value : int
            The value of the property. The value should be in the range
            of the specified property.

        Raises
        ------
        ValueError
            If cam_property is not a valid camera property
        ImagingSourceError
            If cam_property is not available or camera returned error
        """

        cam_property = CameraProperty.get(cam_property)
        if not self.__available_cam_properties:
            self.check_available_cam_properties()
        if not self.__available_cam_properties.get(cam_property.name):
            raise ImagingSourceError(
                f"Property {cam_property.name} not available"
                )
        self._dll_set_camera_property(self.__handle, cam_property.value,
                                      value)

    # TODO: check return values
    _dll_get_camera_property_range = _dll.IC_CameraPropertyGetRange
    _dll_get_camera_property_range.argtypes = (GrabberHandlePtr, c_int,
                                               POINTER(c_long),
                                               POINTER(c_long))
    _dll_get_camera_property_range.errcheck = check_dll_return()

    def get_camera_property_range(self, cam_property):
        """Retrieve the minimum and maximum value of a property.

        Parameters
        ----------
        property : str, int, or CameraProperty
            The property whose range is to be read

        Returns
        -------
        min, max : int
            The minimum and maximum value of the property.

        Raises
        ------
        ValueError
            If cam_property is not a valid camera property
        RuntimeError
            If cam_property is not available
        """
        cam_property = CameraProperty.get(cam_property)
        if not self.__available_cam_properties:
            self.check_available_cam_properties()
        if not self.__available_cam_properties.get(cam_property.name):
            raise ImagingSourceError(
                f"Property {cam_property.name} not available"
                )
        prop_min, prop_max = c_long(), c_long()

        self._dll_get_camera_property_range(self.__handle, cam_property.value,
                                            prop_min, prop_max)
        return prop_min, prop_max

    _dll_is_property_available =  _dll.IC_IsPropertyAvailable
    _dll_is_property_available.argtypes = GrabberHandlePtr, c_char_p, c_char_p
    _dll_is_property_available.errcheck = check_dll_return()
    # TODO: exclude PROPERTY_ITEM_NOT_AVAILABLE, PROPERTY_ELEMENT_NOT_AVAILABLE,
    # NOT_ALL_PROPERTIES_RESTORED (== PROPERTY_ITEM_NOT_AVAILABLE),
    # DEVICE_NOT_FOUND (== PROPERTY_ELEMENT_NOT_AVAILABLE).
    # Perhaps also add a valid_returns tuple to the error checker!

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

    def find_vcd_properties(self):
        """Return available VCD property/element/interface.

        Returns
        -------
        properties : dict
            keys are property names, values are dict. Each sub-dict
            key is an element name, values are VCDPropertyInterface
            attributes.
        """
        vcd_properties = {}
        self._dll_discover_vcd_properties(self.__handle, store_property,
                                          vcd_properties)
        for prop_name, elements in vcd_properties.items():
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


        print("\n#######################################################\n",
              vcd_properties)



    def is_vcd_property_available(self, property_name, element_name=None):
        """Check if a (property, element) pair exists.

        Parameters
        ----------
        property_name : str or bytes
            Name of the property to be checked
        element_name : str or bytes, optional
            Property element, e.g., "Value", "Auto". If not given
            or None, only the existence of the property is checked.
            Default is None.

        Returns
        -------
        available : bool
            True if the (property, element) pair exists.
        """
        # C-Returns
        # ---------
        # SUCCESS Success,  NO_HANDLE, NO_DEVICE,
        # PROPERTY_ELEMENT_NOT_AVAILABLE,
        # PROPERTY_ELEMENT_WRONG_INTERFACE
        available = self._dll_is_property_available(
            self.__handle, _to_bytes(property_name), _to_bytes(element_name)
            )
        return available == DLLReturns.SUCCESS.value

    # #####################    VCD PROPERTIES    ######################
    #
    # The way properties are 'set', 'gotten' or 'range-gotten' is
    # always specified via the <method> keyword argument. Values
    # given are adjusted accordingly.

    # TODO: it is undocumented which properties are available.
    # Perhaps using VCDPropertyInspector of IC Imaging Control
    # would be a way to check whether enumerate_properties works?
    # A decent read: https://www.theimagingsource.com/support/documentation
    #                        /ic-imaging-control-cpp/tech_VCDProperties.htm

    # Value as float: getter, range getter, setter.
    _dll_get_property_value_float = _dll.IC_GetPropertyAbsoluteValue
    _dll_get_property_value_float.errcheck = check_dll_return()
    _dll_get_property_value_float.argtypes = (
        GrabberHandlePtr, c_char_p,c_char_p, c_float_p
        )
    _dll_get_property_value_range_float = _dll.IC_GetPropertyAbsoluteValueRange
    _dll_get_property_value_range_float.errcheck = check_dll_return()
    _dll_get_property_value_range_float.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_float_p, c_float_p
        )
    _dll_set_property_value_float = _dll.IC_SetPropertyAbsoluteValue
    _dll_set_property_value_float.errcheck = check_dll_return()
    _dll_set_property_value_float.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_float
        )

    # Value as int: getter, range getter, setter.
    _dll_get_property_value_int = _dll.IC_GetPropertyValue
    _dll_get_property_value_int.errcheck = check_dll_return()
    _dll_get_property_value_int.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, POINTER(c_long)
        )
    _dll_get_property_value_range_int = _dll.IC_GetPropertyValueRange
    _dll_get_property_value_range_int.errcheck = check_dll_return()
    _dll_get_property_value_range_int.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, POINTER(c_long), POINTER(c_long)
        )
    _dll_set_property_value_int = _dll.IC_SetPropertyValue
    _dll_set_property_value_int.errcheck = check_dll_return()
    _dll_set_property_value_int.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_int
        )

    # On/Off (== bool): getter, setter. [NO RANGE OBVIOUSLY]
    _dll_get_on_off_property = _dll.IC_GetPropertySwitch
    _dll_get_on_off_property.errcheck = check_dll_return()
    _dll_get_on_off_property.argtypes= (
        GrabberHandlePtr, c_char_p, c_char_p, POINTER(c_long)
        )
    _dll_set_on_off_property = _dll.IC_SetPropertySwitch
    _dll_set_on_off_property.errcheck = check_dll_return()
    _dll_set_on_off_property.argtypes= (
        GrabberHandlePtr, c_char_p, c_char_p, c_int
        )

    # TODO: exclude DEVICE_NOT_FOUND, NOT_ALL_PROPERTIES_RESTORED,
    # NO_PROPERTYSET, DEFAULT_WINDOW_SIZE_SET
    def get_vcd_property(self, property_name, element_name=None, method=float):
        """Return the current value of a property/element pair.

        Parameters
        ----------
        property_name : str or bytes
            The name of the property whose value is to be changed
        element_name : str, bytes or None, optional
            The name of the setting to be accessed (examples: 'Value',
            'Auto'). If None or not given, 'Value' is used. Default is
            None.
        method : {float, int, bool}, optional
            Which method to use to set the property. Typically, float
            gives finer control over the values, but may not be a valid
            way of setting some properties. bool is reserved for on/off
            properties. Default is float.

        Returns
        -------
        property_value : int, float or SwitchProperty
            The value of the property/element.
        """
        if method is float:
            getter = self._dll_get_property_value_float
            property_value = c_float()
        elif method is int:
            getter = self._dll_get_property_value_int
            property_value = c_long()
        elif method is bool:
            getter = self._dll_get_on_off_property
            property_value = c_long()
        else:
            raise ValueError(
                f"Unexpected method {method} for getting VCD "
                f"property/element {property_name}/{element_name}"
                )
        getter(self.__handle, _to_bytes(property_name),
               _to_bytes(element_name), property_value)

        if method is not bool:
            return property_value.value
        return SwitchProperty.get(property_value.value)

    def set_vcd_property(self, property_name, element_name,
                         value, method=float):
        """Set the value of a VCD property/element pair.

        Parameters
        ----------
        property_name : str or bytes
            The name of the property whose value is to be changed
        element_name : str or bytes
            The name of the setting to be accessed (examples: 'Value',
            'Auto').
        value : object
            The value to be set for the property/element. For on/off
            settings (method == bool) accepts strings 'on'/'off' (not
            case sensitive), SwitchProperty, or any object with a truth
            value. Otherwise it should be a number or a string that
            evaluates to a number.
        method : {float, int, bool}, optional
            Which method to use to set the property. Typically, float
            gives finer control over the values, but may not be a valid
            way of setting some properties. bool is reserved for on/off
            properties. Default is float.

        Raises
        ------
        ImagingSourceError
            If errors occur in setting properties.
        """
        # By exclusion, it looks like the int method is using the
        # DirectShow interface, i.e., only setting integer values
        # to properties. Differently from the setter (take 16-bit
        # int), this writes to a 32-bits long pointer.
        #
        #
        # doc: "The interface IVCDRangeProperty, for example, provides
        #      access to the property that is equivalent to the property
        #      access of DirectShow."
        #
        # tisgrabber.h suggests: "For a list of properties and elements
        # use the VCDPropertyInspector of IC Imaging Control."
        #
        # C-Returns
        # ---------
        # ret_val : int
        #     SUCCESS, NO_HANDLE (invalid grabber handle), NO_DEVICE
        #     (no device open), PROPERTY_ITEM_NOT_AVAILABLE (the
        #     property item is not available),
        #     PROPERTY_ELEMENT_NOT_AVAILABLE (the element of the
        #     property is not available),
        #     PROPERTY_ELEMENT_WRONG_INTERFACE (the property element
        #     exists but it does not have the requested attribute,
        #     e.g., b'Exposure' b'Auto' has no 'range').

        if method is float:
            setter = self._dll_set_property_value_float
            value = float(value)
        elif method is int:
            setter = self._dll_set_property_value_int
            value = int(value)
        elif method is bool:
            setter = self._dll_set_on_off_property
            value = SwitchProperty.get(value).value
        else:
            raise ValueError(
                f"Unexpected method {method} for setting VCD "
                f"property/element {property_name}/{element_name}"
                )
        setter(self.__handle, _to_bytes(property_name),
               _to_bytes(element_name), value)

    def get_vcd_property_range(self, property_name, element_name=None,
                               method=float):
        """Return minimum and maximum values for a property/element pair.

        Parameters
        ----------
        property_name : str or bytes
            The name of the property whose value is to be changed
        element_name : str, bytes or None, optional
            The name of the setting to be accessed (examples: 'Value',
            'Auto'). If None or not given, 'Value' is used. Default is
            None.

        Returns
        -------
        property_min, property_max : int
            The minimum and maximum values of the property/element.
        """
        # C-Returns
        # ---------
        # SUCCESS, NO_HANDLE, NO_DEVICE, PROPERTY_ITEM_NOT_AVAILABLE,
        # PROPERTY_ELEMENT_NOT_AVAILABLE, PROPERTY_ELEMENT_WRONG_INTERFACE
        if method is float:
            range_getter = self._dll_get_property_value_range_float
            prop_min, prop_max = c_float(), c_float()
        elif method is int:
            range_getter = self._dll_get_property_value_range_int
            prop_min, prop_max = c_long(), c_long()
        else:
            raise ValueError(
                f"Unexpected method {method} for retrieving the range "
                f"of VCD property/element {property_name}/{element_name}"
                )
        range_getter(self.__handle, _to_bytes(property_name),
                     _to_bytes(element_name), prop_min, prop_max)
        return prop_min.value, prop_max.value

    _dll_trigger_property_once = _dll.IC_PropertyOnePush
    _dll_trigger_property_once.argtypes = GrabberHandlePtr, c_char_p, c_char_p
    _dll_trigger_property_once.errcheck = check_dll_return()

    def trigger_vcd_property_once(self, property_name, element_name):
        """Trigger a property/element automatically once now.

        Parameters
        ----------
        property_name : str or bytes
            The name of the property to be triggered (e.g., 'Trigger')
        element_name : str, bytes, or None, optional
            The name of the setting to be accessed (e.g., 'Software
            Trigger').
        """
        self._dll_trigger_property_once(
            self.__handle, _to_bytes(property_name), _to_bytes(element_name)
            )

    # The py_object at the end is the same that will later be passed
    # on to the frame-ready callback as the last argument.
    _dll_set_frame_ready_callback = _dll.IC_SetFrameReadyCallback
    _dll_set_frame_ready_callback.errcheck = check_dll_return()
    _dll_set_frame_ready_callback.argtypes = (
        GrabberHandlePtr, FrameReadyCallbackType, py_object
        )

    # TODO: it is not entirely clear how the callback should actually look.
    # printing stuff caused snap to fail; passing a list in which stuff
    # is written failed. Perhaps I need a deeper structure than a list,
    # like a dict? Then the CallbackUserData class from Bernhard may make
    # some sense to exist.
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
            on_frame_ready.

        Raises
        -------
        TypeError
            If on_frame_ready is not callable
        ValueError
            If it looks like on_frame_ready was not decorated
            with @FrameReadyCallbackType.
        """
        # Make sure the signature fits:
        try:
            callback_signature = inspect.signature(on_frame_ready)
        except TypeError as err:
            raise TypeError(f"Frame-ready callback {on_frame_ready.__name__} "
                            "is not a callable") from err
        if (on_frame_ready.__class__.__name__ != 'CFunctionType'
            or len(callback_signature.parameters) != 4):
            raise ValueError(
                f"Frame-ready callback was not decorated with "
                f"@{self.__class__.__name__}.FrameReadyCallbackType. This "
                " is needed to ensure appropriate type checking/conversions"
                )
        self._dll_set_frame_ready_callback(self.__handle, on_frame_ready,
                                           py_obj_for_callback)

    _dll_set_continuous_mode = _dll.IC_SetContinuousMode
    _dll_set_continuous_mode.argtypes = GrabberHandlePtr, c_int
    _dll_set_continuous_mode.errcheck = check_dll_return()

    # TODO: it is unclear whether this will return frames while in live
    # mode, and especially WHERE
    def set_frame_acquisition_mode(self, mode):
        """Set the mode for returning images.

        This method will fail if the device is currently streaming
        in live-view mode.

        Parameters
        ----------
        mode : StreamMode
            If CONTINUOUS, each frame acquired by the camera
            will be passed to the frame-ready callback function
            automatically. If SNAP, images are returned only
            when snap_live_image() is called.
        """
        try:
            mode = StreamMode(mode)
        except ValueError as err:
            raise ValueError(f"Invalid stream mode {mode}")
        self._dll_set_continuous_mode(self.__handle, mode.value)

    # TODO: signature was wrong: c_void_p instead of GrabberHandlePtr,
    # and also had a c_char_p that does not fit with the .h
    _dll_close_device = _dll.IC_CloseVideoCaptureDevice
    _dll_close_device.restype = None  # Returns void
    _dll_close_device.argtypes = (GrabberHandlePtr,)
    # _dll_close_device.errcheck = check_dll_return() # TODO: Returns some random positive number??

    def close_device(self):
        """Close the currently open video-capture device."""
        self._dll_close_device(self.__handle)

    _dll_is_trigger_available = _dll.IC_IsTriggerAvailable
    _dll_is_trigger_available.argtypes = (GrabberHandlePtr,)
    _dll_is_trigger_available.errcheck = (
        check_dll_return(exclude_errors=('ERROR',))  # == not available
        )

    # TODO: is this only hardware trigger, or also software?
    def is_trigger_available(self):
        """Return whether the device supports triggering."""
        available = self._dll_is_trigger_available(self.__handle)
        return available == DLLReturns.SUCCESS.value

    _dll_enable_trigger = _dll.IC_EnableTrigger
    _dll_enable_trigger.argtypes = GrabberHandlePtr, c_int
    _dll_enable_trigger.errcheck = check_dll_return()

    def enable_trigger(self, enabled=True):
        """Enable or disable capability to respond to a trigger signal.

        Trigger signals may either be via hardware or via software.

        Parameters
        ----------
        enabled : bool, optional
            If True, the device will afterwards be responsive to a
            (hardware or software) trigger signal. Default is True.

        Returns
        -------
        None
        """
        _enabled = SwitchProperty(bool(enabled)).value
        self._dll_enable_trigger(self.__handle, _enabled)

    _dll_send_software_trigger = _dll.IC_SoftwareTrigger
    _dll_send_software_trigger.argtypes = (GrabberHandlePtr,)
    _dll_send_software_trigger.errcheck = check_dll_return()

    def send_software_trigger(self):
        """Trigger the camera now via software."""
        self._dll_send_software_trigger(self.__handle)

    _dll_get_trigger_modes = _dll.IC_GetTriggerModes
    _dll_get_trigger_modes.argtypes = GrabberHandlePtr, c_char_p, c_int
    _dll_get_trigger_modes.errcheck = check_dll_return('>0')

    @property
    def available_trigger_modes(self):
        """Return a list of trigger modes available (each a string)."""
        # Perhaps this only refers to hardware trigger modes?
        # I get 0 returned on the DMK 33GX265 [specs] 46910583
        # while is_trigger_available() gives True
        n_max_modes = 20
        modes_c = c_char_p()
        n_modes = self._dll_get_trigger_modes(
            self.__handle, modes_c, n_max_modes
            )
        # TODO: unclear how to read the modes.
        print(modes_c)
        # return [mi.value.decode() for mi in modes_c if mi.value]

    _dll_reset = _dll.IC_ResetProperties
    _dll_reset.argtypes = (GrabberHandlePtr,)
    _dll_reset.errcheck = check_dll_return()

    def reset(self):
        """Reset the camera to default properties."""
        self._dll_reset(self.__handle)


class TIS_CAM:

    def GetImage(self):
        BildDaten = self.get_image_description()[:4]
        lWidth = BildDaten[0]
        lHeight = BildDaten[1]
        iBitsPerPixel = BildDaten[2] // 8 # // 8 was commented out??

        buffer_size = lWidth*lHeight*iBitsPerPixel*c_sizeof(c_uint8)
        img_ptr = self.get_image_ptr()

        Bild = c_cast(img_ptr, POINTER(c_ubyte * buffer_size))


        img = np.ndarray(buffer = Bild.contents,
                     dtype = np.uint8,
                     shape = (lHeight,
                              lWidth,
                              iBitsPerPixel))
        return img

    def GetImageEx(self):
        """ Return a numpy array with the image data tyes
        If the sink is Y16 or RGB64 (not supported yet), the dtype in the array is uint16, otherwise it is uint8
        """
        BildDaten = self.get_image_description()[:4]
        lWidth = BildDaten[0]
        lHeight = BildDaten[1]
        iBytesPerPixel = BildDaten[2] // 8
        buffer_size = lWidth*lHeight*iBytesPerPixel*c_sizeof(c_uint8)
        img_ptr = self.get_image_ptr()

        Bild = c_cast(img_ptr, POINTER(c_ubyte * buffer_size))
        pixeltype = np.uint8
        if  BildDaten[3] == 4: #SinkFormat.Y16:
            pixeltype = np.uint16
            iBytesPerPixel = 1

        # Look at numpy.ctypeslib.as_array(obj, shape=None)
        img = np.ndarray(buffer = Bild.contents,
                     dtype = pixeltype,
                     shape = (lHeight,
                              lWidth,
                              iBytesPerPixel))
        return img
