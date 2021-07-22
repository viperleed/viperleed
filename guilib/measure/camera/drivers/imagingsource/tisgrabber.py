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
    if isinstance(string, str)
        return string.encode()
    raise TypeError(f"Invalid string type {type(string).__name__!r}. "
                    "Expected, 'str', 'bytes', or 'bytearray")


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
            raise RuntimeError("Camera should never return "
                               "SinkFormat.MEGA")
        return sink_format


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



# TODO: look at the documentation at https://github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/Python/Open%20Camera%2C%20Grab%20Image%20to%20OpenCV
# TODO: also at (ip-config) https://github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/c%23/Open%20GigE%20Camera%20by%20IP%20Address

# TODO: snake_case, cleanup. Take also a look at the github
# TODO: try removing some of the many functions, especially those we will never use

class GrabberDLL:
    """Interface class that defines python equivalents to the DLL calls."""

    # ctypes.POINTER returns a class
    GrabberHandlePtr = POINTER(GrabberHandle)

    # To prevent creating bindings to the dll functions each time
    # methods are called, each method that calls a dll function
    # has ctypes definition of the function used right above it.
    # Note: ctypes assumes restype == c_int. Thus restype is
    # omitted below where this is the case.

    if sys.maxsize > 2**32:
        __dll = WinDLL("tisgrabber_x64.dll")
    else:
        __dll = WinDLL("tisgrabber.dll")

    __initalized = False

    def __init__(self):
        """Initialize grabber instance."""
        self.__handle = self.create_grabber()
        self.__available_cam_properties = {}

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
            __dll_init = cls.__dll.IC_InitLibrary
            __dll_init.argtypes = c_char_p
            __dll_init.errcheck = check_dll_return()

            __dll_init(_to_bytes(license_key))
            cls.__initalized = True

    # Call init_library with None straight away. It was done the
    # same way in the original tisgrabber.py. This guarantees it
    # is called only once the first time the class is created.
    init_library()

    __dll_create_grabber = __dll.IC_CreateGrabber
    __dll_create_grabber.restype = GrabberHandlePtr
    __dll_create_grabber.argtypes = None
    __dll_create_grabber.errcheck = check_dll_return()

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
        return self.__dll_create_grabber()

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

    __dll_get_frame_rate =  __dll.IC_GetFrameRate
    __dll_get_frame_rate.restype = c_float
    __dll_get_frame_rate.argtypes = (GrabberHandlePtr,)
    __dll_get_frame_rate.errcheck = check_dll_return('>0')

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
        return self.__dll_get_frame_rate(self.__handle).value

    __dll_set_frame_rate = __dll.IC_SetFrameRate
    __dll_set_frame_rate.argtypes = GrabberHandlePtr, c_float
    __dll_set_frame_rate.errcheck = check_dll_return(
        # NOT_AVAILABLE would be an error, but NOT_IN_LIVEMODE
        # seems more reasonable to return. Perhaps better check
        # live mode explicitly and leave NOT_AVAILABLE.
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
        self.__dll_set_frame_rate(self.__handle, rate)

    # TODO: unclear what 'color format' means. Is it == video_format??
    #       is it == sink_type ??
    # TODO: eventually it would be nice to align the returns to the
    #       signature I want in CameraABC (i.e., return the no. of
    #       color channels rather than the format)
    __dll_get_image_description = __dll.IC_GetImageDescription
    __dll_get_image_description.argtypes = (
        GrabberHandlePtr, POINTER(c_long), POINTER(c_long),
        POINTER(c_int), POINTER(c_int)
        )
    __dll_get_image_description.errcheck = check_dll_return()

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
        width, height = c_long(), c_long()
        bits_per_pixel, color_format_c = c_int(), c_int()

        self.__dll_get_image_description(self.__handle, width, height,
                                         bits_per_pixel, color_format_c)
        bytes_per_pixel = bits_per_pixel.value // 8
        try:
            color_format = SinkFormat.get(color_format_c.value)
        except ValueError as err:
            raise ImagingSourceError(err.args[0] + " received.") from err

        return width.value, height.value, bytes_per_pixel, color_format

    __dll_n_devices = __dll.IC_GetDeviceCount
    __dll_n_devices.argtypes = None
    __dll_n_devices.errcheck = check_dll_return('>=0')

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
        return self.__dll_n_devices()

    __dll_get_sink_format = __dll.IC_GetFormat
    __dll_get_sink_format.argtypes = (GrabberHandlePtr,)

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
        # TODO: needs live mode?
        # TODO: old implementation was returning RGB24 in case
        #       of 'unknown' format
        sink_format_c = self.__dll_get_sink_format(self.__handle)
        try:
            sink_format = SinkFormat.get(sink_format_c.value)
        except ValueError as err:
            raise ImagingSourceError(err.args[0] + " received")
        except ImagingSourceError as err:
            raise ImagingSourceError("Camera should be in live mode. Call "
                                     ".start_live(), then retry.") from err

        return sink_format

    __dll_set_sink_format = __dll.IC_SetFormat
    __dll_set_sink_format.argtypes = GrabberHandlePtr, c_int
    __dll_set_sink_format.errcheck = check_dll_return()

    @sink_format.setter
    def sink_format(self, sink_format):
        """Set the sink type.

        This method must be called before images can be snapped.
        The sink type describes the format of the buffer where
        snapped images are stored. Notice that the sink type may
        differ from the video format.

        Parameters
        ----------
        sink_type : str, int, or SinkFormat
            One of the values in SinkFormat. Note that UYVY can
            only be used in conjunction with a UYVY video format
        """
        sink_format = SinkFormat.get(sink_format)
        self.__dll_set_sink_format(self.__handle, sink_format.value)

    __dll_is_device_valid = __dll.IC_IsDevValid
    __dll_is_device_valid.argtypes = (GrabberHandlePtr,)
    __dll_is_device_valid.errcheck = check_dll_return()

    @property
    def valid(self):  # TODO: maybe I don't need this public
        """Return whether there is a valid open device."""
        # TODO: after checking possible returns, get rid of
        # the error checker if the only thing that can happen
        # is a generic ERROR. The purpose of the error checker
        # is just printing out the return values.
        try:
            ret_val = self.__dll_is_device_valid(self.__handle)
        except ImagingSourceError:
            # No device open
            pass
        return  ret_val == DLLReturns.SUCCESS.value

    __dll_video_format_width = __dll.IC_GetVideoFormatWidth
    __dll_video_format_height = __dll.IC_GetVideoFormatHeight
    __dll_video_format_width.argtypes = (GrabberHandlePtr,)
    __dll_video_format_height.argtypes = (GrabberHandlePtr,)
    __dll_video_format_width.errcheck = check_dll_return()
    __dll_video_format_height.errcheck = check_dll_return()

    @property
    def video_format_shape(self):
        """Return (width, height) (pixels) of the current video format."""
        width = self.__dll_video_format_width(self.__handle)
        height = self.__dll_video_format_height(self.__handle)
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

    __dll_set_video_format = __dll.IC_SetVideoFormat
    __dll_set_video_format.argtypes = GrabberHandlePtr, c_char_p
    ____dll_set_video_format.errcheck = check_dll_return()

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
        self.__dll_set_video_format(self.__handle,
                                    _to_bytes(video_format))

    # TODO: check which errors
    __dll_open_by_unique_name = __dll.IC_OpenDevByUniqueName
    __dll_open_by_unique_name.argtypes = GrabberHandlePtr, c_char_p
    __dll_open_by_unique_name.errcheck = check_dll_return()

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
        self.__dll_open_by_unique_name(self.__handle,
                                       _to_bytes(unique_name))

    # TODO: may it return other errors, like NO_HANDLE and NO_DEVICE?
    __dll_snap_image = __dll.IC_SnapImage
    __dll_snap_image.argtypes = GrabberHandlePtr, c_int
    __dll_snap_image.errcheck = check_dll_return(
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
        self.__dll_snap_image(self.__handle, timeout)

    __dll_start_live = __dll.IC_StartLive
    __dll_start_live.argtypes = GrabberHandlePtr, c_int
    __dll_start_live.errcheck = check_dll_return()

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
        self.__dll_start_live(self.__handle, bool(show_video))

    # TODO: unclear which errors it may return
    __dll_stop_live = __dll.IC_StopLive
    __dll_stop_live.argtypes = (GrabberHandlePtr,)
    __dll_stop_live.errcheck = check_dll_return()

    def stop_live(self):
        """Stop live mode."""
        self.__dll_stop_live(self.__handle)

    __dll_show_device_selection_dialog = __dll.IC_ShowDeviceSelectionDialog
    __dll_show_device_selection_dialog.restype = GrabberHandlePtr
    __dll_show_device_selection_dialog.argtypes = (GrabberHandlePtr,)

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
        self.__handle = self.__dll_show_device_selection_dialog(self.__handle)

    # TODO: restype did not match signature. Also, the old python
    # implementation assigned the return to __handle. Should be
    # checked.
    # TODO: implemented just for testing. will be later removed.
    __dll_show_property_dialog = __dll.IC_ShowPropertyDialog
    __dll_show_property_dialog.restype = c_int  # was GrabberHandlePtr
    __dll_show_property_dialog.argtypes = (GrabberHandlePtr,)
    __dll_show_property_dialog.errcheck = check_dll_return()

    def temp_property_dialog(self):
        """Show the VCDProperty dialog.

        Requires a video-capture device to be open.

        Returns
        -------
        ret_val : int
            SUCCESS on success, NO_HANDLE if grabber_handle is invalid,
            NO_DEVICE if no device was open, ERROR otherwise.
        """
        self.__dll_show_property_dialog(self.__handle)

    __dll_unique_name_from_index = __dll.IC_GetUniqueNamefromList
    __dll_unique_name_from_index.restype = c_char_p
    __dll_unique_name_from_index.argtypes = (c_int,)
    __dll_unique_name_from_index.errcheck = check_dll_return('pointer')

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
        name = self.__dll_unique_name_from_index(index)
        if name is not None:
            return name.decode()
        return ''

    __dll_video_format_by_index = __dll.IC_GetVideoFormat
    __dll_video_format_by_index.restype = c_char_p
    __dll_video_format_by_index..argtypes = GrabberHandlePtr, c_int
    __dll_video_format_by_index.errcheck = check_dll_return('pointer')

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
        video_format = self.__dll_video_format_by_index(self.__handle, index)
        if video_format:
            return video_format
        return ''

    __dll_n_video_formats = __dll.IC_GetVideoFormatCount
    __dll_n_video_formats.argtypes = (GrabberHandlePtr,)
    __dll_n_video_formats.errcheck = check_dll_return('>=0')

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
        return self.__dll_n_video_formats(self.__handle)

    __dll_get_image_ptr = __dll.IC_GetImagePtr
    __dll_get_image_ptr.restype = c_void_p  # void to also allow None
    __dll_get_image_ptr.argtypes = (GrabberHandlePtr,)
    __dll_get_image_ptr.errcheck = check_dll_return('pointer')

    def __get_image_ptr(self):
        """Get a pointer to the pixel data of the last snapped image.

        Returns
        -------
        pointer : ctypes.c_void_p
            A pointer to the the first byte of the bottommost line
            of the image (images are saved from bottom to top!).
            Returns None in case of error.
        """
        return self.__dll_get_image_ptr(self.__handle)

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
            available = self.__dll_is_cam_property_available(
                self.__handle, cam_property.value
                )
            self.__available_cam_properties[cam_property.name] = (
                available.value == DLLReturns.SUCCESS.value
                )
        return self.__available_cam_properties

    # TODO: unclear what it returns exactly.
    # TODO: check also whether setting exposure this way is the
    #       same as using set_property_value_int and absolute value
    #       If it is the same, probably better keep only one to
    #       limit overcrowding of class.
    __dll_get_camera_property = __dll.IC_GetCameraProperty
    __dll_get_camera_property.argtypes = (GrabberHandlePtr, c_int,
                                          POINTER(c_long))
    __dll_get_camera_property.errcheck = check_dll_return()

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
        self.__dll_get_camera_property(self.__handle, cam_property.value, prop)
        return prop.value

    __dll_set_camera_property = __dll.IC_SetCameraProperty
    __dll_set_camera_property.argtypes = GrabberHandlePtr, c_int, c_long
    __dll_set_camera_property.errcheck = check_dll_return()

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
        self.__dll_set_camera_property(self.__handle, cam_property.value,
                                       value)

    # TODO: check return values
    __dll_get_camera_property_range = __dll.IC_CameraPropertyGetRange
    __dll_get_camera_property_range.argtypes = (GrabberHandlePtr, c_int,
                                                POINTER(c_long),
                                                POINTER(c_long))
    __dll_get_camera_property_range.errcheck = check_dll_return()

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
            raise ImagingSourceError(f"Property {cam_property.name} not available")
        prop_min, prop_max = c_long(), c_long()

        self.__dll_get_camera_property(self.__handle, cam_property.value,
                                       prop_min, prop_max)
        return prop_min, prop_max

    __dll_is_property_available =  __dll.IC_IsPropertyAvailable
    __dll_is_property_available.argtypes = (GrabberHandlePtr, c_char_p,
                                            c_char_p)
    __dll_is_property_available.errcheck = check_dll_return()

    def is_property_available(self, property_name, element_name=None):
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
        available = self.__dll_is_property_available(
            self.__handle, _to_bytes(property_name), _to_bytes(element_name)
            )
        return available.value == DLLReturns.SUCCESS.value

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
    __dll_get_property_value_float = __dll.IC_GetPropertyAbsoluteValue
    __dll_get_property_value_float.errcheck = check_dll_return()
    __dll_get_property_value_float.argtypes = (
        GrabberHandlePtr, c_char_p,c_char_p, c_float_p
        )
    __dll_get_property_value_range_float = IC_GetPropertyAbsoluteValueRange
    __dll_get_property_value_range_float.errcheck = check_dll_return()
    __dll_get_property_value_range_float.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_float_p, c_float_p
        )
    __dll_set_property_value_float = __dll.IC_SetPropertyAbsoluteValue
    __dll_set_property_value_float.errcheck = check_dll_return()
    __dll_set_property_value_float.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_float
        )

    # Value as int: getter, range getter, setter.
    __dll_get_property_value_int = __dll.IC_GetPropertyValue
    __dll_get_property_value_int.errcheck = check_dll_return()
    __dll_get_property_value_int.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, POINTER(c_long)
        )
    __dll_get_property_value_range_int = __dll.IC_GetPropertyValueRange
    __dll_get_property_value_range_int.errcheck = check_dll_return()
    __dll_get_property_value_range_int.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, POINTER(c_long), POINTER(c_long)
        )
    __dll_set_property_value_int = __dll.IC_SetPropertyValue
    __dll_set_property_value_int.errcheck = check_dll_return()
    __dll_set_property_value_int.argtypes = (
        GrabberHandlePtr, c_char_p, c_char_p, c_int
        )

    # On/Off (== bool): getter, setter. [NO RANGE OBVIOUSLY]
        __dll_get_on_off_property = __dll.IC_GetPropertySwitch
    __dll_get_on_off_property.errcheck = check_dll_return()
    __dll_get_on_off_property.argtypes= (
        GrabberHandlePtr, c_char_p, c_char_p, POINTER(c_long)
        )
    __dll_set_on_off_property = __dll.IC_SetPropertySwitch
    __dll_set_on_off_property.errcheck = check_dll_return()
    __dll_set_on_off_property.argtypes= (
        GrabberHandlePtr, c_char_p, c_char_p, c_int
        )

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
            getter = self.__dll_get_property_value_float
            property_value = c_float()
        elif method is int:
            getter = self.__dll_get_property_value_int
            property_value = c_long()
        elif method is bool:
            getter = self.__dll_get_on_off_property
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
            setter = self.__dll_set_property_value_float
            value = float(value)
        elif method is int:
            setter = self.__dll_set_property_value_int
            value = int(value)
        elif method is bool:
            setter = self.__dll_set_on_off_property
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
            range_getter = self.__dll_get_property_value_range_float
            prop_min, prop_max = c_float(), c_float()
        elif method is int:
            range_getter = self.__dll_get_property_value_range_int
            prop_min, prop_max = c_long(), c_long()
        else:
            raise ValueError(
                f"Unexpected method {method} for retrieving the range "
                f"of VCD property/element {property_name}/{element_name}"
                )
        range_getter(self.__handle, _to_bytes(property_name),
                     _to_bytes(element_name), prop_min, prop_max)
        return prop_min.value, prop_max.value

    __dll_trigger_property_once = __dll.IC_PropertyOnePush
    __dll_trigger_property_once.argtypes = (GrabberHandlePtr, c_char_p,
                                                c_char_p)
    __dll_trigger_property_once.errcheck = check_dll_return()

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
        self.__dll_trigger_property_once.(
            self.__handle, _to_bytes(property_name), _to_bytes(element_name)
            )

    # Prototype of frame-ready callbacks:
    FrameReadyCallBackType = CFUNCTYPE(
        c_void_p,          # return type
        GrabberHandlePtr,  # HGRABBER; funnily was c_int!
        POINTER(c_ubyte),  # pointer to first byte of image data
        c_ulong,           # current frame number
        py_object          # other data (python) passed along
        )
    
    # A little trick to make functions decorated with
    # @FrameReadyCallBackType inspect-able with inspect.signature
    # NB: other inspect features may not work, as there is no
    # explicit support for CFUNCTYPE() objects.
    def __example_cb(__handle, img_buffer, frame_no, py_obj_for_callback):
        pass
    FrameReadyCallBackType.__signature__ = (
        inspect.Signature.from_callable(__example_cb)
        )

    # The py_object at the end is the same that will later be passed
    # on to the frame-ready callback as the last argument.
    __dll_set_frame_ready_callback = __dll.IC_SetFrameReadyCallback
    __dll_set_frame_ready_callback.errcheck = check_dll_return()
    __dll_set_frame_ready_callback.argtypes = (
        GrabberHandlePtr, FrameReadyCallBackType, py_object
        )
    
    def set_frame_ready_callback(self, on_frame_ready, py_obj_for_callback):
        """Define the frame-ready callback function.
        
        Parameters
        ----------
        on_frame_ready : callable
            Should be decorated with @FrameReadyCallBackType to ensure
            proper typing. This is the function that will be called
            whenever a new frame has been acquired (in 'continuous'
            mode) or whenever an image is 'snapped'. The callable
            should have signature:
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
            with @FrameReadyCallBackType.
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
                f"@{self.__class__.__name__}.FrameReadyCallBackType. This "
                " is needed to ensure appropriate type checking/conversions"
                )
        self.__dll_set_frame_ready_callback(self.__handle, on_frame_ready,
                                            py_obj_for_callback)

    SetContinuousMode = __dll.IC_SetContinuousMode

    # TODO: signature was wrong: c_void_p instead of GrabberHandlePtr
    OpenVideoCaptureDevice = __dll.IC_OpenVideoCaptureDevice
    OpenVideoCaptureDevice.restype = c_int
    OpenVideoCaptureDevice.argtypes = GrabberHandlePtr, c_char_p

    # TODO: signature was wrong: c_void_p instead of GrabberHandlePtr
    CloseVideoCaptureDevice = __dll.IC_CloseVideoCaptureDevice
    CloseVideoCaptureDevice.restype = c_int
    CloseVideoCaptureDevice.argtypes = GrabberHandlePtr, c_char_p

# ############################################################################

    IsTriggerAvailable = __dll.IC_IsTriggerAvailable
    IsTriggerAvailable.restype = c_int
    IsTriggerAvailable.argtypes = (GrabberHandlePtr,)

    EnableTrigger = __dll.IC_EnableTrigger
    EnableTrigger.restype = c_int
    EnableTrigger.argtypes = GrabberHandlePtr, c_int

    SoftwareTrigger = __dll.IC_SoftwareTrigger
    SoftwareTrigger.restype = c_int
    SoftwareTrigger.argtypes = (GrabberHandlePtr,)

# ############################################################################

    ResetProperties = __dll.IC_ResetProperties
    ResetProperties.restype = c_int
    ResetProperties.argtypes = (GrabberHandlePtr,)

    GetTriggerModes = __dll.IC_GetTriggerModes
    GetTriggerModes.restype = c_int
    GetTriggerModes.argtypes = (GrabberHandlePtr, c_char_p,
                                c_int)

# ############################################################################


# TODO: do error checking either with a @decorator or with
# errcheck from ctypes

class TIS_CAM:
    @property
    def callback_registered(self):
        return self._callback_registered

    def __init__(self):
        self._handle = GrabberDLL.create_grabber()
        self._callback_registered = False
        self._frame = {'num'    :   -1,
                       'ready'  :   False}


    def SetFrameReadyCallback(self, CallbackFunction, data):
        """ Set a callback function, which is called, when a new frame arrives.

        CallbackFunction : The callback function

        data : a self defined class with user data.
        """
        return GrabberDLL.SetFrameReadyCallback( self._handle, CallbackFunction, data )

    def SetContinuousMode(self, Mode):
        ''' Determines, whether new frames are automatically copied into memory.

        :param Mode: If 0, all frames are copied automatically into memory. This is recommened, if the camera runs in trigger mode.
                      If 1, then snapImages must be called to get a frame into memory.
        :return: None
        '''
        return GrabberDLL.SetContinuousMode(self._handle, Mode)

    def open(self, unique_device_name):
        """ Open a device

        unique_device_name : The name and serial number of the device to be opened. The device name and serial number are separated by a space.
        """
        test = GrabberDLL.open_by_unique_name(self._handle,
                                              self.s(unique_device_name))

        return test

    def save_device_state_to_file(self, FileName):
        return GrabberDLL.save_device_state_to_file(self._handle, self.s(FileName))

    def load_device_state_from_file(self,FileName):
        self._handle = GrabberDLL.load_device_state_from_file(self._handle,self.s(FileName))

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

    def PropertyOnePush(self, Property, Element ):
        error = GrabberDLL.PropertyOnePush(self._handle,
                                                self.s(Property),
                                                self.s(Element ))
        return error

    def SetPropertyAbsoluteValue(self, Property, Element, Value ):
        error = GrabberDLL.SetPropertyAbsoluteValue(self._handle,
                                                self.s(Property),
                                                self.s(Element),
                                                Value)
        return error

    def GetPropertyAbsoluteValue(self, Property, Element,Value ):
        """ Get a property value of absolute values interface, e.g. seconds or dB.
        Example code:
        ExposureTime=[0]
        Camera.GetPropertyAbsoluteValue("Exposure","Value", ExposureTime)
        print("Exposure time in secods: ", ExposureTime[0])

        :param Property: Name of the property, e.g. Gain, Exposure
        :param Element: Name of the element, e.g. "Value"
        :param Value: Object, that receives the value of the property
        :returns: 0 on success
        """
        lValue = c_float()
        error = GrabberDLL.GetPropertyAbsoluteValue(self._handle,
                                                self.s(Property),
                                                self.s(Element),
                                                lValue)
        Value[0] = lValue.value
        return error

    def openVideoCaptureDevice(self, DeviceName):
        ''' Open the device specified by DeviceName
        :param DeviceName: Name of the device , e.g. "DFK 72AUC02"
        :returns: 1 on success, 0 otherwise.
        '''
        return GrabberDLL.OpenVideoCaptureDevice(self._handle, self.s(DeviceName))

    def enableAutoCameraProperty(self, property, onoff):
        return GrabberDLL.EnableAutoCameraProperty(self._handle,property,onoff)

    def enableVideoAutoProperty(self, property, onoff):
        return GrabberDLL.EnableVideoAutoProperty(self._handle,property,onoff)

    def IsTriggerAvailable(self):
        error = GrabberDLL.IsTriggerAvailable(self._handle)

        return error

    def EnableTrigger(self, onoff):
        error = GrabberDLL.EnableTrigger(self._handle, onoff)

        return error

    def SoftwareTrigger(self):
        error = GrabberDLL.SoftwareTrigger(self._handle)

        return error

    def ResetProperties(self):

        error = GrabberDLL.ResetProperties(self._handle)

        return error
#        def GetTriggerModes(self, ModeList, Size):
#            error = GrabberDLL.GetTriggerModes(self._handle,
#                                                   ModeList,
#                                                   Size)
#            return error

        #    int AC IC_GetTriggerModes( HGRABBER hGrabber,  char *szModeList, int iSize  );