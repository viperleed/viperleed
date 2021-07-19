"""Module grabber of camera.drivers.imagingsource.

========================================
   ViPErLEED Graphical User Interface
========================================

This is an extension of the original tisgrabber.py module available
in the original form at:
github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/
          Python/Open%20Camera%2C%20Grab%20Image%20to%20OpenCV

Edits include code and code-style improvements.

Author/editor: Michele Riva
Date edits started: 2020-07-18

The original tisgrabber.py module was:
Created on Mon Nov 21 09:44:40 2016
@author: Daniel Vassmer, Stefan_Geissler
"""
from enum import Enum
import ctypes as C
import os
from pathlib import Path
import sys
import numpy as np

dll_path = Path(__file__).parent.resolve()
print(f"{dll_path=}")
try:
    os.add_dll_directory(dll_path)
except AttributeError:
    pass
sys.path.append(dll_path)


class SinkFormats(Enum):
    Y800 = 0    # monochrome 1-byte;    top-left
    RGB24 = 1   # 3 bytes: B G R;       top-left
    RGB32 = 2   # 4 bytes: B G R alpha; top-left; alpha unused == 0
    UYVY = 3    # 4-bytes color;        top-left; see fourcc.org; bit complicated (U&V shared by neighbor pixels, y for both)
    Y16 = 4     # gray 2-bytes;         top-left
    NONE = 5    # used only as return value
    MEGA = 65536  # Borland C++ 6 compatibility
    # RGB64     # 8 bytes: 2B 2G 2R 2A; top-left; seems unsupported


class CameraProperty(Enum):
    PAN = 0
	TILT = 1
	ROLL = 2
	ZOOM = 3
	EXPOSURE = 4
	IRIS = 5
	FOCUS = 6
	MEGA = 65536  # Borland C++ 6 compatibility

    @classmethod
	def from_value_or_name(cls, cam_property):
        """Return an appropriate enum.

        Parameters
        ----------
        cam_property : str, int, or CameraProperty
            Some reference to a camera property

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
                cam_property = getattr(cls, cam_property)
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


class DLLReturns(Enum):
    """Enum of return values from the firmware DLL.

    Attributes
    ----------
    SUCCESS
        Function has been performed without an error.
    ERROR
        Function encountered an unspecified error.
    NO_HANDLE
        Grabber handle was not created before calling a
        function that needs it. Call create_grabber().
    NO_DEVICE
        Function requires an open device, but no device is open.
        Call OpenVideoCaptureDevice().
    NOT_AVAILABLE
        The specified property is not available in the device
    NO_PROPERTYSET
        The property set was not queried for the current grabber
        handle. Please check, whether IC_QueryPropertySet() was called
        once before using the function.
    DEFAULT_WINDOW_SIZE_SET
        Setting of a custom live display window size failed, because
        IC_SetDefaultWindowPosition(false) was not called
    NOT_IN_LIVEMODE
        Live mode is needed, but the video-capture device is not in
        live mode. Call start_live().
    PROPERTY_ITEM_NOT_AVAILABLE
        The device does not support the requested property, or the name
        of a property is wrong.
    PROPERTY_ELEMENT_NOT_AVAILABLE
        The device does not support the requested property, or the name
        of a property is wrong.
    PROPERTY_ELEMENT_WRONG_INTERFACE
        The property element selected exists, but it does not have the
        requested attribute (e.g., 'Exposure' 'Auto' has no 'range')
    INDEX_OUT_OF_RANGE
        The index passed is not in the expected range.
    WRONG_XML_FORMAT
        The XML file given does not contain data, or has invalid format
    WRONG_INCOMPATIBLE_XML
        The XML file given does not contain compatible XML data
    NOT_ALL_PROPERTIES_RESTORED
        Camera was opened, but some of the properties for which a new
        setting was passed was not accepted.
    DEVICE_NOT_FOUND
        The device specified in the XML was not found. E.g., same model
        but different serial number, or no camera connected at all.
    FILE_NOT_FOUND
        The file path given does not exist
    """
    SUCCESS = 1
    ERROR = 0
    NO_HANDLE = -1
    NO_DEVICE = -2
    NOT_AVAILABLE = -3
    NO_PROPERTYSET = -3
    DEFAULT_WINDOW_SIZE_SET = -3
    NOT_IN_LIVEMODE = -3
    PROPERTY_ITEM_NOT_AVAILABLE = -4
    PROPERTY_ELEMENT_NOT_AVAILABLE = -5
    PROPERTY_ELEMENT_WRONG_INTERFACE = -6
    INDEX_OUT_OF_RANGE = -7
    WRONG_XML_FORMAT = -1
    WRONG_INCOMPATIBLE_XML = -3
    NOT_ALL_PROPERTIES_RESTORED = -4
    DEVICE_NOT_FOUND = -5
    FILE_NOT_FOUND = 35


ImageFileTypes = {'BMP':0, 'JPEG':1}


class GrabberHandle(C.Structure):
    _fields_ = [('unused', C.c_int)]



# TODO: look at the documentation at https://github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/Python/Open%20Camera%2C%20Grab%20Image%20to%20OpenCV
# TODO: also at (ip-config) https://github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/c%23/Open%20GigE%20Camera%20by%20IP%20Address

# TODO: snake_case, cleanup. Take also a look at the github
# TODO: try out the .errcheck of ctypes
# TODO: try removing some of the many functions, especially those we will never use

class GrabberDLL:
    """Interface class that defines python equivalents to the DLL calls."""

    # ctypes.POINTER returns a class
    GrabberHandlePtr = C.POINTER(GrabberHandle)

    # To prevent create bindings to the dll functions each time
    # methods are called, each method that calls a dll function
    # has ctypes definition of the function used right above it.

    if sys.maxsize > 2**32 :
        __dll = C.windll.LoadLibrary("tisgrabber_x64.dll")
    else:
        __dll = C.windll.LoadLibrary("tisgrabber.dll")

    # def __init__(self, **keyargs):
        # """Initialize the Albatross from the keyword arguments."""
        # self.__dict__.update(keyargs)

    def __init__(self):
        """Initialize grabber instance."""
        self.__handle = self.create_grabber()

    @staticmethod
    def __to_bytes(string):
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

    def init_library(self, license_key):
        """Initialize the IC Imaging Control dll library.

        This function must be called only once before any other
        functions of the dll are called.

        Parameters
        ----------
        license_key : str or None
            The license key. Pass None if using a 'trial' version.

        Returns
        -------
        ret_val : int
            SUCCESS on success, ERROR if license key is wrong or for
            other generic errors.
        """
        return self.__dll.IC_InitLibrary(self.__to_bytes(license_key))

    # Call init_library with None straight away. It was done the
    # same way in the original tisgrabber.py
    init_library(None)

    __dll_create_grabber = __dll.IC_CreateGrabber
    __dll_create_grabber.restype = GrabberHandlePtr
    __dll_create_grabber.argtypes = None

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

    __dll_n_devices = __dll.IC_GetDeviceCount
    __dll_n_devices.restype = C.c_int
    __dll_n_devices.argtypes = None

    # TODO: an alternative doing the loop in C would be to use IC_ListDevices
    # but it is unclear whether it would return unique names or just device
    # names. Perhaps worth implementing and testing. Then we would not need
    # n_devices and unique_name_from_index.
    @property
    def devices(self):
        """Return the unique names of available devices.

        Returns
        -------
        unique_names : list
            List of unique device names, each as bytes (UTF-8).
        """
        return [self.unique_name_from_index(i)
                for i in range(self.n_devices)]

    __dll_get_frame_rate =  __dll.IC_GetFrameRate
    __dll_get_frame_rate.restype = C.c_float
    __dll_get_frame_ratee.argtypes = (GrabberHandlePtr,)

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
        # TODO: check it is ok. The original implementation was first
        # defining frame_rate_c = ctypes.c_long(), then overwriting
        # it with the return value
        frame_rate_c = self.__dll_get_frame_rate(self.__handle)
        return frame_rate_c.value

    __dll_set_frame_rate = __dll.IC_SetFrameRate
    __dll_set_frame_rate.restype = C.c_int
    __dll_set_frame_rate.argtypes = GrabberHandlePtr, C.c_float

    @frame_rate.setter
    def frame_rate(self, rate):
        """Set a new frame rate.

        A new frame rate can not be set while in live mode. Call
        stop_live(), then try again.

        Parameters
        ----------
        frame_rate : float
            The new frame rate.

        Returns
        -------
        ret_val : DLLReturns.value
            SUCCESS, NOT_AVAILABLE (not supported by the device),
            NO_HANDLE (invalid grabber handle), NO_DEVICE (no device
            opened), NOT_IN_LIVEMODE (device is in live mode).
        """
        # TODO: if live, stop, then set, then restart
        return self.__dll_set_frame_rate(self.__handle, rate)

    # TODO: unclear what 'color format' means. Is it == video_format??
    #       is it == sink_type ??
    # TODO: eventually it would be nice to align the returns to the
    #       signature I want in CameraABC (i.e., return the no. of
    #       color channels rather than the format)
    __dll_get_image_description = __dll.IC_GetImageDescription
    __dll_get_image_description.restype = C.c_int
    __dll_get_image_description.argtypes = (GrabberHandlePtr,
                                            C.POINTER(C.c_long),
                                            C.POINTER(C.c_long),
                                            C.POINTER(C.c_int),
                                            C.POINTER(C.c_int))

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
        format : SinkFormats
            The current color format

        C-Returns
        ---------
        ret_val : int
            SUCCESS on success, ERROR if something went wrong.
        """
        width = C.c_long()
        height = C.c_long()
        bits_per_pixel = C.c_int()
        color_format_c = C.c_int()

        self.__dllget_image_description(self.__handle, width, height,
                                        bits_per_pixel, color_format_c)
        bytes_per_pixel = bits_per_pixel.value // 8
        try:
            color_format = SinkFormats(color_format_c.value)
        except ValueError as err:
            raise RuntimeError("Unknown color format identifier "
                               f"{color_format_c} received") from err

        return width.value, height.value, bytes_per_pixel, color_format

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

    # TODO: unclear what this returns in case of error
    __dll_get_sink_format = __dll.IC_GetFormat
    __dll_get_sink_format.restype = C.c_int
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
            One of the values in SinkFormats. Returns !!UNCLEAR WHAT!!
            if not in live mode or if an error occurred.
        """
        # TODO: needs live mode?
        # TODO: old implementation was returning RGB24 in case
        #       of 'unknown' format
        sink_format_c = self.__dll_get_sink_format(self.__handle)
        try:
            sink_format = SinkFormats(sink_format_c.value)
        except ValueError as err:
            raise RuntimeError("Camera returned unknown image "
                               f"format {sink_format_c.value}")
        return sink_format

    __dll_set_sink_format = __dll.IC_SetFormat
    __dll_set_sink_format.restype = C.c_int
    __dll_set_sink_format.argtypes = GrabberHandlePtr, C.c_int

    @sink_format.setter
    def sink_format(self, sink_format):
        """Set the sink type.

        This method must be called before images can be snapped.
        The sink type describes the format of the buffer where
        snapped images are stored. Notice that the sink type may
        differ from the video format.

        Parameters
        ----------
        sink_type : SinkFormats
            One of the values in SinkFormats. Note that UYVY can
            only be used in conjunction with a UYVY video format

        Returns
        -------
        ret_val : int
            SUCCESS on success, ERROR if something went wrong.
        """
        try:
            sink_format = SinkFormats(sink_format)
        except ValueError as err:
            raise ValueError(f"Unknown sink format {sink_format}") from err
        return self.__dll_set_sink_format(self.__handle, sink_format.value)

    __dll_is_device_valid = __dll.IC_IsDevValid
    __dll_is_device_valid.restype = C.c_int
    __dll_is_device_valid.argtypes = (GrabberHandlePtr,)

    # TODO: check errors, make the return a bool
    @property
    def valid(self):  # TODO: maybe I don't need this public
        """Return whether there is a valid open device.

        Returns
        -------
        ret_val : int
            SUCCESS if a valid device is open, ERROR otherwise.
        """
        return sef.__dll_is_device_valid(self.__handle)

    @property
    def video_format_size(self):
        """Return (width, height) (pixels) of the current video format."""
        return self.__video_format_width(), self.__video_format_height()

    # TODO: perhaps would be nice to run the loop at the c level. This
    # can be done with IC_ListVideoFormats. Then one would not need
    # __n_video_formats() nor __video_format_by_index(idx)
    @property
    def video_formats(self):
        """Return a list of the video formats available for a device."""
        n_formats = self.__n_video_formats()
        return [self.__video_format_by_index(i)
                for i in range(self.__n_video_formats())]

    __dll_set_video_format = __dll.IC_SetVideoFormat
    __dll_set_video_format.restype = C.c_int
    __dll_set_video_format.argtypes = GrabberHandlePtr, C.c_char_p

    # TODO: funnily enough there is no way to get the video format
    # back as a string. The only way to fetch info about the
    # video format if via get_image_description()
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
        return self.__dll_set_video_format(self.__handle,
                                           self.__to_bytes(video_format))

    __dll_open_by_unique_name = __dll.IC_OpenDevByUniqueName
    __dll_open_by_unique_name.restype = C.c_int
    __dll_open_by_unique_name.argtypes = GrabberHandlePtr, C.c_char_p

    def open_by_unique_name(self, unique_name):
        """Open a video capture by using its unique name.

        Use IC_GetUniqueName() to retrieve the unique name of an
        open device. unique_name_from_index(index) can be used
        before the device is open to retrieve the unique name.

        Parameters
        ----------
        unique_name : bytes
            Unique name of the device to be opened

        Returns
        -------
        ret_val : int
            Not entirely clear what it may return. Probably:
            SUCCESS, ERROR, NO_HANDLE. Maybe DEVICE_NOT_FOUND.
        """
        return self.__dll_open_by_unique_name(self.__handle,
                                              self.__to_bytes(unique_name))


    __dll_snap_image = __dll.IC_SnapImage
    __dll_snap_image.restype = C.c_int
    __dll_snap_image.argtypes = GrabberHandlePtr, C.c_int

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
        ret_val : int
            SUCCESS if an image has been snapped, NOT_IN_LIVEMODE if
            the device is not in live mode, i.e., start_live() was
            not called, ERROR if something else went wrong.
        """
        return self.__dll_snap_image(self.__handle, timeout)

    __dll_start_live = __dll.IC_StartLive
    __dll_start_live.restype = C.c_int
    __dll_start_live.argtypes = GrabberHandlePtr, C.c_int

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
        return self.__dll_start_live(self.__handle, bool(show_video))

    __dll_stop_live = __dll.IC_StopLive
    __dll_stop_live.restype = C.c_int
    __dll_stop_live.argtypes = (GrabberHandlePtr,)

    def stop_live(self):
        """Stop live mode.

        Returns
        -------
        ret_val : int
            Unclear what it returns!
        """
        return self.__dll_stop_live(self.__handle)

    __dll_show_device_selection_dialog = __dll.IC_ShowDeviceSelectionDialog
    __dll_show_device_selection_dialog.restype = GrabberHandlePtr
    __dll_show_device_selection_dialog.argtypes = (GrabberHandlePtr,)

    # TODO: implemented just for testing. will be later removed.
    def temp_selection_dialog(self):
        """Show a device-selection dialog.

        This dialog allows to select the video capture device, the
        video norm, video format, input channel and frame rate.

        Returns
        -------
        None.
        """
        self.__handle = __dll_show_device_selection_dialog(self.__handle)

    # TODO: restype did not match signature. Also, the old python
    # implementation assigned the return to __handle. Should be
    # checked.
    __dll_show_property_dialog = __dll.IC_ShowPropertyDialog
    __dll_show_property_dialog.restype = C.c_int  # was GrabberHandlePtr
    __dll_show_property_dialog.argtypes = (GrabberHandlePtr,)

    def temp_property_dialog(self):
        """Show the VCDProperty dialog.

        Requires a video-capture device to be open.

        Returns
        -------
        ret_val : int
            SUCCESS on success, NO_HANDLE if grabber_handle is invalid,
            NO_DEVICE if no device was open, ERROR otherwise.
        """
        return self.__dll_show_property_dialog(self.__handle)

    __dll_unique_name_from_index = __dll.IC_GetUniqueNamefromList
    __dll_unique_name_from_index.restype = C.c_char_p
    __dll_unique_name_from_index.argtypes = (C.c_int,)

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
        device_name : bytes or None
            The unique name of the device (UTF-8). None on failure.
        """
        return self.__dll_unique_name_from_index(index)

    __dll_video_format_by_index = __dll.IC_GetVideoFormat
    __dll_video_format_by_index.restype = C.c_char_p
    __dll_video_format_by_index..argtypes = GrabberHandlePtr, C.c_int

    def __video_format_by_index(self, index):
        """Return the video format specified by index as bytes.

        Parameters
        ----------
        index : int
            Index of the video format to be used. Must be between 0 and
            n_video_formats().

        Returns
        -------
        video_format : bytes or None
            The name of the specified video format as a UTF-8-encoded
            string. Returns None if an error occurred.
        """
        return self.__dll_video_format_by_index(self.__handle, index)

    __dll_n_video_formats = __dll.IC_GetVideoFormatCount
    __dll_n_video_formats.restype = C.c_int
    __dll_n_video_formats.argtypes = (GrabberHandlePtr,)

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

    __dll_video_format_width = __dll.IC_GetVideoFormatWidth
    __dll_video_format_width.restype = C.c_int
    __dll_video_format_width.argtypes = (GrabberHandlePtr,)

    def __video_format_width(self):
        """Return the width of the current video format.

        Likely requires a device to be open. Probably returns
        errors similar to __n_video_formats().

        Returns
        -------
        width : int
            The width of the video format in pixels.
        """
        return self.__dll_video_format_width(self.__handle)

    __dll_video_format_height = __dll.IC_GetVideoFormatHeight
    __dll_video_format_height.restype = C.c_int
    __dll_video_format_height.argtypes = (GrabberHandlePtr,)

    def __video_format_height(self):
        """Return the height of the current video format.

        Likely requires a device to be open. Probably returns
        errors similar to __n_video_formats().

        Returns
        -------
        height : int
            The height of the video format in pixels."""
        return self.__dll_video_format_height(self.__handle)

    __dll_get_image_ptr = __dll.IC_GetImagePtr
    __dll_get_image_ptr.restype = C.c_void_p  # void to also allow None
    __dll_get_image_ptr.argtypes = (GrabberHandlePtr,)

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

    # TODO: unclear which format the file name should have.
    #       would str(Path.resolve()) work? returns separators as '\\'
    __dll_load_settings_from_file = __dll.IC_LoadDeviceStateFromFile
    __dll_load_settings_from_file.restype = GrabberHandlePtr
    __dll_load_settings_from_file.argtypes = GrabberHandlePtr, C.c_char_p

    def load_settings_from_file(self, filename):
        """Load a device-settings file.

        On success the device is opened automatically.

        Parameters
        ----------
        file_name : bytes
            UTF-8-encoded path to the settings file to be loaded.

        Returns
        -------
        None.
        """
        self.__handle = self.__dll_load_settings_from_file(
            self.__handle, self.__to_bytes(filename)
            )

    # TODO: similar question as load_settings_from_file
    __dll_save_settings_to_file = __dll.IC_SaveDeviceStateToFile
    __dll_save_settings_to_file.restype = C.c_int
    __dll_save_settings_to_file.argtypes = GrabberHandlePtr, C.c_char_p

    def save_settings(self, filename):
        """Save the state of a video capture device to file.

        Parameters
        ----------
        file_name : bytes
            UTF-8-encoded path to the settings file where
            settings should be saved.

        Returns
        -------
        ret_val : int
            SUCCESS if successful, ERROR otherwise.
        """
        return self.__dll_save_settings_to_file(self.__handle,
                                                self.__to_bytes(filename))

    # TODO: unclear what it returns exactly.
    __dll_get_camera_property = __dll.IC_GetCameraProperty
    __dll_get_camera_property.restype = C.c_int
    __dll_get_camera_property.argtypes = (GrabberHandlePtr, C.c_int,
                                          C.POINTER(C.c_long))

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

        C-Returns
        ---------
        ret_val : int
            SUCCESS or ERROR. Probably also NOT_AVAILABLE, and perhaps
            PROPERTY_ITEM_NOT_AVAILABLE, PROPERTY_ELEMENT_NOT_AVAILABLE,
            PROPERTY_ELEMENT_WRONG_INTERFACE
        """
        cam_property = CameraProperty.from_value_or_name(cam_property)
        prop = C.c_long()
        ret_val = self.__dll_get_camera_property(self.__handle,
                                                 cam_property.value,
                                                 prop)
        return prop.value

    __dll_set_camera_property = __dll.IC_SetCameraProperty
    __dll_set_camera_property.restype = C.c_int
    __dll_set_camera_property.argtypes = (GrabberHandlePtr, C.c_int,
                                          C.c_long)

    def set_camera_property(self, cam_property, value):
        """Set a camera property.

        Parameters
        ----------
        property : str, int, or CameraProperty
            The property to be set
        value : int
            The value of the property. The value should be in the range
            of the specified property.

        Returns
        -------
        ret_val : int
            SUCCESS or ERROR. Returns error if the property could not
            be set, e.g., if the value is out of range or if it is
            currently set to 'auto'.
        """
        cam_property = CameraProperty.from_value_or_name(cam_property)
        return self.__dll_get_camera_property(self.__handle,
                                              cam_property.value,
                                              value)

    # TODO: check return values
    __dll_get_camera_property_range = __dll.IC_CameraPropertyGetRange
    __dll_get_camera_property_range.restype = C.c_int
    __dll_get_camera_property_range.argtypes = (GrabberHandlePtr, C.c_int,
                                                C.POINTER(C.c_long),
                                                C.POINTER(C.c_long))

    def __get_camera_property_range(self, cam_property):
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

        C-Returns
        ---------
        ret_val : int
            SUCCESS or ERROR. Perhaps others?
        """
        cam_property = CameraProperty.from_value_or_name(cam_property)
        prop_min = C.c_long()
        prop_max = C.c_long()

        ret_val = self.__dll_get_camera_property(
            self.__handle, cam_property.value, prop_min, prop_max
            )
        return prop_min, prop_max

    # TODO: this has to be tested. Would be nice to have it working,
    #       at least as an alternative to VCDPropertyInspector.
    # Judging from .h EnumCallback takes a name (of an enum or property?)
    # and some data, and returns int. Unclear id the int returned is the
    # value of the property enum or what. It is unclear also if this can
    # be a python function (probably yes)
    EnumCallback = C.CFUNCTYPE(C.c_int, C.c_char_p, C.c_void_p)

    # TODO: new from dll. Make sure it works.
    enumerate_properties = __dll.IC_enumProperties
    enumerate_properties.restype = C.c_int
    enumerate_properties.argtypes = (GrabberHandlePtr, EnumCallback,
                                     C.c_void_p)
    enumerate_properties.__doc__ = (
        """Enumerate properties of an open device.

        This method needs testing as it was not originally wrapped.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object.
        callback : EnumCallback
            Function that will be called on enums? Has signature
            (char* name, void* data) and returns int
        callback_data : int or None (== c_void_p)
            Data passed to the callback function

        Returns
        -------
        ret_val : int
            SUCCESS, ERROR or other errors."""
        )

    # TODO: it is undocumented which properties are available.
    # Perhaps using VCDPropertyInspector of IC Imaging Control
    # would be a way to check whether enumerate_properties works?
    # A decent read: https://www.theimagingsource.com/support/documentation
    #                        /ic-imaging-control-cpp/tech_VCDProperties.htm
    set_property_value = __dll.IC_SetPropertyValue
    set_property_value.restype = C.c_int
    set_property_value.argtypes = (GrabberHandlePtr, C.c_char_p,
                                   C.c_char_p, C.c_int)
    set_property_value.__doc__ = (
        """Set the value of a property.
        
        By exclusion, it looks like this method is using the
        DirectShow interface, i.e., only setting integer values
        to properties. However, it seems it is only using 'int'
        (16-bit) rather than the 'long' (32-bit) that DirectShow
        is supposed to use.
        
        doc: "The interface IVCDRangeProperty, for example, provides
             access to the property that is equivalent to the property
             access of DirectShow."

        tisgrabber.h suggests: "For a list of properties and elements
        use the VCDPropertyInspector of IC Imaging Control."

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object.
        property_name : bytes
            The name of the property whose value is to be changed
        element_of_property : bytes or None
            The name of the setting to be accessed (examples: b'Value',
            b'Auto'). If None, b'Value' is used.
        value : int
            The value to be set for the property.

        Returns
        -------
            ret_val : int
                SUCCESS, NO_HANDLE (invalid grabber handle), NO_DEVICE
                (no device open), PROPERTY_ITEM_NOT_AVAILABLE (the
                property item is not available),
                PROPERTY_ELEMENT_NOT_AVAILABLE (the element of the
                property is not available),
                PROPERTY_ELEMENT_WRONG_INTERFACE (the property element
                exists but it does not have the requested attribute,
                e.g., b'Exposure' b'Auto' has no 'range')."""
        )

    get_property_value = __dll.IC_GetPropertyValue
    get_property_value.restype = C.c_int
    get_property_value.argtypes = (GrabberHandlePtr, C.c_char_p,
                                   C.c_char_p, C.POINTER(C.c_long))
    get_property_value.__doc__ = (
        """Return the current value of a property.

        tisgrabber.h suggests: "For a list of properties and elements
        use the VCDPropertyInspector of IC Imaging Control."
        """


# ############################################################################
    SetPropertySwitch = __dll.IC_SetPropertySwitch
    SetPropertySwitch.restype = C.c_int
    SetPropertySwitch.argtypes= (GrabberHandlePtr, C.c_char_p,
                                 C.c_char_p, C.c_int)

    GetPropertySwitch = __dll.IC_GetPropertySwitch
    GetPropertySwitch.restype = C.c_int
    GetPropertySwitch.argtypes= (GrabberHandlePtr, C.c_char_p,
                                 C.c_char_p, C.POINTER(C.c_long))

# ############################################################################

    IsPropertyAvailable =  __dll.IC_IsPropertyAvailable
    IsPropertyAvailable.restype =  C.c_int
    IsPropertyAvailable.argtypes=  (GrabberHandlePtr, C.c_char_p,
                                    C.c_char_p)

    IsCameraPropertyAvailable = __dll.IC_IsCameraPropertyAvailable
    IsCameraPropertyAvailable.restype = C.c_int
    IsCameraPropertyAvailable.argtypes= (GrabberHandlePtr, C.c_char_p)

    PropertyOnePush = __dll.IC_PropertyOnePush
    PropertyOnePush.restype = C.c_int
    PropertyOnePush.argtypes = (GrabberHandlePtr, C.c_char_p,
                                C.c_char_p)


    SetPropertyAbsoluteValue = __dll.IC_SetPropertyAbsoluteValue
    SetPropertyAbsoluteValue.restype = C.c_int
    SetPropertyAbsoluteValue.argtypes = (GrabberHandlePtr, C.c_char_p,
                                         C.c_char_p, C.c_float)

    GetPropertyAbsoluteValue = __dll.IC_GetPropertyAbsoluteValue
    GetPropertyAbsoluteValue.restype = C.c_int
    GetPropertyAbsoluteValue.argtypes = (GrabberHandlePtr, C.c_char_p,
                                         C.c_char_p, C.POINTER(C.c_float))

    GetPropertyValueRange = __dll.IC_GetPropertyValueRange
    GetPropertyValueRange.restype = C.c_int
    GetPropertyValueRange.argtypes= (GrabberHandlePtr, C.c_char_p,
                                     C.c_char_p, C.POINTER(C.c_long),
                                     C.POINTER(C.c_long))



# ############################################################################
    EnableAutoCameraProperty = __dll.IC_EnableAutoCameraProperty
    EnableAutoCameraProperty.restype = C.c_int
    EnableAutoCameraProperty.argtypes = (GrabberHandlePtr, C.c_int,
                                         C.c_int)
    EnableVideoAutoProperty = __dll.IC_EnableAutoVideoProperty
    EnableVideoAutoProperty.restype = C.c_int
    EnableVideoAutoProperty.argtypes = (GrabberHandlePtr, C.c_int,
                                        C.c_int)


# ############################################################################

    # definition of the frameready callback
    FRAME_READY_CALLBACK = C.CFUNCTYPE(C.c_void_p, C.c_int,
                                     C.POINTER(C.c_ubyte),
                                     C.c_ulong, C.py_object)

    # set callback function
    SetFrameReadyCallback = __dll.IC_SetFrameReadyCallback
    SetFrameReadyCallback.restype = C.c_int
    SetFrameReadyCallback.argtypes = (GrabberHandlePtr, FRAME_READY_CALLBACK,
                                      C.py_object)

    SetContinuousMode = __dll.IC_SetContinuousMode

    SaveImage = __dll.IC_SaveImage
    SaveImage.restype = C.c_int
    SaveImage.argtypes = (C.c_void_p, C.c_char_p, C.c_int, C.c_int)

    OpenVideoCaptureDevice = __dll.IC_OpenVideoCaptureDevice
    OpenVideoCaptureDevice.restype = C.c_int
    OpenVideoCaptureDevice.argtypes = C.c_void_p, C.c_char_p

    CloseVideoCaptureDevice = __dll.IC_CloseVideoCaptureDevice
    CloseVideoCaptureDevice.restype = C.c_int
    CloseVideoCaptureDevice.argtypes = C.c_void_p, C.c_char_p

# ############################################################################

    IsTriggerAvailable = __dll.IC_IsTriggerAvailable
    IsTriggerAvailable.restype = C.c_int
    IsTriggerAvailable.argtypes = (GrabberHandlePtr,)

    EnableTrigger = __dll.IC_EnableTrigger
    EnableTrigger.restype = C.c_int
    EnableTrigger.argtypes = GrabberHandlePtr, C.c_int

    SoftwareTrigger = __dll.IC_SoftwareTrigger
    SoftwareTrigger.restype = C.c_int
    SoftwareTrigger.argtypes = (GrabberHandlePtr,)

# ############################################################################

    ResetProperties = __dll.IC_ResetProperties
    ResetProperties.restype = C.c_int
    ResetProperties.argtypes = (GrabberHandlePtr,)

    GetTriggerModes = __dll.IC_GetTriggerModes
    GetTriggerModes.restype = C.c_int
    GetTriggerModes.argtypes = (GrabberHandlePtr, C.c_char_p,
                                C.c_int)

#    int AC IC_SoftwareTrigger(HGRABBER hGrabber);

#int AC IC_EnableTrigger( HGRABBER hGrabber, int iEnable );
#int AC IC_IsTriggerAvailable( HGRABBER hGrabber );
#this is above functions

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

    def s(self,strin):
        if sys.version[0] == "2":
            return strin
        if type(strin) == "byte":
            return strin
        return strin.encode("utf-8")

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

    def is_device_valid(self):
        return GrabberDLL.is_device_valid(self._handle)

    def set_window_handle(self, Hwnd):
        return GrabberDLL.set_window_handle(self._handle, Hwnd)

    def save_device_state_to_file(self, FileName):
        return GrabberDLL.save_device_state_to_file(self._handle, self.s(FileName))

    def load_device_state_from_file(self,FileName):
        self._handle = GrabberDLL.load_device_state_from_file(self._handle,self.s(FileName))

    def get_input_channels(self):
        n_channels = GrabberDLL.get_input_channel_count(self._handle)
        input_channels = [GrabberDLL.get_input_channel(self._handle, i)
                          for i in range(n_channels)]
        return input_channels

    def get_video_norms(self):
        n_norms = GrabberDLL.get_video_norm_count(self._handle)
        video_norms = [GrabberDLL.get_video_norm(self._handle, i)
                       for i in range(n_norms)]
        return video_norms

    def get_image_description(self):
        lWidth=C.c_long()
        lHeight= C.c_long()
        iBitsPerPixel=C.c_int()
        COLORFORMAT=C.c_int()

        Error = GrabberDLL.get_image_description(self._handle, lWidth,
                                    lHeight,iBitsPerPixel,COLORFORMAT)
        return (lWidth.value,lHeight.value,iBitsPerPixel.value,COLORFORMAT.value)

    def get_image_ptr(self):
        ImagePtr = GrabberDLL.get_image_ptr(self._handle)

        return ImagePtr

    def GetImage(self):
        BildDaten = self.get_image_description()[:4]
        lWidth = BildDaten[0]
        lHeight = BildDaten[1]
        iBitsPerPixel = BildDaten[2] // 8 # // 8 was commented out??

        buffer_size = lWidth*lHeight*iBitsPerPixel*C.sizeof(C.c_uint8)
        img_ptr = self.get_image_ptr()

        Bild = C.cast(img_ptr, C.POINTER(C.c_ubyte * buffer_size))


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
        buffer_size = lWidth*lHeight*iBytesPerPixel*C.sizeof(C.c_uint8)
        img_ptr = self.get_image_ptr()

        Bild = C.cast(img_ptr, C.POINTER(C.c_ubyte * buffer_size))
        pixeltype = np.uint8
        if  BildDaten[3] == 4: #SinkFormats.Y16:
            pixeltype = np.uint16
            iBytesPerPixel = 1

        # Look at numpy.ctypeslib.as_array(obj, shape=None)
        img = np.ndarray(buffer = Bild.contents,
                     dtype = pixeltype,
                     shape = (lHeight,
                              lWidth,
                              iBytesPerPixel))
        return img


    def get_camera_property(self,iProperty):
        lFocusPos = C.c_long()
        Error = GrabberDLL.get_camera_property(self._handle,iProperty, lFocusPos)
        return (lFocusPos.value)

    def set_camera_property(self,iProperty,iValue):
        Error = GrabberDLL.set_camera_property(self._handle,iProperty, iValue)
        return (Error)

    def set_property_value(self, Property, Element, Value ):
        error = GrabberDLL.set_property_value(self._handle,
                                                self.s(Property),
                                                self.s(Element),
                                                Value)
        return error


    def get_property_value(self, Property, Element ):
        Value = C.c_long()
        error = GrabberDLL.get_property_value(self._handle,
                                                self.s(Property),
                                                self.s(Element),
                                                Value)
        return Value.value

    def IsCameraPropertyAvailable(self, Property):

        error = GrabberDLL.IsCameraPropertyAvailable(self._handle,
                                                   self.s(Property))
        return error

    def PropertyAvailable(self, Property):
        Null = None
        error = GrabberDLL.IsPropertyAvailable(self._handle,
                                                        self.s(Property),
                                                        Null)
        return error


    def SetPropertySwitch(self, Property, Element, Value):
        error = GrabberDLL.SetPropertySwitch(self._handle,
                                                  self.s(Property),
                                                  self.s(Element),
                                                  Value)
        return error

    def GetPropertySwitch(self, Property, Element, Value):
        lValue = C.c_long()
        error = GrabberDLL.GetPropertySwitch(self._handle,
                                                  self.s(Property),
                                                  self.s(Element),
                                                  lValue)
        Value[0] = lValue.value
        return error

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
        lValue = C.c_float()
        error = GrabberDLL.GetPropertyAbsoluteValue(self._handle,
                                                self.s(Property),
                                                self.s(Element),
                                                lValue)
        Value[0] = lValue.value
        return error

    def GetPropertyValueRange(self, Property, Element, Value_min, Value_max):
        minValue = C.c_long()
        maxValue = C.c_long()

        error = GrabberDLL.GetPropertyValueRange(self._handle,
                                                     self.s(Property),
                                                     self.s(Element),
                                                     minValue,
                                                     maxValue)
        Value_min[0] = minValue.value
        Value_max[0] = maxValue.value
        return error

    def SaveImage(self,FileName, FileType, Quality=75):
        ''' Saves the last snapped image. Can by of type BMP or JPEG.
        :param FileName : Name of the mage file
        :param FileType : Determines file type, can be "JPEG" or "BMP"
        :param Quality : If file typ is JPEG, the qualitly can be given from 1 to 100.
        :return: Error code
        '''
        #return GrabberDLL.SaveImage(self._handle, self.s(FileName), IC.ImageFileTypes[self.s(FileType)],Quality)

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