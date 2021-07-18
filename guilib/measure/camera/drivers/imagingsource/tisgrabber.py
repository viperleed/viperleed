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

dll_path = Path(__file__).parent.resolve
print(f"{dll_path=}")
try:
    os.add_dll_directory(dll_path)
except AttributeError:
    pass
sys.path.append(dll_path)


class SinkFormats(Enum):
    Y800 = 0
    RGB24 = 1
    RGB32 = 2
    UYVY = 3
    Y16 = 4
    NONE = 5  # used as return value
    MEGA = 65536  # Borland C++ 6 compatibility


class CameraProperty(Enum):
    PAN = 0
	TILT = 1
	ROLL = 2
	ZOOM = 3
	EXPOSURE = 4
	IRIS = 5
	FOCUS = 6
	MEGA = 65536  # Borland C++ 6 compatibility


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

    if sys.maxsize > 2**32 :
        __dll = C.windll.LoadLibrary("tisgrabber_x64.dll")
    else:
        __dll = C.windll.LoadLibrary("tisgrabber.dll")

    def __init__(self, **keyargs):
        """Initialize the Albatross from the keyword arguments."""
        self.__dict__.update(keyargs)

    # ctypes.POINTER returns a class
    GrabberHandlePtr = C.POINTER(GrabberHandle)

    init_library = __dll.IC_InitLibrary
    init_library.__doc__ = (
        """Initialize the ICImagingControl class library.

        This function must be called only once before any other
        functions of this library are called.

        Parameters
        ----------
        license_key : str or None
            The license key. Pass None if using a 'trial' version.

        Returns
        -------
        ret_val : {DLLReturns.SUCCESS.value, DLLReturns.ERROR.value}
            ERROR if license key is wrong, or for other generic errors"""
        )
    # Call init_library with None straight away. It was done the
    # same way in the original tisgrabber.py
    init_library(None)

    get_device_count =  __dll.IC_GetDeviceCount
    get_device_count.restype = C.c_int
    get_device_count.argtypes = None
    get_device_count.__doc__ = (
        """Get the number of the currently available devices.

        This function creates an internal array of all connected video-
        capture devices. With each call to this function, this array is
        rebuilt. The name and the unique name can be retrieved from the
        internal array using the function IC_GetDevice() and the method
        get_unique_name_from_list(). Both are useful to for retrieve
        device names for opening devices.

        Returns
        -------
        n_devices : int
            The number of devices found. DLLReturns.NO_HANDLE.value is
            returned in case of error."""
        )

    get_unique_name_from_list = __dll.IC_GetUniqueNamefromList
    get_unique_name_from_list.restype = C.c_char_p
    get_unique_name_from_list.argtypes = (C.c_int,)
    get_unique_name_from_list.__doc__ = (
        """Get unique name of a device given its index.

        The unique device name consist of the device name and its
        serial number. It allows to differentiate between more devices
        of the same type connected to the computer. The unique device
        name is passed to the method open_device_by_unique_name().

        Parameters
        ----------
        index : int
            The index of the device whose name is to be returned. It
            must be between 0 and get_device_count().

        Returns
        -------
        device_name : str or None
            The unique name of the device. Returns None on failure."""
        )

    create_grabber = __dll.IC_CreateGrabber
    create_grabber.restype = GrabberHandlePtr
    create_grabber.argtypes = None
    create_grabber.__doc__ = (
        """Return a new grabber handle.

        Unused grabber handles should be release by calling
        IC_ReleaseGrabber. There is no need to release a handle
        when closing a device, as handles can be reused.

        Returns
        -------
        handle : GrabberHandlePtr
            New handle. Can be used to open devices."""
        )

    open_device_by_unique_name = __dll.IC_OpenDevByUniqueName
    open_device_by_unique_name.restype = C.c_int
    open_device_by_unique_name.argtypes = GrabberHandlePtr, C.c_char_p
    open_device_by_unique_name.__doc__ = (
        """Open a video capture by using its unique name.

        Use IC_GetUniqueName() to retrieve the unique name of an
        open device. get_unique_name_from_list(index) can be used
        before the device is open.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        unique_name : bytes
            Unique name of the device to be opened

        Returns
        -------
        ret_val : int
            Not entirely clear what it may return. Probably:
            SUCCESS, ERROR, NO_HANDLE. Maybe DEVICE_NOT_FOUND."""
        )

    set_video_format = __dll.IC_SetVideoFormat
    set_video_format.restype = C.c_int
    set_video_format.argtypes = GrabberHandlePtr, C.c_char_p
    set_video_format.__doc__ = (
        """Set a video format for the current video-capture device.

        The video format must be supported by the current device.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        video_format : bytes
            A UTF-8-encoded representation of the desired video format.

        Returns
        -------
        ret_val : {DLLReturns.SUCCESS.value, DLLReturns.ERROR.value}
            ERROR if something went wrong."""
        )

    set_frame_rate = __dll.IC_SetFrameRate
    set_frame_rate.restype = C.c_int
    set_frame_rate.argtypes = GrabberHandlePtr, C.c_float
    set_frame_rate.__doc__ = (
        """Set a new frame rate.

        A new frame rate can not be set while in live mode. Call
        stop_live(), then try again.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        frame_rate : float
            The new frame rate.

        Returns
        -------
        ret_val : DLLReturns.value
            SUCCESS, NOT_AVAILABLE (not supported by the device),
            NO_HANDLE (invalid grabber handle), NO_DEVICE (no device
            opened), NOT_IN_LIVEMODE (device is in live mode)."""
        )

    get_video_format_width = __dll.IC_GetVideoFormatWidth
    get_video_format_width.restype = C.c_int
    get_video_format_width.argtypes = (GrabberHandlePtr,)
    get_video_format_width.__doc__ = (
        """Return the width of the video format.

        Likely requires a device to be open. Probably returns
        errors similar to get_video_format_count().

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        width : int
            The width of the video format in pixels."""
        )

    get_video_format_height = __dll.IC_GetVideoFormatHeight
    get_video_format_height.restype = C.c_int
    get_video_format_height.argtypes = (GrabberHandlePtr,)
    get_video_format_height.__doc__ = (
        """Return the height of the video format.

        Likely requires a device to be open. Probably returns
        errors similar to get_video_format_count().

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        height : int
            The height of the video format in pixels."""
        )

    get_video_format_count = __dll.IC_GetVideoFormatCount
    get_video_format_count.restype = C.c_int
    get_video_format_count.argtypes = (GrabberHandlePtr,)
    get_video_format_count.__doc__ = (
        """Return the number of video formats supported.

        Requires a video-capture device to be open.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        n_formats : int
            The number of video formats supported by the open device
            if successful. Otherwise NO_DEVICE (no device open) or
            NO_HANDLE (invalid grabber_handle)."""
        )

    get_video_format = __dll.IC_GetVideoFormat
    get_video_format.restype = C.c_char_p
    get_video_format.argtypes = GrabberHandlePtr, C.c_int
    get_video_format.__doc__ = (
        """Return the video format specified by index as a string.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        index : int
            Index of the video format to be used. Must be between 0 and
            get_video_format_count().  Call get_video_format_count()
            before this method, otherwise it will always fail.

        Returns
        -------
        video_format : bytes or None
            The name of the specified video format. Returns None
            if an error occurred."""
        )

    # TODO: unclear what this returns in case of error
    get_sink_format = __dll.IC_GetFormat
    get_sink_format.restype = C.c_int
    get_sink_format.argtypes = (GrabberHandlePtr,)
    get_sink_format.__doc__ = (
        """Return the color format of the sink type currently set.

        The sink type describes the format of the buffer where snapped
        images are stored. The method returns a valid value only after
        IC_PreprareLive() or start_live() are called. Notice that the
        sink type may differ from the video format.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        color_format : int
            One of the values in SinkFormats. Returns !!UNCLEAR WHAT!!
            if not in live mode or if an error occurred."""
        )

    set_sink_format = __dll.IC_SetFormat
    set_sink_format.restype = C.c_int
    set_sink_format.argtypes = GrabberHandlePtr, C.c_int
    set_sink_format.__doc__ = (
        """Set the sink type.

        This method must be called before images can be snapped.
        The sink type describes the format of the buffer where
        snapped images are stored. Notice that the sink type may
        differ from the video format.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        sink_type : int
            One of the values in SinkFormats. Note that UYVY can
            only be used in conjunction with a UYVY video format

        Returns
        -------
        ret_val : int
            SUCCESS on success, ERROR if something went wrong."""
        )

    # TODO: it is unclear what an input channel is! ANALOG INPUT. Probably not useful.
    get_input_channel_count = __dll.IC_GetInputChannelCount
    get_input_channel_count.restype = C.c_int
    get_input_channel_count.argtypes = (GrabberHandlePtr,)
    get_input_channel_count.__doc__ = (
        """Return the number of available input channels for a device.

        Requires a video-capture device to be open.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        n_channels : int
            The number of channels on success, otherwise NO_DEVICE
            (no video capture device open) or NO_HANDLE (invalid
            grabber_handle)."""
        )

    # TODO: probably not useful
    get_input_channel = __dll.IC_GetInputChannel
    get_input_channel.restype = C.c_char_p
    get_input_channel.argtypes = GrabberHandlePtr, C.c_int
    get_input_channel.__doc__ = (
        """Return the input channel at index as a string.

        A device should be open for this method to succeed. This method
        will fail unless get_input_channel_count() was called before.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        index : int
            Index of the channel to be returned. index must be between
            0 and get_input_channel_count().

        Returns
        -------
        channel : bytes or None
            The name of the specified input channel (UTF-8). Returns
            None if an error occurred."""
        )

    # TODO: video norm == Analog/PAL/NTSC/SECAM: https://docs.microsoft.com/en-us/windows/win32/api/strmif/ne-strmif-analogvideostandard)
    get_video_norm_count = __dll.IC_GetVideoNormCount
    get_video_norm_count.restype = C.c_int
    get_video_norm_count.argtypes = (GrabberHandlePtr,)
    get_video_norm_count.__doc__ = (
        """Return the number of the video norms for the current device.

        A device should be open for this method to succeed.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        n_norms : int
            The number of video norms on success, otherwise NO_DEVICE
            (no video capture device open) or NO_HANDLE (invalid
            grabber_handle)."""
        )

    get_video_norm = __dll.IC_GetVideoNorm
    get_video_norm.restype = C.c_char_p
    get_video_norm.argtypes = GrabberHandlePtr, C.c_int
    get_video_norm.__doc__ = (
        """Return the video norm at index as a string.

        This method will fail unless get_video_norm_count() was
        called before.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        index : int
            Index of the video norm to be returned. index must be
            between 0 and get_video_norm_count().

        Returns
        -------
        video_norm : bytes or None
            The name of the specified video norm (UTF-8). Returns
            None if an error occurred."""
        )

    start_live = __dll.IC_StartLive
    start_live.restype = C.c_int
    start_live.argtypes = GrabberHandlePtr, C.c_int
    start_live.__doc__ = (
        """Start camera in live mode.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        show_video : {0, 1}
            If 0, the video is not shown, but frames are delivered
            (i.e., callbacks can be used). If 1, a video is also shown

        Returns
        -------
        ret_val : int
            SUCCESS on success, ERROR if something went wrong."""
        )

    stop_live = __dll.IC_StopLive
    stop_live.restype = C.c_int
    stop_live.argtypes = (GrabberHandlePtr,)
    stop_live.__doc__ = (
        """Stop live mode.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        ret_val : int
            Unclear what it returns!
        """
        )

    set_window_handle = __dll.IC_SetHWnd
    set_window_handle.restype = C.c_int
    set_window_handle.argtypes = GrabberHandlePtr, C.c_int
    set_window_handle.__doc__ = (
        """Assign a window handle to display the video in.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        window_handle : int
            The handle of the window where a live video will be shown.

        Returns
        -------
        ret_val : int
            SUCCESS if an image has been snapped, ERROR if something
            went wrong."""
        )

    snap_image = __dll.IC_SnapImage
    snap_image.restype = C.c_int
    snap_image.argtypes = GrabberHandlePtr, C.c_int
    snap_image.__doc__ = (
        """Snap an image.

        The video capture device must be set to live mode and a sink
        type has to be set before this call. The format of the snapped
        images depends on the selected sink type.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        timeout : int
            Time interval in milliseconds after which the device will
            time out. A value of -1 corresponds to no time-out.

        Returns
        -------
        ret_val : int
            SUCCESS if an image has been snapped, NOT_IN_LIVEMODE if
            the device is not in live mode, i.e., start_live() was
            not called, ERROR if something else went wrong."""
        )

    # TODO: unclear what 'color format' means. Is it == video_format??
    #       is it == sink_type ??
    get_image_description = __dll.IC_GetImageDescription
    get_image_description.restype = C.c_int
    get_image_description.argtypes = (GrabberHandlePtr,
                                      C.POINTER(C.c_long),
                                      C.POINTER(C.c_long),
                                      C.POINTER(C.c_int),
                                      C.POINTER(C.c_int))
    get_image_description.__doc__ = (
        """Fetch properties of the current video format and sink type.

        Note that the pointers below can also be passed as pure
        c_int/c_long, as ctypes cares of passing byref().

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object
        width, height : ctypes.pointer(ctypes.c_long) or ctypes.c_long
            Width and height of the image in pixels will be
            stored in these 1-element buffers
        bits_per_pixel : ctypes.pointer(ctypes.c_int) or ctypes.c_int
            The number of bits for each pixel, accounting for color
            channels, will be stored in this 1-element buffer
        format : ctypes.pointer(ctypes.c_int) or ctypes.c_int
            The current color format value (??) is stored in this
            1-element buffer.

        Returns
        -------
        ret_val : int
            SUCCESS on success, ERROR if something went wrong."""
        )

    get_image_ptr = __dll.IC_GetImagePtr
    get_image_ptr.restype = C.c_void_p  # void to also allow None
    get_image_ptr.argtypes = (GrabberHandlePtr,)
    get_image_ptr.__doc__ = (
        """Get a pointer to the pixel data of the last snapped image.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        pointer : ctypes.c_void_p
            A pointer to the the first byte of the bottommost line
            of the image (images are saved from bottom to top!).
            Returns None in case of error."""
        )

    show_device_selection_dialog = __dll.IC_ShowDeviceSelectionDialog
    show_device_selection_dialog.restype = GrabberHandlePtr
    show_device_selection_dialog.argtypes = (GrabberHandlePtr,)
    show_device_selection_dialog.__doc__ = (
        """Show the device selection dialog.

        This dialogs enables to select the video capture device, the
        video norm, video format, input channel and frame rate.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        grabber_handle : GrabberHandlePtr
            The same handle passed if valid, a new one otherwise."""
        )

    # TODO: restype did not match signature
    show_property_dialog = __dll.IC_ShowPropertyDialog
    show_property_dialog.restype = C.c_int  # was GrabberHandlePtr
    show_property_dialog.argtypes = (GrabberHandlePtr,)
    show_property_dialog.__doc__ = (
        """Show the VCDProperty dialog.

        Requires a video-capture device to be open.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        ret_val : int
            SUCCESS on success, NO_HANDLE if grabber_handle is invalid,
            NO_DEVICE if no device was open, ERROR otherwise."""
        )

    is_device_valid = __dll.IC_IsDevValid
    is_device_valid.restype = C.c_int
    is_device_valid.argtypes = (GrabberHandlePtr,)
    is_device_valid.__doc__ = (
        """Return whether there is a valid open device.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object

        Returns
        -------
        ret_val : int
            SUCCESS if a valid device is open, ERROR otherwise."""
        )

    load_device_state_from_file = __dll.IC_LoadDeviceStateFromFile
    load_device_state_from_file.restype = GrabberHandlePtr
    load_device_state_from_file.argtypes = GrabberHandlePtr, C.c_char_p
    load_device_state_from_file.__doc__ = (
        """Load a device-settings file.

        On success the device is opened automatically.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr or None
            Handle to grabber object. If None, a new a new handle is
            created and returned. This should be released by a call
            to IC_ReleaseGrabber() when it is no longer needed.
        file_name : bytes
            UTF-8-encoded path to the settings file to be loaded.

        Returns
        -------
        grabber_handle : GrabberHandlePtr
            The same handle if it was valid, a new one if not."""
        )

    save_device_state_to_file = __dll.IC_SaveDeviceStateToFile
    save_device_state_to_file.restype = C.c_int
    save_device_state_to_file.argtypes = GrabberHandlePtr, C.c_char_p
    save_device_state_to_file.__doc__ = (
        """Save the state of a video capture device to file.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object.
        file_name : bytes
            UTF-8-encoded path to the settings file where
            settings should be saved.

        Returns
        -------
        ret_val : int
            SUCCESS if successful, ERROR otherwise."""
        )

    # TODO: unclear what it returns exactly.
    get_camera_property = __dll.IC_GetCameraProperty
    get_camera_property.restype = C.c_int
    get_camera_property.argtypes = (GrabberHandlePtr, C.c_int,
                                    C.POINTER(C.c_long))
    get_camera_property.__doc__ = (
        """Retrieve the value of a camera property.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object.
        property : CameraProperty.value
            The property whose value is to be read
        value : ctypes.pointer(ctypes.c_long) or ctypes.long
            The value of the property will be stored in this buffer

        Returns
        -------
        ret_val : int
            SUCCESS or ERROR. Probably also NOT_AVAILABLE, and perhaps
            PROPERTY_ITEM_NOT_AVAILABLE, PROPERTY_ELEMENT_NOT_AVAILABLE,
            PROPERTY_ELEMENT_WRONG_INTERFACE
        """
        )

    set_camera_property = __dll.IC_SetCameraProperty
    set_camera_property.restype = C.c_int
    set_camera_property.argtypes = (GrabberHandlePtr, C.c_int, C.c_long)
    set_camera_property.__doc__ = (
        """Set a camera property.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object.
        property : CameraProperty.value
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
        )

    # TODO: call in TIS_CAM; check return values
    get_camera_property_range = __dll.IC_CameraPropertyGetRange
    get_camera_property_range.restype = C.c_int
    get_camera_property_range.argtypes = (GrabberHandlePtr, C.c_int,
                                          C.POINTER(C.c_long),
                                          C.POINTER(C.c_long))
    get_camera_property_range.__doc__ = (
        """Retrieve the minimum and maximum value of a property.

        Parameters
        ----------
        grabber_handle : GrabberHandlePtr
            Handle to grabber object.
        property : CameraProperty.value
            The property whose range is to be read
        min, max : ctypes.pointer(ctypes.c_long) or c_long
            The minimum and maximum value of the property will be
            stored in these 1-element buffers.

        Returns
        -------
        ret_val : int
            SUCCESS or ERROR. Perhaps others?"""
        )
    
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

    GetFrameRate =  __dll.IC_GetFrameRate
    GetFrameRate.restype = C.c_float
    GetFrameRate.argtypes = (GrabberHandlePtr,)

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

class TIS_CAM(object):
    @property
    def callback_registered(self):
        return self._callback_registered

    def __init__(self):

        self._handle = C.POINTER(GrabberHandle)
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

    def open(self,unique_device_name):
        """ Open a device

        unique_device_name : The name and serial number of the device to be opened. The device name and serial number are separated by a space.
        """
        test = GrabberDLL.open_device_by_unique_name(self._handle,
                                                   self.s(unique_device_name))

        return test

    def show_device_selection_dialog(self):
        self._handle = GrabberDLL.show_device_selection_dialog(self._handle)

    def show_property_dialog(self):
        self._handle = GrabberDLL.show_property_dialog(self._handle)

    def is_device_valid(self):
        return GrabberDLL.is_device_valid(self._handle)

    def set_window_handle(self, Hwnd):
        return GrabberDLL.set_window_handle(self._handle, Hwnd)

    def save_device_state_to_file(self, FileName):
        return GrabberDLL.save_device_state_to_file(self._handle, self.s(FileName))

    def load_device_state_from_file(self,FileName):
        self._handle = GrabberDLL.load_device_state_from_file(self._handle,self.s(FileName))


    def SetVideoFormat(self,Format):
        return GrabberDLL.set_video_format(self._handle, self.s(Format))

    def SetFrameRate(self,FPS):
        return GrabberDLL.set_frame_rate(self._handle, FPS)

    def get_video_format_width(self):
        return GrabberDLL.get_video_format_width(self._handle)

    def get_video_format_height(self):
        return GrabberDLL.get_video_format_height(self._handle)


    def GetDevices(self):
        self._Devices=[]
        iDevices = GrabberDLL.get_device_count()
        for i in range(iDevices):
            self._Devices.append(GrabberDLL.get_unique_name_from_list(i))
        return self._Devices


    def GetVideoFormats(self):
        self._Properties=[]
        iVideoFormats = GrabberDLL.get_video_format_count(self._handle)
        for i in range(iVideoFormats):
            self._Properties.append(GrabberDLL.get_video_format(self._handle,i))
        return self._Properties

    def GetInputChannels(self):
        self.InputChannels=[]
        InputChannelscount = GrabberDLL.get_input_channel_count(self._handle)
        for i in range (InputChannelscount):
            self.InputChannels.append(GrabberDLL.get_input_channel(self._handle,i))
        return self.InputChannels

    def GetVideoNormCount(self):
        self.get_video_norm=[]
        GetVideoNorm_Count=GrabberDLL.get_video_norm_count(self._handle)
        for i in range(GetVideoNorm_Count):
            self.get_video_norm.append(GrabberDLL.get_video_norm(self._handle, i))
        return self.get_video_norm


    def set_sink_format(self, Format):
        ''' set_sink_format
        Sets the pixel format in memory
        @param Format Sinkformat enumeration
        '''
        return GrabberDLL.set_sink_format(self._handle, Format.value)

    def get_sink_format(self):
        val = GrabberDLL.get_sink_format(self._handle)
        if val == 0:
            return SinkFormats.Y800
        if val == 2:
            return SinkFormats.RGB32
        if val == 1:
            return SinkFormats.RGB24
        if val == 3:
            return SinkFormats.UYVY
        if val == 4:
            return SinkFormats.Y16
        return SinkFormats.RGB24


    def start_live(self, showlive = 1):
        """
        Start the live video stream.

        showlive: 1 : a live video is shown, 0 : the live video is not shown.
        """
        Error = GrabberDLL.start_live(self._handle, showlive)
        return Error

    def stop_live(self):
        """
        Stop the live video.
        """
        Error = GrabberDLL.stop_live(self._handle)
        return Error


    def snap_image(self):
        Error = GrabberDLL.snap_image(self._handle, 2000)
        return Error


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
        iBitsPerPixel = BildDaten[2] // 8

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

    def GetFrameRate(self):
        Value = C.c_long()
        Value =  GrabberDLL.GetFrameRate(self._handle)

        return Value

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