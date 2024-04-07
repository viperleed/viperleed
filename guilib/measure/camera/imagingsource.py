"""Module imagingsource of viperleed.guilib.measure.camera.

===============================================
      ViPErLEED Graphical User Interface
===============================================

Defines the ImagingSourceCamera class. This is a concrete subclass of
CameraABC meant to handle (GigE only?) cameras from The Imaging Source.

Created: 2021-10-21
Author: Michele Riva
"""
import re
from time import perf_counter as timer
from math import ceil
from ctypes import POINTER, c_ubyte, cast as c_cast

import numpy as np
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure.camera.drivers.imagingsource import (
    ISCamera as ImagingSourceDriver, FrameReadyCallbackType,
    ImagingSourceError, SinkFormat,
    )
from viperleed.guilib.measure.camera import (
    abc,
    imagingsourcecalibration as is_calib
    )
from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.widgets.mappedcombobox import MappedComboBox


_CUSTOM_NAME_RE = re.compile(r"\[.*\]")


# TODO: test roi offset increments with models other than IMX265


# pylint: disable=useless-param-doc,useless-type-doc
# I prefer to keep the documentation of the (unused) arguments
# to remember the signature required by the driver.
@FrameReadyCallbackType
def on_frame_ready_(__grabber_handle, image_start_pixel,
                    __frame_number, process_info):
    """Prepare a frame to be processed.

    This is the callback function that will be automatically called
    whenever a new frame arrives at the PC. According to the ctypes
    documentation, this callback is executed in a separate python
    thread.

    Parameters
    ----------
    __grabber_handle : imagingsource.GrabberHandlePointer
        Unused C pointer to the grabber handle used internally
        by the Imaging Source driver to address the camera.
    image_start_pixel : ctypes.POINTER(c_ubyte)
        Pointer to the first byte of image data, i.e., the MSB??
        of the bottom-left pixel. Image data is stored from bottom
        to top and from left to right.
    __frame_number : ctypes.c_ulong
        Progressive index of current frame since the stream of
        frames was initiated (i.e., since the camera was started)
    process_info : ImageProcessInfo
        The image processing information object

    Emits
    -----
    camera.camera_busy(False)
        After enough frames have been received while estimating the
        optimal frame rate.
    camera.abort_trigger_burst()
        When camera is triggered, supports trigger burst,
        and enough frames have been received.
    camera.frame_ready(image)
        Each time, except while estimating the frame rate. Carries
        a numpy array with a copy of the image data.
    """
    camera = process_info.camera
    if not camera.connected:
        return

    process_info.frame_times.append(timer())
    n_frames_received = len(process_info.frame_times)

    if n_frames_received > 2:
        # pylint: disable=unused-variable
        # n_lost is currently unused but will be as soon as
        # the todo below is taken care of.
        n_lost, best_rate = estimate_frame_loss(process_info.frame_times,
                                                camera.get_frame_rate(),
                                                camera.exposure)
        camera.best_next_rate = best_rate

    if camera.is_finding_best_frame_rate:
        # No need to display images while estimating the best frame rate
        if n_frames_received >= camera.n_frames_estimate:
            camera.driver.stop_delivering_frames()
            camera.is_finding_best_frame_rate = False
            camera.busy = False
            # TODO: report that frames were lost, i.e., exposure
            # time can be made longer.
        return

    # Don't emit frame_ready if we got more frames than expected
    if camera.mode == 'triggered' and n_frames_received > camera.n_frames:
        return

    # Interrupt when all frames arrived. This may happen if the camera
    # was set to 'triggered' and it supports trigger burst, as we keep
    # a little safety margin (only for > 1 frame) in case frames are
    # actually lost. Notice the use of a signal: A direct call makes the
    # program access the camera from two different threads (main, and the
    # thread in which this callback runs) and silently crashes the whole
    # application.
    if (camera.mode == 'triggered'
        and camera.supports_trigger_burst
            and n_frames_received == camera.n_frames > 1):
        camera.abort_trigger_burst.emit()

    # Prepare the actual image
    image = _img_ptr_to_numpy(image_start_pixel, camera)

    # Check that the minimum intensity fits what the camera reports
    # This fixes an issue observed for some Imaging Source cameras
    # for which the lowest bit is set to zero, rather than to 1 as
    # one would expect from the documentation.
    if not camera.has_zero_minimum and not image.min():
        camera.has_zero_minimum = True

    try:
        camera.parent()
    except RuntimeError:
        # Wrapped C++ object deleted
        return

    # Send frame out
    camera.frame_ready.emit(image.copy())
# pylint: enable=useless-param-doc,useless-type-doc


def _img_ptr_to_numpy(image_start_pixel, camera):
    """Return an image as numpy array from a ctypes.POINTER.

    Parameters
    ----------
    image_start_pixel : ctypes.POINTER
        Pointer to the first byte of the first pixel of the image.
    camera : ImagingSourceCamera
        The camera that acquired the image.

    Returns
    -------
    image : numpy.ndarray
        Gray-scale image. The first index runs along rows,
        the second one along pixels within a row.

    Raises
    ------
    ValueError
        If the camera returns images with bit depth larger than 32.
    """
    # Prepare the actual image
    width, height, n_bytes, n_colors = camera.image_info
    *_, color_fmt = camera.driver.image_info

    # (1) shape: notice that height and width are swapped, since
    #            the image is stored row-by-row
    shape = (height, width, n_colors)

    # (2) data type
    if n_bytes == 1:
        dtype = np.uint8
    elif n_bytes == 2:
        dtype = np.uint16
    else:
        raise ValueError(f"Cannot work with images having {n_bytes} bytes "
                         "per color channel")

    # (3) Make a pointer to the whole array, such
    #     that the size is known to numpy later on.
    buff_size = width * height * n_colors * n_bytes * 8
    img_ptr = POINTER(c_ubyte * buff_size)

    # (4.1) create array with all color channels
    img_buffer = c_cast(image_start_pixel, img_ptr).contents
    image = np.ndarray(buffer=img_buffer, dtype=dtype, shape=shape)

    # (4.2) pick only the relevant one (i.e., green for multicolor),
    #       making the image 2D. Notice that this returns a VIEW!
    return image[:,:,color_fmt.green_channel]


def estimate_frame_loss(frame_times, frame_rate, exposure):
    """Estimate whether frames have been lost and the optimal frame rate.

    The optimal frame rate is such that frame loss
    probability is minimized.

    Parameters
    ----------
    frame_times : Sequence
        Sequence of frame arrival times.
    frame_rate : float
        Frame rate used to deliver the frames to the PC.
    exposure : float
        Exposure time in milliseconds for each frame.

    Returns
    -------
    n_lost : int
        Estimate of the number of frames that were lost. A negative value
        indicates that the current frame rate is too low. This will be
        reflected in optimal_frame_rate.
    optimal_frame_rate : float
        Maximum frame rate that will prevent frame loss.

    Raises
    ------
    ValueError
        If this function is called on a zero- or
         one-element-long sequence of frame times.
    """
    n_frames = len(frame_times)
    if n_frames <= 1:
        raise ValueError("Not enough frames to estimate frame loss!")

    # Ideally, a frame should take (slightly more) than dt_ideal
    # seconds to be delivered
    dt_ideal = 1/min(1000/exposure, frame_rate)

    # However, it can take a little longer (or shorter) to do so
    dt_total = frame_times[-1] - frame_times[0]

    # During this time interval we should have received n_total frames
    # rather than n_frames
    n_total = ceil(dt_total/dt_ideal)

    n_lost = n_total - n_frames
    opt_rate = frame_rate * n_frames/n_total

    return n_lost, opt_rate


def _color_format_mapper(formats):
    """Return keys and values from formats. Used in settings dialog."""
    return ((i.display_name, i.name) for i in formats)


# pylint: disable=too-many-public-methods
# pylint counts 33, I count 26. Of those, 17 are actually
# reimplementations of abstract methods that are only supposed
# to be used internally (i.e., they are the driver interface)
class ImagingSourceCamera(abc.CameraABC):
    """Concrete subclass of CameraABC handling Imaging Source Hardware."""

    _mandatory_settings = [
        # pylint: disable=protected-access
        # Needed for extending
        *abc.CameraABC._mandatory_settings,
        ('camera_settings', 'black_level'),
        ]

    abort_trigger_burst = qtc.pyqtSignal()                                      # TODO: could be done with QMetaObject.invokeMethod

    def __init__(self, *args, settings=None, parent=None, **kwargs):
        """Initialize instance.

        Parameters
        ----------
        *args : object
            Unused positional arguments
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
        self.hardware_supported_features.extend(['roi', 'color_format',
                                                 'black_level'])
        self.is_finding_best_frame_rate = False
        self.has_zero_minimum = False
        self.best_next_rate = 1024
        self.__black_level = -1
        self.__burst_count = -1
        self.__extra_delay = []  # Duration of abort_trigger_burst

        self.__has_callback = False

        super().__init__(ImagingSourceDriver(), *args,
                         settings=settings, parent=parent, **kwargs)
        self.abort_trigger_burst.connect(self.__on_abort_trigger_burst)

    @qtc.pyqtSlot()
    def __on_abort_trigger_burst(self):
        """Abort trigger burst."""
        # Calculate how long this takes, as this is slowing
        # down the camera, since it requires a pause/start.
        t_start = self.process_info.frame_times[-1]
        self.driver.abort_trigger_burst()
        self.__extra_delay.append(1000*(timer() - t_start))

    @property
    def calibration_tasks(self):
        """Return a dictionary of CameraCalibrationTask for self.

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
        tasks = super().calibration_tasks
        if not any(isinstance(t, is_calib.DarkLevelCalibration)
                   for t in tasks['bad_pixels']):
            tasks['bad_pixels'].append(is_calib.DarkLevelCalibration(self))
        return tasks

    @property
    def exceptions(self):
        """Return a tuple of camera exceptions."""
        return (ImagingSourceError,)

    @property
    def extra_delay(self):
        """Return the interval spent not measuring when triggered (msec)."""
        # Imaging Source cameras have significant delay only
        # when we have to interrupt an oversized trigger burst
        if (self.mode != 'triggered'
            or self.n_frames <= 1
            or not self.__extra_delay
                or not self.supports_trigger_burst):
            return 0.0
        return np.mean(self.__extra_delay)

    @property
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
        *info, _ = self.driver.image_info
        return info

    @property
    def intensity_limits(self):
        """Return the minimum and maximum value for a pixel.

        Returns
        -------
        pixel_min : int
            Minimum intensity for a pixel
        pixel_max : int
            Maximum intensity for a pixel
        """
        # TODO: here something is still incorrect: The right way to             # TODO: Probably the better way would be to have a CameraCalibrationTask that actually detects these limits from the camera by acquiring a very short dark frame and a somewhat long illuminated frame. The advantage: I would get rid of all the sensor crap, which also needs updating when new models come up.
        # do it is:
        # MAX_BIT = 8*n_bytes
        # if MAX_BIT <= DYN_RANGE --> has_zero_minimum = True
        # MIN_BIT = 0 if has_zero_minimum else MAX_BIT - DYN_RANGE
        # DELTA = 1 if has_zero_minimum else 0
        # PX_MIN = 2**MIN_BIT - DELTA
        # PX_MAX = 2**MAX_BIT - 2**(MAX_BIT - DYN_RANGE) + 2**MIN_BIT - 1

        # For 16-bit with 12-bit range and has_zero_minimum
        # this gives px_max = 2**16 - 2**4 + 2**0 - 1 = 65520

        *_, n_bytes, _ = self.image_info
        min_bit = 8*n_bytes - self.driver.dynamic_range
        if n_bytes == 1 or self.has_zero_minimum:
            pixel_min = 0
        else:
            pixel_min = 2**min_bit
        return pixel_min, 2**(8*n_bytes) - 1

    @property
    def is_running(self):
        """Return whether the camera is currently running."""
        try:
            running = self.driver.is_running
        except ImagingSourceError:
            running = False
        return running

    @property
    def black_level(self):
        """Return the black-level setting of the camera.

        The black level is a measure of a minimum photon intensity
        at pixels. Pixels illuminated with less than this intensity
        will appear in images as having self.intensity_limits[0]
        intensity. Therefore, black_level determines the lower limit
        at which image-intensity histograms are 'cut'.
        """
        try:
            black_level = self.settings.getint('camera_settings',
                                               'black_level', fallback=-2)
        except (TypeError, ValueError):
            black_level = -2

        if black_level == self.__black_level:
            return black_level

        if black_level <= -2:
            # Was not present. Let's read it from the camera
            # and store it in the settings.
            black_level = self.get_black_level()
            self.settings.set('camera_settings', 'black_level',
                              str(black_level))
            self.settings.update_file()

        _min, _max = self.get_black_level_limits()
        if black_level < _min or black_level > _max:
            base.emit_error(
                self, abc.CameraErrors.INVALID_SETTINGS,
                'camera_settings/black_level',
                f"{black_level} [out of range ({_min}, {_max})]",
                )
            return -2

        self.__black_level = black_level
        return black_level

    @property
    def color_format(self):
        """Return the the color format used."""
        color_fmt_s = self.settings.get('camera_settings', 'color_format',
                                        fallback='Y16')
        try:
            color_fmt = SinkFormat.get(color_fmt_s)
        except (ValueError, ImagingSourceError):
            # pylint: disable=redefined-variable-type
            # Probably a bug.
            base.emit_error(
                self, abc.CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                color_fmt_s, 'camera_settings/color_format', 'Y16'
                )
            color_fmt = SinkFormat.Y16
        self.settings.set('camera_settings', 'color_format', color_fmt.name)
        return color_fmt

    @color_format.setter
    def color_format(self, color_fmt):
        """Set a new color format."""
        try:
            color_fmt = SinkFormat.get(color_fmt)
        except (ValueError, ImagingSourceError):
            # pylint: disable=redefined-variable-type
            # Probably a bug.
            base.emit_error(
                self, abc.CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                color_fmt, 'camera_settings/color_format', 'Y16'
                )
            color_fmt = SinkFormat.Y16
        self.settings.set('camera_settings', 'color_format', color_fmt.name)

    @property
    def n_frames_estimate(self):
        """Return no. of frames for frame rate optimization.

        The number of frames for estimates scales with the
        exposure time such that it never takes more than 5
        seconds to perform the initial optimization.

        Returns
        -------
        n_frames : int
            Negative values indicate that there is no need
            to optimize (i.e., exposure is very long).
        """
        if self.exposure > 1000:
            return -1
        if self.exposure >= 500:
            return 4
        if self.exposure >= 200:
            return 7
        if self.exposure >= 100:
            return 11
        return 21

    @property
    def name_clean(self):
        """Return a version of .name suitable for file names."""
        # Remove completely the custom part of the name that
        # is enclosed in square brackets
        return base.as_valid_filename(_CUSTOM_NAME_RE.sub('', self.name))

    @property
    def supports_trigger_burst(self):
        """Return whether the camera allows triggering multiple frames.

        Returns
        -------
        supported : bool
            True if the camera internally supports triggering multiple
            frames, i.e., only one 'trigger_now()' call is necessary to
            deliver all the frames needed.
        """
        return "Burst Count" in self.driver.vcd_properties["Trigger"]

    def close(self):
        """Close the camera device."""
        self.driver.close()

    def get_settings_handler(self):
        """Return a SettingsHandler object for displaying settings.

        Adds:
        - 'camera_settings'/'black_level' (advanced)
        - 'camera_settings'/color_format' (advanced)

        Returns
        -------
        handler : SettingsHandler
            The handler used in a SettingsDialog to display the
            settings of this controller to users.
        """
        handler = super().get_settings_handler()

        # pylint: disable=redefined-variable-type
        # Triggered for _widget. While this is true, it is clear what
        # _widget is used for in each portion of filling the handler

        # Black level
        _widget = qtw.QSpinBox()
        _widget.setRange(*self.get_black_level_limits())
        _widget.setAccelerated(True)
        _tip = (
            "<nobr>Dark Level, Black Level, or Brightness is a measure of"
            "</nobr> minimum photon intensity at pixels. Pixels illuminated "
            "with less than this intensity will appear in images as "
            f"having minimum intensity (= {self.intensity_limits[0]}). "
            "Therefore, it determines the lower limit at which image-"
            "intensity histograms are 'cut'. The dark level is <b>optimized "
            "automatically</b> before bad pixels are identified with "
            "<b>Tools->Find bad pixels...</b>"
            )
        handler.add_option('camera_settings', 'black_level',
                           handler_widget=_widget, tooltip=_tip,
                           is_advanced=True, display_name="Dark Level")

        # Color format
        _tip = (
            "<nobr>Color format defines the bit depth of the images</nobr> "
            "acquired. It should be set to <b>monochrome</b>, and using a "
            "number of bits <b>as large as possible</b>."
            )
        _widget = MappedComboBox(_color_format_mapper)
        _widget.addItems(self.driver.video_formats_available)
        _widget.set_ = _widget.set_current_data
        _widget.get_ = _widget.currentData
        _widget.notify_ = _widget.currentIndexChanged
        handler.add_option('camera_settings', 'color_format',
                           handler_widget=_widget, tooltip=_tip,
                           is_advanced=True)
        return handler

    def list_devices(self):
        """Return a list of available devices.

        Returns
        -------
        devices : list of DeviceInfo
            Information for each of the detected Imaging Source cameras.
            For each item, only .unique_name is set, i.e., there is no
            .more information.
        """
        # Use empty dictionaries as there is no
        # additional information to pass along.
        return [base.DeviceInfo(name) for name in self.driver.devices]

    def open(self):
        """Open the camera device.

        After execution of this method the camera is ready
        to deliver frames.

        Returns
        -------
        successful : bool
            True if the device was opened successfully.
        """
        try:
            self.driver.open(self.name)
        except ImagingSourceError:
            return False

        # Without the next line the camera does not
        # work for no obvious reason (returns errors).
        self.set_roi(no_roi=True)
        self.driver.enable_auto_properties(False)
        if not self.__has_callback:
            self.set_callback(on_frame_ready_)
            self.__has_callback = True
        self.set_roi()
        return True

    @qtc.pyqtSlot()
    def get_black_level(self):
        """Return the black level currently set in the camera."""
        return self.driver.get_vcd_property("Brightness", "Value")

    @qtc.pyqtSlot()
    def set_black_level(self):
        """Set the black level in the camera from settings."""
        _level = self.black_level
        if _level < 0:  # Invalid settings
            return
        self.driver.set_vcd_property("Brightness", "Value", _level)

    def get_black_level_limits(self):
        """Return minimum and maximum values for the black level."""
        return self.driver.get_vcd_property_range("Brightness", "Value")

    def set_color_format(self):
        """Set a color format in the camera."""
        self.driver.sink_format = self.color_format
        _ = self.driver.sink_format  # update also video format

    def get_color_format(self):
        """Return the color format set in the camera."""
        return self.driver.sink_format

    def get_binning(self):
        """Return zero, i.e., unsupported binning.

        Notice that some of the Imaging Source cameras do in principle
        support 2x2 or 3x3 binning, but we want to remove hot pixels
        before binning, i.e., we need raw images back.

        Returns
        -------
        binning : {0}
        """
        return 0

    def get_exposure(self):
        """Return the exposure time (milliseconds) set in the camera."""
        return self.driver.exposure

    def set_exposure(self):
        """Set the exposure time for one frame (from settings)."""
        self.driver.exposure = self.exposure

    def get_exposure_limits(self):
        """Return the the minimum and maximum exposure times supported.

        Returns
        -------
        min_exposure, max_exposure : float
            Shortest and longest exposure times in milliseconds
        """
        return self.driver.exposure_range

    def get_frame_rate(self):
        """Return the number of frames delivered per second."""
        return self.driver.frame_rate

    def get_gain(self):
        """Get the gain (in decibel) from the camera device."""
        return self.driver.gain

    def set_gain(self):
        """Set the gain of the camera in decibel."""
        self.driver.gain = self.gain

    def get_gain_limits(self):
        """Return the the minimum and maximum gains supported.

        Make sure to self.open() before querying the camera.

        Returns
        -------
        min_gain, max_gain : float
            Smallest and largest gain factors in decibel
        """
        return self.driver.gain_range

    def get_mode(self):
        """Return the mode set in the camera.

        Returns
        -------
        mode : {'live', 'triggered'}
            The mode the camera is operating in. Continuous mode
            corresponds to receiving a 'stream' of images at
            well-defined time intervals (i.e., synchronous), while
            'triggered' is asynchronous: the camera returns a frame
            only when asked by self.trigger_now().
        """
        return 'triggered' if self.driver.trigger_enabled else 'live'

    def set_mode(self):
        """Set the camera mode from self.settings."""
        if self.mode == 'triggered':
            self.driver.trigger_enabled = True
        else:
            self.driver.trigger_enabled = False

    def get_n_frames(self):
        """Return zero as the camera does not support frame averaging."""
        return 0

    def get_n_frames_limits(self):
        """Return the minimum and maximum number of frames supported.

        Returns
        -------
        min_n_frames, max_n_frames : int
            Minimum and maximum number of frames usable for averaging
        """
        if self.supports_trigger_burst:
            return self.driver.trigger_burst_count_range
        return super().get_n_frames_limits()

    def get_roi(self):
        """Return the region of interest set in the camera.

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
        roi_x = self.driver.get_vcd_property("Partial scan", "X Offset")
        roi_y = self.driver.get_vcd_property("Partial scan", "Y Offset")
        roi_width, roi_height = self.driver.video_format_shape

        return roi_x, roi_y, roi_width, roi_height

    def set_roi(self, no_roi=False):
        """Set up region of interest in the camera (from settings).

        Parameters
        ----------
        no_roi : bool, optional
            If True, set the ROI to the full size of the sensor rather
            than using the value from the settings. Default is True.

        Returns
        -------
        None.
        """
        if no_roi:
            roi = (0, 0, *self.driver.video_format_shape_range[2:])
        else:
            roi = self.roi
        roi_x, roi_y, roi_width, roi_height = roi
        color = self.color_format

        self.driver.set_video_format(
            f"{color.name} ({roi_width}x{roi_height})"
            )
        self.driver.set_vcd_property("Partial scan", "X Offset", roi_x)
        self.driver.set_vcd_property("Partial scan", "Y Offset", roi_y)

    def get_roi_size_limits(self):
        """Return minimum, maximum and granularity of the ROI.

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
        roi_range = self.driver.video_format_shape_range
        roi_increments = self.driver.video_format_shape_increments
        return roi_range[:2], roi_range[2:], roi_increments, (2, 2)

    def start_frame_rate_optimization(self):
        """Start estimation of the best frame rate."""
        # Connect the busy signal here. The callback
        # takes care of making the camera not busy
        # when done with the estimate.
        base.safe_connect(self.camera_busy, self.__start_postponed,
                          type=qtc.Qt.UniqueConnection)

        self.is_finding_best_frame_rate = True
        self.process_info.clear_times()

        self.preparing.emit(True)
        self.busy = True

        # Begin delivering frames at maximum speed.
        # The on_frame_ready_ callback will take care of
        # stopping whenever enough frames arrived and
        # will set the optimal frame rate.
        self.driver.start_delivering_frames(readout_rate=1024)

    def reset(self):
        """Reset the camera to default mode."""
        try:
            # Not actually sure this does reset, but does not
            # matter much: Seems it always returns ERROR.
            # Wrapped in try...except to swallow it.
            self.driver.reset()
        except ImagingSourceError:
            pass

    def set_callback(self, on_frame_ready):
        """Pass a frame-ready callback to the camera driver.

        Parameters
        ----------
        on_frame_ready : FrameReadyCallbackType
            The function that will be called by the camera each time
            a new frame arrives. Takes self.process_info as the last
            argument.
        """
        self.driver.set_frame_ready_callback(on_frame_ready, self.process_info)

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(object)
    def start(self, *_):
        """Start the camera in self.mode."""
        self.n_frames_done = 0
        self.__extra_delay = []

        if self.n_frames_estimate > 0:
            # Optimize now the frame rate to minimize frame losses.
            # Postpone actual starting to after optimization is over.
            # Once done, the camera will be automatically started with
            # a call to __start_postponed() (connected in the next
            # call to the camera_busy signal).
            self.start_frame_rate_optimization()
        else:
            self.best_next_rate = 1024
            self.__start_postponed()

    @qtc.pyqtSlot(bool)
    def __start_postponed(self, *_):
        """Actually start camera after frame-rate estimate is over.

        This method is connected to the camera_busy signal right
        before the camera starts estimating the best frame rate.

        Parameters
        ----------
        *_ : object
            Unused arguments. Necessary to allow connection to
            the camera_busy signal.

        Returns
        -------
        None.
        """
        if self.is_finding_best_frame_rate:
            return

        base.safe_disconnect(self.camera_busy, self.__start_postponed)

        # Call base that starts the processing thread if needed
        super().start()
        if abs(self.best_next_rate - self.driver.frame_rate) > 1:
            self.driver.frame_rate = self.best_next_rate
        self.process_info.clear_times()
        mode = self.mode if self.mode == 'triggered' else 'continuous'
        self.driver.start(mode)
        self.preparing.emit(False)
        self.started.emit()

    @qtc.pyqtSlot()
    def stop(self):
        """Stop the camera."""
        if not super().stop():
            # No need to stop, or cannot stop yet
            return False

        # Now actually stop the driver.
        # One could wrap the next line in if self.is_running,
        # but it seems that there is a bug in the driver that
        # can return True if the device was lost.
        try:
            self.driver.stop()
        except ImagingSourceError:
            pass

        self.stopped.emit()
        return True

    @qtc.pyqtSlot()
    def trigger_now(self):
        """Start acquiring one (or more) frames now.

        Returns
        -------
        successfully_triggered : bool
            True if the camera was successfully triggered.

        Emits
        -----
        error_occurred(CameraErrors.UNSUPPORTED_OPERATION)
        """
        if not super().trigger_now():
            return False
        if abs(self.best_next_rate - self.driver.frame_rate) > 1:
            self.driver.frame_rate = self.best_next_rate
            # Store the value currently used, to avoid trying to
            # set it multiple times if we're already at maximum
            self.best_next_rate = self.driver.frame_rate
        if self.supports_trigger_burst:
            burst_count = self.n_frames
            if burst_count > 1:
                # If more than one trigger frame, keep a little
                # safety margin, in case some frames are lost.
                # NB: the estimate of the best frame rate seems
                # to work pretty well, so only 20% more frames
                # should be enough.
                # TODO: does it make sense to scale the excess with
                #       how far the frame rate is from the maximum?
                burst_count = ceil(burst_count * 1.2)
            if burst_count != self.__burst_count:
                self.driver.trigger_burst_count = burst_count
                self.__burst_count =  self.driver.trigger_burst_count
        self.driver.send_software_trigger()
        return True
# pylint: enable=too-many-public-methods
