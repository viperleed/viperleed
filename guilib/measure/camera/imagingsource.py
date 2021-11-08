"""Module imagainsource of viperleed.guilib.measure.camera.

===============================================
      ViPErLEED Graphical User Interface
===============================================

Defines the ImagingSourceCamera class. This is a concrete subclass of
CameraABC meant to handle (GigE only?) cameras from The Imaging Source.

Created: 2021-10-21
Author: Michele Riva
"""
from time import perf_counter as timer
from math import ceil
from ctypes import POINTER, c_ubyte, cast as c_cast

import numpy as np
from PyQt5 import QtCore as qtc

from viperleed.guilib.measure.camera.drivers.imagingsource import (
    ISCamera as ImagingSourceDriver, FrameReadyCallbackType,
    ImagingSourceError, SinkFormat,
    )
from viperleed.guilib.measure.camera.abc import CameraABC, CameraErrors
from viperleed.guilib.measure.camera.imageprocess import ImageProcessInfo
from viperleed.guilib.measure.hardwarebase import emit_error

@FrameReadyCallbackType
def on_frame_ready(__grabber_handle, image_start_pixel,
                   frame_number, process_info):
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
    frame_number : ctypes.c_ulong
        Progressive index of current frame since the stream of
        frames was initiated (i.e., since the camera was started)
    process_info : ImageProcessInfo
        The image processing information object

    Returns
    -------
    None.

    Emits
    -----
    camera.done_estimating_frame_rate
        After enough frames have been received while estimating the
        optimal frame rate.
    camera.frame_ready
        Each time, except while estimating the frame rate. Carries
        a numpy array with a copy of the image data.
    """
    process_info.frame_times.append(timer())
    camera = process_info.camera
    n_frames_received = len(process_info.frame_times)

    if n_frames_received > 2:
        n_lost, best_rate = estimate_frame_loss(process_info.frame_times,
                                                camera.driver.frame_rate,
                                                camera.exposure)
        camera.best_next_rate = best_rate

    if camera.is_finding_best_frame_rate:
        # No need to display images while estimating the best frame rate
        if n_frames_received >= camera.n_frames_estimate:
            camera.driver.stop_delivering_frames()
            camera.is_finding_best_frame_rate = False
            camera.done_estimating_frame_rate.emit()
        return

    # This may happen if the camera was set to 'triggered' and it
    # supports trigger burst, as we keep a little safety margin
    # in case frames are actually lost. Notice the use of a signal:
    # A direct call makes the program access the camera from two
    # different threads (main, and the thread in which this callback
    # runs) and crashes the whole application.
    if (camera.supports_trigger_burst
            and n_frames_received >= camera.n_frames):
        camera.abort_trigger_burst.emit()

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

    # (3) Make pointer to the whole array, such that the size is known
    #     to numpy later on.
    buff_size = width * height * n_colors * n_bytes * 8
    img_ptr = POINTER(c_ubyte * buff_size)

    # (4.1) create array with all color channels
    img_buffer = c_cast(image_start_pixel, img_ptr).contents
    image = np.ndarray(buffer=img_buffer, dtype=dtype, shape=shape)

    # (4.2) pick only the relevant one (i.e., green for multicolor),
    #       making the image 2D.
    image = image[:,:,color_fmt.green_channel]

    camera.frame_ready.emit(image.copy())


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


class ImagingSourceCamera(CameraABC):
    """Concrete subclass of CameraABC handling Imaging Source Hardware."""

    abort_trigger_burst = qtc.pyqtSignal()
    done_estimating_frame_rate = qtc.pyqtSignal()

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
        self.hardware_supported_features.append('color_format')
        self.is_finding_best_frame_rate = False
        self.best_next_rate = 1024

        self.__has_callback = False

        super().__init__(ImagingSourceDriver(), *args,
                         settings=settings, parent=parent, **kwargs)
        try:
            self.done_estimating_frame_rate.connect(self.__start_postponed,
                                                    qtc.Qt.UniqueConnection)
        except TypeError:
            # Already connected
            pass
        self.abort_trigger_burst.connect(self.driver.abort_trigger_burst)

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
    def is_running(self):
        """Return whether the camera is currently running."""
        try:
            running = self.driver.is_running
        except ImagingSourceError:
            running = False
        return running

    @property
    def color_format(self):
        """Return the the color format used."""
        color_fmt = self.settings.get('camera_settings', 'color_format',
                                      fallback='Y16')
        try:
            color_fmt = SinkFormat.get(color_fmt)
        except (ValueError, ImagingSourceError):
            emit_error(self, CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                       color_fmt, 'camera_settings/color_format', 'Y16')
            color_fmt = SinkFormat.Y16
        self.settings.set('camera_settings', 'color_format', color_fmt.name)
        return color_fmt

    @color_format.setter
    def color_format(self, color_fmt):
        """Set a new color format."""
        try:
            color_fmt = SinkFormat.get(color_fmt)
        except (ValueError, ImagingSourceError):
            emit_error(self, CameraErrors.INVALID_SETTING_WITH_FALLBACK,
                       color_fmt, 'camera_settings/color_format', 'Y16')
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
        if self.exposure >= 1500:
            return -1
        if self.exposure >= 500:
            return 4
        if self.exposure >= 200:
            return 7
        if self.exposure >= 100:
            return 11
        return 21

    @property
    def supports_trigger_burst(self):
        """Return whether the camera allows triggering multiple frames.

        This property should be reimplemented in concrete subclasses.

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

    def list_devices(self):
        """Return a list of available device names.

        Returns
        -------
        devices : list
            Each element is a string representing
            the name of a camera device.
        """
        return self.driver.devices

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
            self.set_callback(on_frame_ready)
            self.__has_callback = True
        self.set_roi()

        return True

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
        return super().get_binning()

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
        return super().get_n_frames()

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
        x_center = self.driver.get_vcd_property("Partial scan", "X Offset")
        y_center = self.driver.get_vcd_property("Partial scan", "Y Offset")
        roi_width, roi_height = self.driver.video_format_shape

        roi_x = round(x_center - roi_width/2)
        roi_y = round(y_center - roi_height/2)

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
            *_, roi_width, roi_height = self.driver.video_format_shape_range
            x_center, y_center = roi_width/2, roi_height/2
        else:
            roi_x, roi_y, roi_width, roi_height = self.roi
            x_center, y_center = (roi_x + roi_width/2, roi_y + roi_height/2)
        x_center, y_center = round(x_center), round(y_center)
        color = self.color_format

        self.driver.set_video_format(
            f"{color.name} ({roi_width}x{roi_height})"
            )
        self.driver.set_vcd_property("Partial scan", "X Offset", x_center)
        self.driver.set_vcd_property("Partial scan", "Y Offset", y_center)

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
        """
        roi_range = self.driver.video_format_shape_range
        roi_increments = self.driver.video_format_shape_increments
        return roi_range[:2], roi_range[2:], roi_increments

    def start_frame_rate_optimization(self):
        """Start estimation of the best frame rate."""
        self.is_finding_best_frame_rate = True
        self.process_info.clear_times()

        # Begin delivering frames at maximum speed.
        # The on_frame_ready callback will take care of
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

        Returns
        -------
        None
        """
        self.driver.set_frame_ready_callback(on_frame_ready,
                                             self.process_info)

    def start(self, *_):
        """Start the camera in self.mode."""
        self.n_frames_done = 0

        if self.n_frames_estimate > 0:
            # Optimize now the frame rate to minimize frame losses.
            # Postpone actual starting to after optimization is over.
            # Once done, the camera will be automatically started with
            # a call to __start_postponed() (connected in __init__ to
            # the done_estimating_frame_rate signal).
            # If we ever decide to not do this, we can simply replace
            # the call above with a direct call to __start_postponed()
            self.start_frame_rate_optimization()
        else:
            self.best_next_rate = 1024
            self.done_estimating_frame_rate.emit()

    def __start_postponed(self):
        # Call base that starts the processing thread if needed
        super().start()
        self.driver.frame_rate = self.best_next_rate
        self.process_info.clear_times()
        mode = self.mode if self.mode == 'triggered' else 'continuous'
        self.driver.start(mode)

    def stop(self):
        """Stop the camera."""
        self.driver.stop()
        super().stop()

    def trigger_now(self):
        """Start acquiring one (or more) frames now.

        Returns
        -------
        None.

        Emits
        -----
        error_occurred(CameraErrors.UNSUPPORTED_OPERATION)
        """
        super().trigger_now()
        self.driver.frame_rate = self.best_next_rate
        if self.supports_trigger_burst:
            trigger_burst = self.n_frames
            if trigger_burst > 1:
                # If more than one trigger frame, keep a little
                # safety margin, in case some frames are lost.
                # NB: the estimate of the best frame rate seems
                # to work pretty well, so only 20% more frames
                # should be enough.
                trigger_burst = ceil(trigger_burst * 1.2)
            self.driver.trigger_burst_count = trigger_burst
        self.driver.send_software_trigger()
