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

    def __init__(self, *args, **kwargs):
        """Initialize instance."""
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
    def set_exposure(self, new_exposure_time):
        """Set the exposure time for one frame.

        Parameters
        ----------
        new_exposure_time : float
            Expoure time for one frame in milliseconds
        """
        return

    @abstractmethod
    def set_gain(self, new_gain):
        """Set the gain of the camera in decibel.

        Parameters
        ----------
        new_gain : float
            Gain in decibel
        """
        return

    @property
    def n_frames(self):
        """Return the number of frames used for averaging."""
        return self.callback_data.process_params['n_frames']

    @n_frames.setter
    def n_frames(self, n_frames):
        """Set the number of frames for averaging.

        Each frame will be acquired for .exposure_time
        milliseconds.

        Parameters
        ----------
        n_frames : int or str
            When a string, it must represent an integer.
            The number of averaging frames is set to one in
            case n_fames is a non-positive integer.
        """
        if isinstance(n_frames, str):
            n_frames = int(n_frames)
        if n_frames <= 0:
            n_frames = 1

        self.callback_data.process_params['n_frames'] = n_frames

    @abstractmethod
    def enable_trigger(self, active):
        """Enable triggering of the camera.

        When triggering is off, the camera should be in
        live-view mode, i.e., returning a stream of images.

        Parameters
        ----------
        active : bool
            False for "live" mode, True for single exposures
        """
        return

    @abstractmethod
    def set_continuous_mode(self, active):
        """
        TODO: how is this any different from the previous one.

        Parameters
        ----------
        active: bool
                False for "live" mode, True for single exposures
        """
        return

    @abstractmethod
    def trigger_now(self):
        """Start acquiring one (or more) frames now.

        Before triggering, set enable_trigger(True)
        """
        return

    @abstractmethod
    def set_roi(self, x0, y0, width, height):
        """Set up region of interest.
        If the camera does not support setting internally a ROI, the
        the ROI setting should be used in the callback function that
        processes the images to crop the images to the ROI.

        Parameters
        ----------
        x0, y0 : int
            Pixel coordinates of the top-left corner.
            (0, 0) is the most top-left pixel of the camera.
        width, height: int
            Width and height (in pixels) of the ROI, before binning.
        """
        return

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


