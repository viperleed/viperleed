# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 16:26:14 2020

@author: e-thi
"""

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


class Camera:
    def __init__(self, *args, **kwargs):
        self.callback_data = CallbackData()
        
    @staticmethod
    def camera_from_ini(camera_class_name):
        """
        Given a camera class name as a string, returns an instance of the class
        if implemented. Raises ValueError otherwise.
        """
        try:
            camera_class = getattr(sys.modules[__name__], camera_class_name)
        except AttributeError:
            raise ValueError(f"camera class {camera_class_name} not found")
        return camera_class()
    
    def initialize(self, *args, **kwargs):
        """
        Initialize camera
        """
        raise NotImplementedError()
    
    def is_busy(self):
        return self.callback_data.busy
    
    def set_exposure(self, time):
        """
        Parameters
        ----------
        time: float
              time for one frame in milliseconds
        """
        raise NotImplementedError()
    
    def set_gain(self, gain):
        """
        Parameters
        ----------
        gain: float
              gain in decibel
        """
        raise NotImplementedError()
    
    def set_n_frames(self, n_frames):
        """
        Parameters
        ----------
        n_frames: int, or str that can be cast to int
                  number of frames to be averaged
        """
        if isinstance(n_frames, str):
            n_frames = int(n_frames)
        if n_frames <= 0:
            n_frames = 1
        
        self.callback_data.process_params['n_frames'] = n_frames
    
    def enable_trigger(self, active):
        """
        Parameters
        ----------
        active: bool
                False for "live" mode, True for single exposures
        """
        raise NotImplementedError()
        
    def set_continuous_mode(self, active):
        """
        Parameters
        ----------
        active: bool
                False for "live" mode, True for single exposures
        """
        raise NotImplementedError()
    
    def trigger_now(self):
        """
        Before triggering, set enable_trigger(True)
        """
        raise NotImplementedError()
    
    def set_roi(self, x0, y0, width, height):
        """
        Set up region of interest. If the camera does not support setting a
        ROI, the code for setting the ROI should go in the callback function
        that processes the images
        
        Parameters
        ----------
        x0, y0: int
                top-left corner
        width, height: ints
        """
        raise NotImplementedError()
    
    def set_binning(self, bin_factor):
        """
        Set up binning in the camera. If the camera does not support binning,
        the code for binning can be inserted in the callback function that
        processes the images
        
        Parameters
        ----------
        bin_factor: int
                    binning will be executed on a (bin_factor, bin_factor) mesh
        """
        raise NotImplementedError()
    
    def set_callback(self, user_function):
        """
        Parameters
        ----------
        user_function: callable
              function that will be called in triggered mode after a frame
              exposure time has elapsed
        """
        raise NotImplementedError()
    
    def snap_image(self):
        """
        In live mode (i.e., self.enable_trigger(False)) returns a snapshot of
        the current frame as a numpy.ndarray
        """
        raise NotImplementedError()
    
    def set_imagepath(self):
        """
        Sets the directory to save the images
        """
        raise NotImplementedError()
    
    def start_camera(self, live_mode_active):
        """
        Starts the camera with live mode (True) or starts the camera without 
        live mode (False)
        Parameters: 
        ----------
        live_mode_active : boolean 
        """
        raise NotImplementedError()
    
    def stop_camera(self):
        """
        Stops the camera and resets all the camera properties to their default
        values. If a property has automation, the automatic will be enabled. WHAT DOES THIS MEAN??
        """
        raise NotImplementedError
        
from drivers import * #has to be down here, otherwise I get circle dependecies(?)

# if __name__ == "__main__":
# #    ImagingSourceDMKCamera()
#     camera_instance = Camera.camera_from_ini("ImagingSourceDMKCamera")


