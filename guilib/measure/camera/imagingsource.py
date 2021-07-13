# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 18:00:17 2020

@author: e-thi
"""
import os
import time
import ctypes

import cv2
# NON STANDARD, needs to be tested against speed of image saving with respect to
# skimage.io, gdal (part of osgeo?), rasterio, perhaps tiffile?, and pillow-simd

import numpy as np

# rather add the directories to sys.path. Also, the generic camera libraries
# should be reorganized into a better directory tree:
# ./
# config.ini goes here, no need to have a folder for one single file
# ./drivers
#   camera.py
#   /imagingsource
#    imagingsource.py
#    (all the current contents of Camera_libraries, excl. camera.py)
# ./temp (folder to create at runtime, and to be cleaned up afterwards;
#         path is hard-coded here, not included in .ini file)
#   all measurement files go here, before packing to a .zip
# ./LEED_data (standard folder for saving the .zip files)
os.chdir('Camera_libraries')

from camera import CameraABC
from camera import bin_image
import tisgrabber as IC

os.chdir('..')


def callback(hGrabber, pBuffer, framenumber, user_data):
    """
    This is the Callback-function, which is called when a software Trigger is
    send to the camera. Main functions are importing the image from the "
    "buffer, do some image post processing (e.g. binning) and saving to the "
    "harddisc"

    Parameters
    ----------
    hGrabber: This is the real pointer to the grabber object. Do not use.
    pBuffer : Pointer to the first pixel's first byte
    framenumber : Number of the frame since the stream started
    user_data : Pointer to user data structure in 'camera.py'
    """
    Imageformat = user_data.camera_lib.GetImageDescription()[:4]
    buffersize = Imageformat[0] * Imageformat[1] * Imageformat[2]
    if buffersize > 0:
        # WHY IS THE NEXT LINE EVEN NEEDED? pBuffer should already be a pointer
        # to a byte array. Probably best to use (copy avoids memory overlaps):
        # np.ctypeslib.as_array(ptr, shape=(8,)).copy()
        image = ctypes.cast(pBuffer,
                            ctypes.POINTER(ctypes.c_ubyte * buffersize))
        image_name = (str(user_data.process_params['filename_prefix'])+
                      '_' +
                      str(user_data.process_params['filename_energy']) +
                      'eV.tif')
        # BUG: user_data.image is never initialized to a buffer of zeros!!
        #      -> with multiple frames and the summation for averaging it sums
        #         random values in memory with the first image. The
        #         initialization of the image array should be done in the camera
        #         class, probably as part of the self.initialize(). Also,
        #         the buffer should be reset to zero after an image is saved,
        #         and each time the video format or ROI is changed
        if user_data.process_params['n_frames'] > 1:
            array_data_type = np.int16
            # array_data_type = np.int32 #if saved on this, python aborts?
            user_data.image += np.uint32(
                np.ndarray(buffer = image.contents,
                           dtype = array_data_type,
                           shape = (Imageformat[1], # Width WRONG
                                    Imageformat[0])))  # Height WRONG
            if (user_data.process_params['counter_frames']
                < user_data.process_params['n_frames']):
                user_data.busy = False
        else:
            user_data.dac_busy = False
            array_data_type = np.uint16
            user_data.image = np.ndarray(buffer = image.contents,
                                          dtype = array_data_type,
                                          shape = (Imageformat[1], # Width
                                                   Imageformat[0]))  # Height
        if (user_data.process_params['counter_frames'] >=
            user_data.process_params['n_frames']):
            user_data.dac_busy = False
            user_data.process_params['counter_frames'] = 0
            if user_data.process_params['n_frames'] > 1:
                user_data.image = (user_data.image /
                                   user_data.process_params['n_frames'])
            if user_data.process_params['binning'] is not None:
                user_data.image = bin_image(
                    user_data.image,
                    user_data.process_params['binning'])
            user_data.busy = False
            user_data.image = np.uint16(user_data.image)
            os.chdir(user_data.process_params['directory'])
            cv2.imwrite(image_name, user_data.image)
            os.chdir('..')
    else:
        print("Error: There is nothing in the Buffer!")
        user_data.busy = False


class ImagingSourceDMKCamera(CameraABC):
    """Class to control a DMK-model camera from Imaging Source."""

    PROP_CAM_AUTO_EXPOSURE = 4
    ENABLE = 1
    DISABLE = 0
    SUCCESS = 1

    RESOLUTION_MAX_WIDTH = 2048
    RESOLUTION_MIN_WIDTH = 528
    RESOLUTION_MULTIPLE_WIDTH = 16
    RESOLUTION_MAX_HEIGHT = 1536
    RESOLUTION_MIN_HEIGHT = 4
    RESOLUTION_MULTIPLE_HEIGHT = 4

    def __init__(self):
        super().__init__()
        self.camera_lib = IC.TIS_CAM()
        self.callback_data.camera_lib = self.camera_lib
        self.callback_function = IC.TIS_GrabberDLL.FRAMEREADYCALLBACK(callback)

    def fail(self, funct, *args, **kwargs):
        """Run self.camera_lib.funct(*args, **kwargs).
        
        Returns
        -------
        fail : bool
            True if the action failed
        """
        return getattr(self.camera_lib, funct)(*args, **kwargs) != self.SUCCESS

    def initialize(self, camera_config):
        
        # TODO: the next lines should rather go into a dedicated
        #       portion of the GUI where stuff is loaded from
        #       the config file. Part of this can go to a general
        #       .settings property of the ABC.
        imagingsource_devices = self.camera_lib.GetDevices()
        device_name = camera_config['device_name']

        # Check whether there is a camera with the device name from the config
        # file
        indices = [i for i, s in enumerate(imagingsource_devices)
                     if device_name in s.decode()]
        if len(indices) == 0:
            raise RuntimeError("No Camera with Name \"%s\" has been found in "
                               "Devicelist!" % str(device_name))
        elif len(indices) > 1:
            print("WARNING: There are more devices with the name \"%s\", program "
                  "takes the first one..." % str(device_name))
        # END PART TO MOVE

        # try connecting
        if self.fail('openVideoCaptureDevice', camera_config['device_name']):   # device_name mandatory
            raise RuntimeError("Impossible to connect to Camera " +
                               device_name)

        # set up the video format
        video_format = camera_config['video_format'].split()[0]                 # video_format seems IS specific --> super()
        try: 
            video_format = getattr(IC.SinkFormats, video_format)
        except AttributeError as err:
            raise RuntimeError(
                f"Invalid video format {video_format}"
                )
            
        if self.fail('SetVideoFormat', video_format):
            raise RuntimeError("Error setting format" + video_format)
        # if 'Y16' in video_format:
            # if self.fail('SetFormat', IC.SinkFormats.Y16):
                # raise RuntimeError("Error setting format Y16")

        # set up the region of interest
        if camera_config['set_roi'] == 'True':                                  # roi optional, stored in config as x, y, w, h, fallback to None
            self.set_roi(camera_config['roi_x'], camera_config['roi_y'])

        # set binning
        self.set_binning(camera_config['bin_factor'])                           # binning optional fallback to 1

        # set exposure, gain, number of frames, frame rate,image path,
        # nameprefix
        self.set_exposure(camera_config['exposure_time'])                       # exposure_time mandatory
        self.set_gain(camera_config['gain'])                                    # gain mandatory?
        self.set_n_frames(camera_config['n_frames'])                            # optional, fallback to 1
        self.set_framerate(camera_config['frame_rate'])                         # probably skip this altogether and stick to exposure_time only
        self.set_imagepath(camera_config['image_path'],                         # useless
                           camera_config['file_name'])

        self.set_callback()                                                     # probably IS specific
        self.set_continuous_mode(camera_config['live_mode'])                    # mandatory
        self.enable_trigger(camera_config['live_mode'])                         # unclear how this differs from continuous on/off

        self.start_camera(camera_config['live_mode'])
        print('Camera Initialization: DONE!')

    def set_exposure(self, time_exposure):
        """
        Parameters
        ----------
        time: float, or string that can be cast to float
              time for one frame in milliseconds
        """
        if isinstance(time_exposure, str):
            time_musecs = float(time_exposure)
        else:
            time_musecs = time_exposure
        time_musecs *= 1e3
        time_musecs = int(time_musecs)

        time_min = [0]
        time_max = [0]
        if self.fail('GetPropertyValueRange', 'Exposure', 'Value',
                     time_min, time_max):
            raise RuntimeError("Failed to receive the exposure time range from camera.")

        if time_musecs < time_min[0] or time_musecs > time_max[0]:
            raise ValueError("Camera doesn't support an exposure time of %f "
                             " milliseconds. Supported exposure times are from "
                             " %i to %i milliseconds" %
                             (float(time_exposure),
                              int(time_min[0]*1e-3),
                              int(time_max[0]*1e-3)))

        if self.fail('enableAutoCameraProperty', self.PROP_CAM_AUTO_EXPOSURE,
                     self.DISABLE):
            raise RuntimeError("Failed disabling auto-exposure")
        if self.fail('SetPropertyValue', "Exposure", "Value", time_musecs):
            raise ValueError("Failed setting exposure to %s ms" % str(time))

    def set_gain(self, gain):
        """
        Parameters
        ----------
        gain: float
              gain in decibel
        """
        if isinstance(gain, str):
            gain_db = round(float(gain))
        else:
            gain_db = round(gain)

        gain_min = [0]
        gain_max = [0]
        if self.fail('GetPropertyValueRange', 'Gain', 'Value',
                     gain_min, gain_max):
            raise RuntimeError("Failed to receive the gain range from camera.")

        if gain_db < gain_min[0] or gain_db > gain_max[0]:
            raise ValueError("Camera doesn't support gain of %i. The camera "
                             "supports gain from %i to %i" %
                             (gain_db,
                              int(gain_min[0]),
                              int(gain_max[0])))

        if self.fail('SetPropertySwitch', "Gain", "Auto", 0):
            raise RuntimeError("Failed disabling auto-gain")

        if self.fail('SetPropertyValue', "Gain", "Value", gain_db):
            raise ValueError("Failed setting gain to %s dB" % str(gain))

    def set_framerate(self, rate):
        """
        Parameters
        ----------
        rate: float
              frame rate in frames per second
        """
        if isinstance(rate, str):
            rate_fps = float(rate)
        else:
            rate_fps = rate
        if self.fail('SetFrameRate', rate_fps):
            raise ValueError("Failed setting frame rate to %s fps" % str(rate))

    def enable_trigger(self, active):
        """
        Parameters
        ----------
        active: string
                False for "live" mode, True for single exposures
        """
        if active == 'False':
            if self.fail('EnableTrigger', self.ENABLE):
                raise ValueError("Failed to enable trigger mode")
        else:
            if self.fail('EnableTrigger', self.DISABLE):
                raise ValueError("Failed to disable trigger mode")

    def set_continuous_mode(self, active):  # THIS IS COUNTERINTUITIVE
        """
        If live mode is on, the Function needs 'False'.
        If this function is set to 'True' the ContinuousMode is set to DISABLE.
        This means, that every trigger the image is processed automatically in
        the callback function. It is needed for the trigger mode. Otherwise the
        images has to be handled separately each snap.

        Parameters
        ----------
        active: string
                False for "live" mode, True for single exposures
        """
        if active == 'False':
            if self.fail('SetContinuousMode', self.DISABLE):
                raise ValueError("Failed to disable continuous mode")
        else:
            if self.fail('SetContinuousMode', self.ENABLE):
                raise ValueError("Failed to enable continuous mode")

    def trigger_now(self):
        """
        Before triggering, set enable_trigger(True) and
        set_continuous_mode(True).

        Returns:
        ----------
        True: when trigger is sent successfully
        False: when camera is still busy or image post processing in the
               callback is still not finished
        """
        while(self.callback_data.busy):
            pass  # this is super resource-intensive as it will clog the CPU
        while(self.callback_data.process_params['counter_frames'] <
              self.callback_data.process_params['n_frames']):
            if(not self.callback_data.busy):
                self.callback_data.busy = True #For next trigger
                self.callback_data.dac_busy = True #For next dac value
                if self.fail('SoftwareTrigger'):
                    self.callback_data.busy = False
                    raise RuntimeError("Sending software trigger failed.")
                self.callback_data.process_params['counter_frames'] += 1
            else:
                print("Camera busy, next try...")
                time.sleep(1)

    # this does not have the correct signature from the camera class! It should
    # take (and set) also the width and height of the ROI
    def set_roi(self, x0, y0):
        """
        Set up region of interest. Here it disables the auto-center and moves
        the frame from the top-left corner. Allowed values for x (horizontal)
        are from 0 - 1520 and for y from 0 - 1532.
        region of interest.

        Parameters
        ----------
        x0, y0: ints; in pixel
                top-left corner
        """

        # width = int(re.findall(r'\d+', video_format)[1])
        # height = int(re.findall(r'\d+', video_format)[2])

        if isinstance(x0, str):
            x0_pixel = int(x0)
        elif isinstance(y0, str):
            y0_pixel = int(y0)
        else:
            x0_pixel = x0
            y0_pixel = y0
            # width_pixel = width
            # height_pixel = height

        # if (width_pixel < self.RESOLUTION_MIN_WIDTH or
        #     width_pixel > self.RESOLUTION_MAX_WIDTH or
        #     width_pixel > self.RESOLUTION_MAX_WIDTH or
        #     height_pixel < self.RESOLUTION_MIN_HEIGHT or
        #     height_pixel > self.RESOLUTION_MIN_HEIGHT):

            # raise ValueError("Camera doesn't support Video Format, Width "
            #                    "of Video has to be in the Range of %i and %i "
            #                    "and has to be a multiple of %i, the height "
            #                    "of the video has to be in the range of "
            #                    "%i and %i and has to be a muliple of %i." %
            #                    (self.RESOLUTION_MIN_WIDTH,
            #                     self.RESOLUTION_MAX_WIDTH,
            #                     self.RESOLUTION_MULTIPLE_WIDTH,
            #                     self.RESOLUTION_MIN_HEIGHT,
            #                     self.RESOLUTION_MAX_HEIGHT,
            #                     self.RESOLUTION_MULTIPLE_HEIGHT))

        min_offset = [0]
        max_offset = [0]
        if self.fail('GetPropertyValueRange', 'Partial scan', 'X Offset',
                     min_offset, max_offset):
            raise RuntimeError("Camera failed to return X Offset.")

        if x0_pixel < int(min_offset[0]) or x0_pixel > int(max_offset[0]):
            raise ValueError("Camera doesn't support %i pixel offset in "
                             "horizontal direction, allowed offset in "
                             "horizontal direction is from %i to %i pixel" %
                             (x0_pixel,
                              int(min_offset[0]),
                              int(max_offset[0])))

        if self.fail('GetPropertyValueRange', 'Partial scan', 'Y Offset',
                     min_offset, max_offset):
            raise RuntimeError("Camera failed to return Y Offset.")

        if y0_pixel < min_offset[0] or y0_pixel > max_offset[0]:
            raise ValueError("Camera doesn't support %i pixel offset in "
                             "horizontal direction, allowed offset in "
                             "horizontal direction is from %i to %i"
                             % (y0_pixel, min_offset[0],max_offset[0]))

        # if self.camera_lib.SetVideoFormat(video_format) != self.SUCCESS:
        #     raise RuntimeError("Error setting format" + video_format)
        # if 'Y16' in video_format:
        #     if self.camera_lib.SetFormat(
        #             IC.SinkFormats.Y16) != self.SUCCESS:
        #         raise RuntimeError("Error setting format Y16")

        if self.fail('SetPropertySwitch', 'Partial scan', 'Auto-center', 0):
            raise RuntimeError("Error disabling Auto-center")

        if self.fail('SetPropertyValue', 'Partial scan', 'X Offset', x0_pixel):
            raise RuntimeError("Error setting horizontal offset of "
                               + x0_pixel)

        if self.fail('SetPropertyValue', 'Partial scan', 'Y Offset', y0_pixel):
            raise RuntimeError("Error setting horizontal offset of "
                               + y0_pixel)

    def set_binning(self, bin_factor):
        """
        Set up binning in the camera. Here it post processes the image in the
        Callback function after the image is taken by the trigger

        Parameters
        ----------
        bin_factor: int
                    binning will be executed on a (bin_factor, bin_factor) mesh
        """
        if isinstance(bin_factor, str) :
            bin_factor_int = int(bin_factor)
        else:
            bin_factor_int = int(bin_factor)
        if bin_factor_int > 0:
            self.callback_data.process_params['binning'] = bin_factor_int

    def set_callback(self):
        """
        Parameters
        ----------
        user_function: callable
              function that will be called in triggered mode after a frame
              exposure time has elapsed
        """
        if self.fail('SetFrameReadyCallback',
                     self.callback_function, self.callback_data):
            raise ValueError("Failed to pass the function pointer and the user"
                             "data to the library")

    def snap_image(self):
        """
        In live mode (i.e., self.enable_trigger(False)) returns a snapshot of
        the current frame as a numpy.ndarray, when snap was successful, or
        False when not successful

        Returns:
        ----------
        image: np.ndarray, dtype=np.uint16
        boolean: False
        """
        if self.fail('SnapImage'):
            return False
        image = self.camera_lib.GetImageEx()  # THIS IS NEVER RETURNED!!
        cv2.imshow('Window', image)

    def set_imagepath(self, image_path, image_prefix):
        """
        Sets the directory to save the images in the callback_userdatas
        and sets the name prefix from configuration file
        """
        self.callback_data.process_params['directory'] = image_path
        self.callback_data.process_params['filename_prefix'] = image_prefix

    def start_camera(self, live_mode):
        """
        Starts the camera with live mode (True) or starts the camera without
        live mode (False)

        Parameters:
        ----------
        live_mode : string
        """
        if live_mode == 'True':
            if self.fail('StartLive', self.ENABLE):
                self.camera_lib.StopLive()
                raise IOError("Unable to start camera, replug camera"
                              " to the computer")
        else:
            directory = self.callback_data.process_params['directory']
            if(os.path.exists(directory) == 0):
                os.mkdir(directory)
                print("New directory '%s' was created" % directory)
            if self.fail('StartLive', self.DISABLE):
                self.camera_lib.StopLive()
                raise IOError("Unable to start camera, replug camera"
                              " to the computer")

    def stop_camera(self):
        """
        Stops the camera and resets all the camera properties to their default
        values. If a property has automation, the automatic will be enabled.
        """
        cv2.destroyWindow('Window')
        # if self.camera_lib.ResetProperties() != self.SUCCESS:
        #     self.camera_lib.StopLive()
        #     raise IOError("Unable to reset the camera")
        self.camera_lib.StopLive()  # NO ERROR CHECKING HERE
        # if self.camera_lib.StopLive() != self.SUCCESS:
        #     raise IOError("Unable to stop the camera, replug camera to the "
        #                   "computer")
        if self.fail('SetPropertySwitch',
                     'Partial scan', 'Auto-center', self.ENABLE):
            raise RuntimeError("Unable to set camera back to auto-center")
        if self.fail('enableAutoCameraProperty',
                     self.PROP_CAM_AUTO_EXPOSURE, self.ENABLE):
            raise RuntimeError(
                "Unable to set camera back to autoexposure mode")
        if self.fail('SetPropertySwitch', 'Gain', 'Auto', self.ENABLE):
            raise RuntimeError("Unable to set camera back to autogain mode")
