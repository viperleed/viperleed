"""Test functionality of imaging source camera.

Author: Michele Riva
Created: 2021-07-23
"""
from pathlib import Path
import sys
import unittest

vpr_path = Path(__file__).resolve().parents[2]
if vpr_path not in sys.path:
    sys.path.append(str(vpr_path))

from viperleed.guilib.measure.camera.drivers.imagingsource.tisgrabber import (
    WindowsCamera, DLLReturns, ImagingSourceError, SinkFormat, CameraProperty
    )


class TestDeviceNotopen(unittest.TestCase):
    """Test appropriate behavior when no device is open"""

    def setUp(self):
        """Initialize the camera DLL, but do not open a device."""
        self.camera = WindowsCamera()

    def tearDown(self):
        """Clean up test class."""
        self.camera = None

    def test_n_devices(self):
        """Test that no devices are found."""
        # I'm not sure yet whether I'll make n_devices private
        if hasattr(self.camera, 'n_devices'):
            assert self.camera.n_devices == 0

    def test_imaging_source_error(self):
        """Test ImagingSourceError is raised when accessing multiple."""
        with self.assertRaises(ImagingSourceError) as err:
            rate = self.camera.frame_rate
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.frame_rate = 30
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            info = self.camera.image_info
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            sink = self.camera.sink_format
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.sink_format = SinkFormat.Y16
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            img_shape = self.camera.video_format_shape
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            fmts = self.camera.video_formats
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.set_video_format('_some_invalid_format_zz')
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.open_by_unique_name('__AbcdEfg')
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.snap_live_image()
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.start_live()
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.stop_live()
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.unique_name_from_index(1000)
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.check_available_cam_properties()
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            cam_prop = self.camera.get_camera_property('EXPOSURE')
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            self.camera.set_camera_property('EXPOSURE', 3)
            print(err.exception)
        with self.assertRaises(ImagingSourceError) as err:
            p_min, p_max = self.camera.get_camera_property_range('EXPOSURE')
            print(err.exception)


if __name__ == '__main__':
    unittest.main()