"""Package camera of viperleed.gui.measure.

Defines classes and functions to handle cameras.

Modules
-------
abc
    Abstract base classes that outline the common interface for all
    cameras.
badpixels
    Objects for finding bad-pixel information and reading/writing
    bad-pixel files.
cameracalibration
    Abstract base classes for performing calibration tasks on cameras.
imageprocess
    Functionality for post-processing frames acquired from cameras.
imagingsource
    Concrete camera class for handling cameras produced by The
    Imaging Source GmbH.
tifffile
    Functionality for reading/writing images in TIFF format.

Packages
--------
drivers
    Functionality for low-level interaction with camera hardware.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-12'
__license__ = 'GPLv3+'

try:
    from viperleed.gui.measure.camera.imagingsource import ImagingSourceCamera
except ImportError:
    # Probably we are in the wrong environment
    pass
