"""Module camera of viperleed.guilib.measure.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-12
Author: Michele Riva
Author: Florian Doerr

This module contains classes and functions to handle cameras
"""

try:
    from viperleed.guilib.measure.camera.imagingsource import (
        ImagingSourceCamera
        )
except ImportError:
    # Probably we are in the wrong environment
    pass
