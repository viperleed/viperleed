"""Package camera of viperleed.gui.measure.

Defines classes and functions to handle cameras.
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
