"""Package imagingsource of viperleed.gui.camera.drivers.

Defines the driver for cameras produced by The Imaging Source.

Modules
-------
models
    Collects information about known camera models.
properties
    Functionality for accessing settings of cameras.
tisgrabber
    The actual driver for low-level interaction with the hardware.
winerrors
    Functionality for reporting Python exceptions when the low-level
    driver code returns error codes.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-12'
__license__ = 'GPLv3+'

import sys

if not sys.platform.startswith('win'):
    # TODO: it would be better to have a custom exception for this,
    # but we should first remove the various __init__ imports in
    # gui. When we do have a custom exception, we should also fix
    # the except block in camera.__init__.
    raise ImportError('Imaging Source cameras are currently '
                      'unsupported on non-Windows platforms.')

from viperleed.gui.measure.camera.drivers.imagingsource.tisgrabber import (
    FrameReadyCallbackType,
    SinkFormat,
    WindowsCamera as ISCamera,
    )
from viperleed.gui.measure.camera.drivers.imagingsource.winerrors import (
    ImagingSourceError,
    )
