"""Module camera.drivers.imagingsource of viperleed.guilib.measure.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-12
Author: Michele Riva
Author: Florian Doerr

This module contains drivers for cameras produced by The Imaging Source.
"""

import sys

if not sys.platform.startswith('win'):
    # TODO: it would be better to have a custom exception for this,
    # but we should first remove the various __init__ imports in
    # guilib. When we do have a custom exception, we should also fix
    # the except block in camera.__init__.
    raise ImportError('Imaging Source cameras are currently '
                      'unsupported on non-Windows platforms.')

from viperleed.guilib.measure.camera.drivers.imagingsource.tisgrabber import (
    FrameReadyCallbackType,
    SinkFormat,
    WindowsCamera as ISCamera,
    )
from viperleed.guilib.measure.camera.drivers.imagingsource.winerrors import (
    ImagingSourceError,
    )
