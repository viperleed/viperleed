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

if sys.platform.startswith('win'):
    from viperleed.guilib.measure.camera.drivers.imagingsource.tisgrabber import (
        WindowsCamera as ISCamera, FrameReadyCallbackType, SinkFormat
        )
    from viperleed.guilib.measure.camera.drivers.imagingsource.winerrors import (
        ImagingSourceError
        )
else:
    raise EnvironmentError("Imaging Source cameras are currently "
                           "unsupported on non-Windows platforms.")
