"""Module camera.drivers.imagingsource of viperleed.???

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-12
Author: Michele Riva
Author: Florian Doerr

This module contains drivers for cameras produced by The Imaging Source.
"""

import os
from viperleed.guilib.measure.camera.drivers.imagingsource.tisgrabber import (
    WindowsCamera, FrameReadyCallbackType, SinkFormat
    )
from viperleed.guilib.measure.camera.drivers.imagingsource.winerrors import (
    ImagingSourceError
    )

if 'nt' in os.name:
    ISCamera = WindowsCamera
else:
    raise EnvironmentError("Imaging Source cameras are unsupported on "
                           "non-Windows platforms.")
