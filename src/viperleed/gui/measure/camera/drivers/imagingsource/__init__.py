"""Package imagingsource of viperleed.gui.camera.drivers.

Defines the driver for cameras produced by The Imaging Source.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-12'
__license__ = 'GPLv3+'

import os
from viperleed.gui.measure.camera.drivers.imagingsource.tisgrabber import (
    FrameReadyCallbackType,
    SinkFormat,
    WindowsCamera,
    )
from viperleed.gui.measure.camera.drivers.imagingsource.winerrors import (
    ImagingSourceError,
    )

if 'nt' in os.name:
    ISCamera = WindowsCamera
else:
    raise EnvironmentError("Imaging Source cameras are unsupported on "
                           "non-Windows platforms.")
