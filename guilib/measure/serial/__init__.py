"""Package serial of viperleed.?????.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-08
Author: Michele Riva
Author: Florian Doerr

This package contains classes to communicate with LEED controllers via serial.
"""

import sys

# Here import all the concrete reimplementations of serialabc.SerialABC
# so that they are available to get_serial()
from viperleed.guilib.measure.serial.viperleedserial import ViPErLEEDSerial


def get_serial(class_name):
    """Return the serial class given its name."""
    try:
        cls = getattr(sys.modules[__name__], class_name)
    except AttributeError as err:
        raise ValueError(
            f"No serial class named {class_name} found."
            ) from err
    return cls
