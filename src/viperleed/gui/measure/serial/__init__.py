"""Package serial of viperleed.gui.measure.

This package contains classes to communicate with LEED controllers via
the serial line.

Modules
-------
abc
    Abstract base classes defining the common interface for all
    serial-communication implementations.
viperleedserial
    Classes for interacting with hardware via the ViPErLEED
    communication protocol.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-08'
__license__ = 'GPLv3+'

# Here import all the concrete subclasses of abc.SerialABC
# so that they are available to class_by_name()
from viperleed.gui.measure.serial.viperleedserial import ViPErLEEDSerial
