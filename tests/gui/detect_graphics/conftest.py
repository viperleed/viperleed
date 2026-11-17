"""Shared definitions for the tests of module detect_graphics."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'

from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6

BOTH = (PYQT6, PYQT5)  # In the default order of preference
MODULE = 'viperleed.gui.detect_graphics'
