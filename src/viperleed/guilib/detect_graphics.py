"""Module detect_graphics of viperleed.guilib.

Defines functionality to detect whether the current system is capable
of displaying widgets.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-22'
__license__ = 'GPLv3+'

from multiprocessing import Process
import sys

try:
    from PyQt5.QtWidgets import QApplication
except ImportError:
    _HAS_PYQT = False
else:
    _HAS_PYQT = True

# How long (in seconds) to wait before deciding there's no graphics?
_TIMEOUT = 1


def _init_dummy_qapp():
    """Initialize a dummy QApplication instance.

    This function is only meant to be used for checking whether the
    current platform has graphical capabilities. Calling this will
    block forever on a system that has no graphical capability.

    Returns
    -------
    None.
    """
    _ = QApplication([])
    sys.exit(0)


def has_graphics():
    """Return whether the system has graphical capabilities.

    This function will block for one second on systems that do not
    have graphical capabilities.

    Returns
    -------
    bool
    """
    if not _HAS_PYQT:
        return False
    proc = Process(target=_init_dummy_qapp)
    proc.start()
    proc.join(_TIMEOUT)
    if proc.is_alive():
        # The QApplication is still trying to start up. This means we
        # have no display.
        proc.terminate()
        return False
    return True


def has_pyqt():
    """Return whether PyQt5 is installed."""
    return _HAS_PYQT
