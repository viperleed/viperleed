"""Module events of viperleed.gui.qt_compat.

Functions that hide the differences between PyQt5 and PyQt6 in the
API of Qt events, where no method works with both bindings.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import get_binding


def get_global_position(event):
    """Return the position of a mouse `event` in screen coordinates.

    PyQt5 only has QMouseEvent.globalPos (a QPoint), PyQt6 only has
    QMouseEvent.globalPosition (a QPointF).

    Parameters
    ----------
    event : QMouseEvent
        The event whose position is requested.

    Returns
    -------
    position : QPoint
        The global position of `event`.
    """
    if get_binding() == PYQT5:
        return event.globalPos()
    return event.globalPosition().toPoint()
