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

from viperleed.gui.qt_compat import QtCore as qtc
from viperleed.gui.qt_compat import QtGui as qtg
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import get_binding


def get_global_position(event):
    """Return the position of a mouse `event` in screen coordinates.

    PyQt5 only has QMouseEvent.globalPos (a QPoint), PyQt6 only has
    QMouseEvent.globalPosition (a QPointF). Qt 6 still has globalPos,
    deprecated, but PyQt6 does not expose it.

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


def is_wheel_from_touchpad(event):
    """Return whether a wheel `event` comes from a touchpad.

    PyQt5 only has QWheelEvent.source, which is
    MouseEventSynthesizedBySystem for touchpads. PyQt6 only has
    QInputEvent.deviceType, a QInputDevice.DeviceType flag that
    includes TouchPad for touchpads.

    Parameters
    ----------
    event : QWheelEvent
        The event to be checked.

    Returns
    -------
    from_touchpad : bool
        Whether `event` comes from a touchpad, rather than from
        a mouse wheel.
    """
    if get_binding() == PYQT5:
        by_system = qtc.Qt.MouseEventSource.MouseEventSynthesizedBySystem
        return event.source() == by_system
    touchpad = qtg.QInputDevice.DeviceType.TouchPad
    return bool(event.deviceType() & touchpad)


# This mirrors the signature of the QMouseEvent constructor
# pylint: disable-next=too-many-arguments,too-many-positional-arguments
def make_mouse_event(event_type, local_pos, window_pos, screen_pos,
                     button, buttons, modifiers, source):
    """Return a new QMouseEvent.

    PyQt5 has a QMouseEvent constructor that takes the source of the
    event, PyQt6 does not: there, the source cannot be set, as Qt 6
    derives it from the input device. Positions may be given as QPoint
    or QPointF, as PyQt6 only accepts the latter.

    Parameters
    ----------
    event_type : QEvent.Type
        The type of the event (e.g., MouseButtonPress).
    local_pos, window_pos, screen_pos : QPoint or QPointF
        The position of the event relative to the receiver, to its
        window, and to the screen.
    button : Qt.MouseButton
        The button that caused the event.
    buttons : Qt.MouseButton
        The state of all buttons when the event occurred.
    modifiers : Qt.KeyboardModifier
        The state of the keyboard modifiers when the event occurred.
    source : Qt.MouseEventSource
        The source of the event. Only used with PyQt5.

    Returns
    -------
    event : QMouseEvent
        The new event.
    """
    positions = (qtc.QPointF(pos) for pos in (local_pos, window_pos,
                                              screen_pos))
    args = (event_type, *positions, button, buttons, modifiers)
    if get_binding() == PYQT5:
        return qtg.QMouseEvent(*args, source)
    return qtg.QMouseEvent(*args)
