"""Tests for module events of viperleed.gui.qt_compat."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

from pytest_cases import parametrize

from viperleed.gui.qt_compat import events
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6


@parametrize(binding=(PYQT5, PYQT6))
def test_get_global_position_calls(binding, mocker):
    """Check that the method of the selected binding is used."""
    mocker.patch.object(events, 'get_binding', return_value=binding)
    event = mocker.MagicMock()
    position = events.get_global_position(event)
    if binding == PYQT5:
        assert position is event.globalPos.return_value
        event.globalPosition.assert_not_called()
    else:
        to_point = event.globalPosition.return_value.toPoint
        assert position is to_point.return_value
        event.globalPos.assert_not_called()


def test_get_global_position_real_event(qtbot):
    """Check the global position of a real mouse click."""
    # Import here, as qtbot skips the test without Qt
    # pylint: disable=import-outside-toplevel
    from viperleed.gui.qt_compat import QtCore as qtc  # noqa: PLC0415
    from viperleed.gui.qt_compat import QtWidgets as qtw  # noqa: PLC0415

    class _Widget(qtw.QWidget):  # pylint: disable=too-few-public-methods
        """A widget that stores the global position of mouse presses."""

        position = None

        def mousePressEvent(self, event):  # pylint: disable=invalid-name
            """Store the global position of `event`."""
            self.position = events.get_global_position(event)

    widget = _Widget()
    qtbot.addWidget(widget)
    widget.resize(100, 100)
    with qtbot.waitExposed(widget):
        widget.show()
    local = qtc.QPoint(30, 40)
    qtbot.mouseClick(widget, qtc.Qt.MouseButton.LeftButton, pos=local)
    assert widget.position == widget.mapToGlobal(local)
