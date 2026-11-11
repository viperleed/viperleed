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

    class _Widget(qtw.QWidget):
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


_TOUCHPAD = 'touchpad'


class TestIsWheelFromTouchpad:
    """Tests for the is_wheel_from_touchpad function."""

    @parametrize(by_system=(True, False))
    def test_pyqt5(self, by_system, mocker):
        """Check that PyQt5 uses the source of the event."""
        mocker.patch.object(events, 'get_binding', return_value=PYQT5)
        qtc = mocker.patch.object(events, '_import_qt').return_value
        source = qtc.Qt.MouseEventSource.MouseEventSynthesizedBySystem
        event = mocker.MagicMock()
        event.source.return_value = (source if by_system
                                     else mocker.MagicMock())
        assert events.is_wheel_from_touchpad(event) == by_system

    @parametrize(device=(_TOUCHPAD, 'mouse', None))
    def test_pyqt6(self, device, mocker):
        """Check that PyQt6 uses the device of the event."""
        mocker.patch.object(events, 'get_binding', return_value=PYQT6)
        qtg = mocker.patch.object(events, '_import_qt').return_value
        touchpad = qtg.QInputDevice.DeviceType.TouchPad
        event = mocker.MagicMock()
        if device is None:
            event.pointingDevice.return_value = None
        else:
            event.pointingDevice.return_value.type.return_value = (
                touchpad if device == _TOUCHPAD else mocker.MagicMock()
                )
        assert events.is_wheel_from_touchpad(event) == (device == _TOUCHPAD)

    def test_real_wheel_event(self, qapp):
        """Check that a wheel event made by the application works."""
        _ = qapp  # Skips without Qt
        # pylint: disable=import-outside-toplevel
        from viperleed.gui.qt_compat import QtCore as qtc  # noqa: PLC0415
        from viperleed.gui.qt_compat import QtGui as qtg  # noqa: PLC0415
        origin = qtc.QPointF(0, 0)
        wheel = qtg.QWheelEvent(origin, origin, qtc.QPoint(0, 0),
                                qtc.QPoint(0, 120),
                                qtc.Qt.MouseButton.NoButton,
                                qtc.Qt.KeyboardModifier.NoModifier,
                                qtc.Qt.ScrollPhase.NoScrollPhase,
                                False)
        assert not events.is_wheel_from_touchpad(wheel)


class TestMakeMouseEvent:
    """Tests for the make_mouse_event function."""

    @parametrize(binding=(PYQT5, PYQT6))
    def test_source_only_pyqt5(self, binding, mocker):
        """Check that the source is passed only with PyQt5."""
        mocker.patch.object(events, 'get_binding', return_value=binding)
        mock_qt = mocker.patch.object(events, '_import_qt').return_value
        source = mocker.MagicMock(name='source')
        events.make_mouse_event('type', 'local', 'window', 'screen',
                                'button', 'buttons', 'modifiers', source)
        point = mock_qt.QPointF.return_value
        expect_args = ['type', point, point, point,
                       'button', 'buttons', 'modifiers']
        if binding == PYQT5:
            expect_args.append(source)
        mock_qt.QMouseEvent.assert_called_once_with(*expect_args)

    def test_real_event(self, qapp):
        """Check the attributes of a real event, with QPoint positions."""
        _ = qapp  # Skips without Qt
        # pylint: disable=import-outside-toplevel
        from viperleed.gui.qt_compat import QtCore as qtc  # noqa: PLC0415
        left = qtc.Qt.MouseButton.LeftButton
        event = events.make_mouse_event(
            qtc.QEvent.Type.MouseButtonPress,
            qtc.QPoint(1, 2), qtc.QPoint(3, 4), qtc.QPoint(5, 6),
            left, left, qtc.Qt.KeyboardModifier.NoModifier,
            qtc.Qt.MouseEventSource.MouseEventSynthesizedByApplication,
            )
        assert event.type() == qtc.QEvent.Type.MouseButtonPress
        assert event.button() == left
        assert events.get_global_position(event) == qtc.QPoint(5, 6)


@parametrize(binding=(PYQT5, PYQT6))
def test_copy_timestamp(binding, mocker):
    """Check that timestamps are copied only with PyQt5."""
    mocker.patch.object(events, 'get_binding', return_value=binding)
    source, target = mocker.MagicMock(), mocker.MagicMock()
    copied = events.copy_timestamp(source, target)
    assert copied == (binding == PYQT5)
    if copied:
        target.setTimestamp.assert_called_once_with(
            source.timestamp.return_value
            )
    else:
        target.setTimestamp.assert_not_called()
