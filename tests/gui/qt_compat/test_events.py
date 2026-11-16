"""Tests for module events of viperleed.gui.qt_compat."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6

BOTH = (PYQT5, PYQT6)
# Values of QInputDevice.DeviceType, a bit mask
_MOUSE = 0x01
_TOUCHPAD = 0x04


@fixture(name='events')
def fixture_events(import_with_qt):
    """Return the events module of qt_compat."""
    return import_with_qt('viperleed.gui.qt_compat.events')


@fixture(name='press_recorder')
def fixture_press_recorder(events, qtbot, qtw):
    """Return a shown widget that stores where the mouse was pressed."""
    class _PressRecorder(qtw.QWidget):
        """A widget storing the global position of mouse presses."""

        position = None

        # pylint: disable-next=invalid-name
        def mousePressEvent(self, event):
            """Store the global position of `event`."""
            self.position = events.get_global_position(event)

    widget = _PressRecorder()
    qtbot.addWidget(widget)
    widget.resize(100, 100)
    with qtbot.waitExposed(widget):
        widget.show()
    return widget


@parametrize(binding=BOTH)
def test_copy_timestamp(binding, events, mocker):
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


class TestGetGlobalPosition:
    """Tests for the get_global_position function."""

    @parametrize(binding=BOTH)
    def test_calls(self, binding, events, mocker):
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

    def test_real_event(self, press_recorder, qtbot, qtc):
        """Check the global position of a real mouse click."""
        local = qtc.QPoint(30, 40)
        qtbot.mouseClick(press_recorder, qtc.Qt.MouseButton.LeftButton,
                         pos=local)
        assert press_recorder.position == press_recorder.mapToGlobal(local)


class TestIsWheelFromTouchpad:
    """Tests for the is_wheel_from_touchpad function."""

    @parametrize(by_system=(True, False))
    def test_pyqt5(self, by_system, events, mocker):
        """Check that PyQt5 uses the source of the event."""
        mocker.patch.object(events, 'get_binding', return_value=PYQT5)
        qtc = mocker.patch.object(events, 'qtc')
        source = qtc.Qt.MouseEventSource.MouseEventSynthesizedBySystem
        event = mocker.MagicMock()
        event.source.return_value = (source if by_system
                                     else mocker.MagicMock())
        assert events.is_wheel_from_touchpad(event) == by_system

    @parametrize('device,expect', (
        (_TOUCHPAD, True),
        (_MOUSE, False),
        (_TOUCHPAD | _MOUSE, True),
        ))
    def test_pyqt6(self, device, expect, events, mocker):
        """Check that PyQt6 uses the device type of the event."""
        mocker.patch.object(events, 'get_binding', return_value=PYQT6)
        qtg = mocker.patch.object(events, 'qtg')
        qtg.QInputDevice.DeviceType.TouchPad = _TOUCHPAD
        event = mocker.MagicMock()
        event.deviceType.return_value = device
        assert events.is_wheel_from_touchpad(event) == expect

    @pytest.mark.usefixtures('qapp')
    def test_real_wheel_event(self, events, qtc, qtg):
        """Check that a wheel event made by the application works."""
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

    @pytest.mark.usefixtures('qapp')
    def test_real_event(self, events, qtc):
        """Check the attributes of a real event, with QPoint positions."""
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

    @parametrize(binding=BOTH)
    def test_source_only_pyqt5(self, binding, events, mocker):
        """Check that the source is passed only with PyQt5."""
        mocker.patch.object(events, 'get_binding', return_value=binding)
        qtc = mocker.patch.object(events, 'qtc')
        qtg = mocker.patch.object(events, 'qtg')
        source = mocker.MagicMock(name='source')
        events.make_mouse_event('type', 'local', 'window', 'screen',
                                'button', 'buttons', 'modifiers', source)
        point = qtc.QPointF.return_value
        expect_args = ['type', point, point, point,
                       'button', 'buttons', 'modifiers']
        if binding == PYQT5:
            expect_args.append(source)
        qtg.QMouseEvent.assert_called_once_with(*expect_args)
