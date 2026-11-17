"""Tests for module decorators of viperleed.gui.widgets."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

from collections import Counter

from pytest_cases import fixture
from pytest_cases import parametrize

_PRESS = 'press'
_RELEASE = 'release'
_WHEEL = 'wheel'


@fixture(name='decorators')
def fixture_decorators(import_with_qt):
    """Return the decorators module."""
    return import_with_qt('viperleed.gui.widgets.decorators')


@fixture(name='make_counter')
def fixture_make_counter(qtc):
    """Return a class of event filters that count mouse events."""
    kinds = {qtc.QEvent.Type.MouseButtonPress: _PRESS,
             qtc.QEvent.Type.MouseButtonRelease: _RELEASE,
             qtc.QEvent.Type.Wheel: _WHEEL}

    class _Counter(qtc.QObject):
        """An event filter counting the mouse events it sees."""

        def __init__(self):
            super().__init__()
            self.counts = Counter()

        # pylint: disable-next=invalid-name
        def eventFilter(self, _, event):
            """Count `event`, if it is a mouse event."""
            kind = kinds.get(event.type())
            if kind:
                self.counts[kind] += 1
            return False

    return _Counter


class TestBroadcastMouse:
    """Tests for the broadcast_mouse decorator."""

    @fixture(name='window')
    def fixture_window(self, decorators, import_with_qt, qtbot, qtw):
        """Return a shown window broadcasting to two receivers.

        The window contains, in this order, a receiver, a label nested
        in a container (both not accepting mouse events), and another
        receiver.
        """
        spinboxes = import_with_qt('viperleed.gui.widgets.spinboxes')

        @decorators.broadcast_mouse
        class _Window(qtw.QWidget):
            """A window broadcasting mouse events."""

        window = _Window()
        layout = qtw.QVBoxLayout(window)
        window.receivers = (spinboxes.TextBox('1.0'),
                            spinboxes.TextBox('2.0'))
        for receiver in window.receivers:
            receiver.setStep(1, 'add')  # Needed for wheel events
        container = qtw.QWidget()
        window.label = qtw.QLabel('Not a receiver')
        qtw.QVBoxLayout(container).addWidget(window.label)
        layout.addWidget(window.receivers[0])
        layout.addWidget(container)
        layout.addWidget(window.receivers[1])
        qtbot.addWidget(window)
        with qtbot.waitExposed(window):
            window.show()
        return window

    @fixture(name='count_events')
    def fixture_count_events(self, window, make_counter):
        """Return the counters of mouse events of the receivers."""
        counters = [make_counter() for _ in window.receivers]
        for receiver, counter in zip(window.receivers, counters):
            receiver.installEventFilter(counter)
        return counters

    @fixture(name='send_wheel')
    def fixture_send_wheel(self, qtc, qtg, qtw):
        """Return a function that sends a wheel event to a widget.

        The event is created here, so its timestamp is zero.
        """
        def _send(widget):
            pos = qtc.QPoint(1, 1)
            wheel = qtg.QWheelEvent(qtc.QPointF(pos),
                                    qtc.QPointF(widget.mapToGlobal(pos)),
                                    qtc.QPoint(0, 0),
                                    qtc.QPoint(0, 120),
                                    qtc.Qt.MouseButton.NoButton,
                                    qtc.Qt.KeyboardModifier.NoModifier,
                                    qtc.Qt.ScrollPhase.NoScrollPhase,
                                    False)
            assert not wheel.timestamp()
            qtw.QApplication.sendEvent(widget, wheel)
        return _send

    def test_first_event_zero_timestamp(self, window, count_events,
                                        send_wheel):
        """Check that a first event with zero timestamp is broadcast."""
        send_wheel(window)
        expect = {_WHEEL: 1}
        assert [c.counts for c in count_events] == [expect, expect]

    @parametrize(target=('window', 'label'))
    # pylint: disable-next=too-many-arguments  # All fixtures
    def test_once_per_receiver(self, target, window, count_events,
                               send_wheel, qtbot, qtc):
        """Check that each receiver gets each mouse event once.

        Mouse events on the label also reach the container and the
        window, as the label does not accept them: the broadcast
        events must not be broadcast again, and neither must those
        passed on to the parents.
        """
        widget = window if target == 'window' else window.label
        pos = qtc.QPoint(1, 1)
        qtbot.mouseClick(widget, qtc.Qt.MouseButton.LeftButton, pos=pos)
        send_wheel(widget)
        expect = {_PRESS: 1, _RELEASE: 1, _WHEEL: 1}
        assert [c.counts for c in count_events] == [expect, expect]
