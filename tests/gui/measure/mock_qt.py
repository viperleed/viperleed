"""Module mock_qt of viperleed.tests.gui.measure.

Contains fake Qt objects for test modules.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-07-28'
__license__ = 'GPLv3+'


class _FakeSignal:
    """A minimal signal-like object."""

    def __init__(self):
        """Initialize fake signal."""
        self.connected = []
        self.emitted = 0
        self.disconnected = []

    def connect(self, slot):
        """Store connected slots."""
        self.connected.append(slot)

    def emit(self, *args, **kwargs):
        """Count emissions and call connected slots."""
        self.emitted += 1
        for slot in self.connected:
            slot(*args, **kwargs)

    def disconnect(self, slot):
        """Store disconnected slots and remove from connected."""
        self.disconnected.append(slot)
        if slot in self.connected:
            self.connected.remove(slot)
