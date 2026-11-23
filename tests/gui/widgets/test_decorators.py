"""Tests for module decorators of viperleed.gui.widgets."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'

from PyQt5 import QtWidgets as qtw

from viperleed.gui.widgets.decorators import broadcast_mouse


class TestBroadcastMouse:
    """Tests for the broadcast_mouse class decorator."""

    def test_decorates_once(self):
        """Check that the class is modified when decorating, only.

        Modifying it again at each instantiation frees, with PyQt5,
        the slots that the children of existing instances use.
        """
        class _Window(qtw.QWidget):
            """A window that broadcasts mouse events."""

        original = dict(vars(_Window))
        decorated = broadcast_mouse(_Window)
        # No wrapper, so that instantiating calls the class itself
        assert decorated is _Window
        patched = dict(vars(_Window))
        assert patched != original
        # pylint: disable-next=magic-value-comparison
        assert 'on_child_destroyed' in patched
