"""Tests for module matricespopup of viperleed.gui.leedsim.widgets."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

import numpy as np
from pytest_cases import fixture


@fixture(name='popup')
def fixture_popup(import_with_qt, show):
    """Return a shown MatricesPopup."""
    module = import_with_qt('viperleed.gui.leedsim.widgets.matricespopup')
    popup = module.MatricesPopup([np.eye(2)], None)
    show(popup)
    return popup


def test_drag(popup, qtbot, qtc, qtg):
    """Check that the pop-up can be dragged."""
    start = popup.pos()
    left = qtc.Qt.MouseButton.LeftButton
    qtbot.mousePress(popup, left, pos=qtc.QPoint(5, 5))
    # QTest.mouseMove does not report pressed buttons: send a move
    # event with the left button pressed via the event handler itself
    local_pos = qtc.QPoint(25, 35)
    move = qtg.QMouseEvent(qtc.QEvent.Type.MouseMove,
                           qtc.QPointF(local_pos),
                           qtc.QPointF(popup.mapToGlobal(local_pos)),
                           left, left, qtc.Qt.KeyboardModifier.NoModifier)
    popup.mouseMoveEvent(move)
    assert popup.pos() - start == qtc.QPoint(20, 30)
