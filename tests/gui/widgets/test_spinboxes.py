"""Tests for module spinboxes of viperleed.gui.widgets."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture


@fixture(name='spinboxes')
def fixture_spinboxes(import_with_qt):
    """Return the spinboxes module."""
    return import_with_qt('viperleed.gui.widgets.spinboxes')


class TestTextBox:
    """Tests for the TextBox class."""

    def test_no_arguments(self, spinboxes):
        """Check that a TextBox can be created without arguments."""
        text_box = spinboxes.TextBox()
        assert not text_box.text()
        assert text_box.parent() is None

    # pylint: disable-next=too-many-arguments  # All fixtures
    def test_wheel(self, spinboxes, show, qtbot, qtc, qtg, qtw):
        """Check that the mouse wheel changes the value of a TextBox."""
        text_box = spinboxes.TextBox('1.0')
        text_box.setStep(1, 'add')
        show(text_box)
        text_box.setFocus()
        qtbot.waitUntil(text_box.hasFocus)
        center = qtc.QPointF(text_box.rect().center())
        wheel = qtg.QWheelEvent(center,
                                qtc.QPointF(text_box.mapToGlobal(
                                    center.toPoint()
                                    )),
                                qtc.QPoint(0, 0),
                                qtc.QPoint(0, 160),  # One step with a mouse
                                qtc.Qt.MouseButton.NoButton,
                                qtc.Qt.KeyboardModifier.NoModifier,
                                qtc.Qt.ScrollPhase.NoScrollPhase,
                                False)
        qtw.QApplication.sendEvent(text_box, wheel)
        assert float(text_box.text()) == pytest.approx(2.0)
