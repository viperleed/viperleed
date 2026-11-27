"""Tests for module hoverannot of viperleed.gui.leedsim.widgets."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from pytest_cases import fixture


class TestHoverAnnot:
    """Tests for the HoverAnnot class."""

    @fixture(name='annotation')
    def fixture_annotation(self, import_with_qt, qtbot):
        """Return a HoverAnnot with its arrow at non-integer points."""
        module = import_with_qt('viperleed.gui.leedsim.widgets.hoverannot')
        annotation = module.HoverAnnot(center=(50.2, 50.7))
        annotation.head_pos = (60.4, 60.6)
        qtbot.addWidget(annotation)
        return annotation

    def test_size_hint(self, annotation):
        """Check that the size hint fits the frame of the text."""
        size = annotation.sizeHint()
        frame = annotation.frame.boundingRect()
        assert size.width() >= frame.width()
        assert size.height() >= frame.height()

    def test_update_position(self, annotation, qtc):
        """Check that the tip of the arrow moves to the arrowhead."""
        annotation.update_position()
        tip = annotation.mapToGlobal(annotation.tip_pos.toPoint())
        assert tip == qtc.QPoint(60, 61)
