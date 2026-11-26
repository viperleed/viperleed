"""Tests for module mainwindow_old of viperleed.gui.leedsim."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from pytest_cases import fixture

_MODULE = 'viperleed.gui.leedsim.mainwindow_old'


class TestLEEDPatternSimulator:
    """Tests for the old LEED-pattern simulator."""

    @fixture(name='layout_simulator')
    def fixture_layout_simulator(self, import_with_qt, show, mocker,
                                 no_blocking_dialogs):
        """Yield a simulator in 'layout' configuration, then close it.

        This configuration is unused: initUI always sets testConfig
        to 'fix'. Ignore this, to test the 'layout' code paths.
        """
        cls = import_with_qt(_MODULE).LEEDPatternSimulator
        mocker.patch.object(cls,
                            'testConfig',
                            property(lambda _: 'layout', lambda *_: None),
                            create=True)
        simulator = cls()
        show(simulator)
        yield simulator
        no_blocking_dialogs.assert_not_called()
        closed = simulator.close()  # Answers Yes to "Are you sure?"
        assert closed

    @fixture(name='simulator')
    def fixture_simulator(self, import_with_qt, no_blocking_dialogs, show):
        """Yield a shown LEEDPatternSimulator, then close it."""
        simulator = import_with_qt(_MODULE).LEEDPatternSimulator()
        show(simulator)
        yield simulator
        no_blocking_dialogs.assert_not_called()
        closed = simulator.close()  # Answers Yes to "Are you sure?"
        assert closed

    def test_cell_shapes_alignment(self, layout_simulator, check_aligned,
                                   qtc, qtw):
        """Check that the cell shapes are aligned, not stretched."""
        cell_shapes = layout_simulator.cellShapes
        grid = layout_simulator.centralWidget().layout()
        layout, = (lay for lay in grid.findChildren(qtw.QBoxLayout)
                   if lay.indexOf(cell_shapes) >= 0)
        check_aligned(layout, cell_shapes,
                      qtc.Qt.AlignmentFlag.AlignHCenter
                      | qtc.Qt.AlignmentFlag.AlignTop)

    def test_move_event(self, simulator, check_annotations_follow_move):
        """Check that the hover annotations follow the moved window."""
        check_annotations_follow_move(simulator, simulator.recSpace.annots)
