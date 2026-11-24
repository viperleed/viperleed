"""Tests for module mainwindow_old of viperleed.gui.leedsim."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from pytest_cases import fixture


class TestLEEDPatternSimulator:
    """Tests for the old LEED-pattern simulator."""

    @fixture(name='simulator')
    def fixture_simulator(self, import_with_qt, no_blocking_dialogs, show):
        """Yield a shown LEEDPatternSimulator, then close it."""
        module = import_with_qt('viperleed.gui.leedsim.mainwindow_old')
        simulator = module.LEEDPatternSimulator()
        show(simulator)
        yield simulator
        no_blocking_dialogs.assert_not_called()
        assert simulator.close()  # Answers Yes to "Are you sure?"

    def test_move_event(self, simulator, check_annotations_follow_move):
        """Check that the hover annotations follow the moved window."""
        check_annotations_follow_move(simulator, simulator.recSpace.annots)
