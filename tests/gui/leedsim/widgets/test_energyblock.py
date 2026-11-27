"""Tests for module energyblock of viperleed.gui.leedsim.widgets."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from pytest_cases import fixture


class TestEnergyBlock:
    """Tests for the EnergyBlock class."""

    @fixture(name='block')
    def fixture_block(self, import_with_qt, qtbot):
        """Return an EnergyBlock."""
        module = import_with_qt('viperleed.gui.leedsim.widgets.energyblock')
        block = module.EnergyBlock()
        qtbot.addWidget(block)
        return block

    def test_limits_alignment(self, block, check_aligned, qtc):
        """Check that the energy limits are aligned, not stretched."""
        check_aligned(block.bottomWidget.layout(), block.limits,
                      qtc.Qt.AlignmentFlag.AlignLeft
                      | qtc.Qt.AlignmentFlag.AlignTop)
