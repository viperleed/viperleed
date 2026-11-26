"""Tests for module mainwindow of viperleed.gui.leedsim."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

from pytest_cases import fixture

from viperleed.gui.leedsim.classes.leedparameters import (
    LEEDParametersList,
    )


class TestLEEDPatternSimulator:
    """Tests for the main paths of the (new) LEED-pattern simulator."""

    @fixture(name='default_parameters')
    def fixture_default_parameters(self, import_with_qt):
        """Return the LEED parameters used when starting anew."""
        dialog = import_with_qt('viperleed.gui.leedsim.dialogs.newfiledialog')
        return LEEDParametersList(dialog.DEFAULT_STARTUP)

    @fixture(name='simulator')
    def fixture_simulator(self, import_with_qt, no_blocking_dialogs,
                          show, qtw):
        """Yield a shown LEEDPatternSimulator, then close it."""
        module = import_with_qt('viperleed.gui.leedsim.mainwindow')
        simulator = module.LEEDPatternSimulator()
        show(simulator)
        yield simulator
        # Do not save changes when closing
        no_blocking_dialogs.return_value = qtw.QMessageBox.StandardButton.No
        assert simulator.close()

    def test_new_file_dialog_done(self, simulator, qtbot, qtc):
        """Check that the new-file dialog can be opened and accepted."""
        # pylint: disable=protected-access
        simulator._on_file_new_or_edit_pressed()
        dialog = simulator._dialogs['file_new']
        qtbot.waitExposed(dialog)
        with qtbot.waitSignal(dialog.editing_finished):
            qtbot.mouseClick(dialog._ctrls['done'],
                             qtc.Qt.MouseButton.LeftButton)
        assert not dialog.isVisible()

    def test_new_pattern(self, simulator, default_parameters, qtbot):
        """Check that a LEED pattern can be created from parameters."""
        # pylint: disable=protected-access
        dialog = simulator._dialogs['file_new']
        with qtbot.waitSignal(dialog.editing_finished):
            dialog.editing_finished.emit(default_parameters)
        assert not simulator.saved
