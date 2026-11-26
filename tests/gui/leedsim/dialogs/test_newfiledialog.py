"""Tests for module newfiledialog of viperleed.gui.leedsim.dialogs."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'

from pytest_cases import fixture


class TestNewFileDialog:
    """Tests for the NewFileDialog class."""

    @fixture(name='dialog')
    def fixture_dialog(self, import_with_qt, show):
        """Return a shown NewFileDialog."""
        module = import_with_qt('viperleed.gui.leedsim.dialogs.newfiledialog')
        dialog = module.NewFileDialog()
        show(dialog)
        return dialog

    def test_energy_changed(self, dialog, qtbot, qtc):
        """Check the maximum energy shown when the user edits it."""
        # pylint: disable-next=protected-access
        e_max = dialog._ctrls['e_max']
        e_max.setFocus()
        qtbot.waitUntil(e_max.hasFocus)
        e_max.selectAll()
        qtbot.keyClicks(e_max, '1234')
        qtbot.keyClick(e_max, qtc.Qt.Key.Key_Return)  # editingFinished
        assert e_max.text() == '1234.0 eV'
        assert e_max.maximumWidth() >= 70
