"""Tests for module newfiledialog_old of viperleed.gui.leedsim.dialogs."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

_MODULE = 'viperleed.gui.leedsim.dialogs.newfiledialog_old'


class TestNewFileDialog:
    """Tests for the NewFileDialog class."""

    @fixture(name='dialog')
    def fixture_dialog(self, import_with_qt, show):
        """Return a shown NewFileDialog."""
        dialog = import_with_qt(_MODULE).NewFileDialog()
        show(dialog)
        return dialog

    @parametrize(key=('Key_Enter', 'Key_Return'))
    @pytest.mark.usefixtures('no_blocking_dialogs')  # Fail, do not hang
    def test_key_clicks_focused_button(self, key, dialog, qtbot, qtc):
        """Check that Enter and Return click the focused Cancel button."""
        button = dialog.cancelBut
        button.setFocus()
        qtbot.waitUntil(button.hasFocus)
        with qtbot.waitSignal(button.clicked, timeout=1000):
            qtbot.keyClick(button, getattr(qtc.Qt.Key, key))
        assert not dialog.isVisible()
