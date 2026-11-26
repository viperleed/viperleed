"""Tests for module errorbox of viperleed.gui.leedsim.dialogs."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

from pytest_cases import fixture

_ERROR_MSG = 'a test error'
_ERROR_TEXT = 'Some error'
_EXTRA_TEXT = ' Extra text.'


@fixture(name='errorbox')
def fixture_errorbox(import_with_qt):
    """Return the errorbox module."""
    return import_with_qt('viperleed.gui.leedsim.dialogs.errorbox')


def test_error_box(errorbox, qtbot, mocker):
    """Check that an ErrorBox is shown only when handling exceptions."""
    box = errorbox.ErrorBox(error_while='testing', text=_ERROR_TEXT)
    qtbot.addWidget(box)
    show_modal = mocker.patch(  # Would block the test
        'PyQt5.QtWidgets.QMessageBox.exec_',
        # Like QDialog.exec_, return an int, not the enum member
        return_value=int(box.StandardButton.Ok),
        )
    assert box.exec_() == box.StandardButton.Ok
    show_modal.assert_not_called()
    try:
        _raise_value_error()
    except ValueError:
        assert box.exec_(_EXTRA_TEXT) == box.StandardButton.Ok
    show_modal.assert_called_once()
    assert box.text() == _ERROR_TEXT + _EXTRA_TEXT
    assert _ERROR_MSG in box.detailedText()


def _raise_value_error():
    """Raise a ValueError."""
    raise ValueError(_ERROR_MSG)
