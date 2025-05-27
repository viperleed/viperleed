"""Tests for module cli of viperleed.gui."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-27'
__license__ = 'GPLv3+'

import importlib
from pathlib import Path
import sys

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.cli import ViPErLEEDGUICLI
from viperleed.gui.cli import commandline_main
from viperleed.gui.cli import is_commandline_mode
from viperleed.gui.cli import gui_main

_MODULE = 'viperleed.gui.cli'


class TestViPErLEEDGUICLI:
    """Tests for the ViPErLEEDGUICLI class."""

    _valid_args = (
        (['--nogui'], {'nogui': True}),
        )

    @parametrize('args,attrs', _valid_args)
    def test_parse_args_valid(self, args, attrs):
        """Check result of successfully parsing command-line arguments."""
        cli = ViPErLEEDGUICLI()
        parsed = cli.parse_cli_args(args)
        for attr_name, attr_value in attrs.items():
            assert getattr(parsed, attr_name) == attr_value

    def test_call_commandline_mode(self, mocker):
        """Check inner calls when running in command-line mode."""
        cli = ViPErLEEDGUICLI()
        mocker.patch(f'{_MODULE}.is_commandline_mode', return_value=True)
        mock_cmd = mocker.patch(f'{_MODULE}.commandline_main')
        mocker.patch.object(cli, 'parse_cli_args')

        result = cli([])
        assert result is mock_cmd.return_value
        mock_cmd.assert_called_once_with()

    def test_call_gui_mode(self, mocker):
        """Check inner calls when running in GUI mode."""
        cli = ViPErLEEDGUICLI()
        mocker.patch(f'{_MODULE}.is_commandline_mode', return_value=False)
        mock_cmd = mocker.patch(f'{_MODULE}.gui_main')
        mocker.patch.object(cli, 'parse_cli_args')

        result = cli([])
        assert result is mock_cmd.return_value
        mock_cmd.assert_called_once_with()


class TestCommandlineMain:
    """Tests for the commandline_main function."""

    _print = {
        'no pyqt': (
            False,
            'Running in command line version because PyQt5 was not found...'
            'Not implemented yet.\n'
            ),
        'with pyqt': (
            True,
            'Running in command line version...'
            'Not implemented yet.\n'
            ),
        }

    @parametrize('has_qt,expect_print', _print.values(), ids=_print)
    def test_print(self, has_qt, expect_print, mocker, capsys):
        """Check expected writing to standard output."""
        mocker.patch(f'{_MODULE}.has_pyqt', return_value=has_qt)
        error = commandline_main()
        assert not error
        stdout = capsys.readouterr().out
        assert stdout == expect_print


class TestIsCommandlineMode:
    """Tests for the is_commandline_mode function."""

    _cmd_line_info = {
        'normal gui': (
            False,
            {'has_pyqt': True,
             'has_graphics': True,
             'nogui': False},
             ),
        'nogui flag': (
            True,
            {'has_pyqt': True,
             'has_graphics': True,
             'nogui': True},
            ),
        'no qt': (
            True,
            {'has_pyqt': False,
             'has_graphics': True,
             'nogui': False},
            ),
        'no graphics': (
            True,
            {'has_pyqt': True,
             'has_graphics': False,
             'nogui': False},
            ),
        }

    @parametrize('expect,mocks', _cmd_line_info.values(), ids=_cmd_line_info)
    def test_result(self, expect, mocks, mocker):
        """Check the expected return value."""
        for mock_name in ('has_pyqt', 'has_graphics'):
            mocker.patch(f'{_MODULE}.{mock_name}',
                         return_value=mocks[mock_name])
        args = mocker.MagicMock(nogui=mocks['nogui'])
        result = is_commandline_mode(args)
        assert result == expect


class TestGuiMain:
    """Tests for the gui_main function."""

    @fixture(name='mocks')
    def fixture_mock_gui_components(self, mocker):
        mocks = {
            'catch_crash': mocker.patch(f'{_MODULE}.catch_gui_crash'),
            'font_path': mocker.patch(f'{_MODULE}.resources_path',
                                      return_value='fake/path/to/fonts'),
            'qtc': mocker.patch(f'{_MODULE}.qtc'),
            'qtg': mocker.patch(f'{_MODULE}.qtg'),
            'qtw': mocker.patch(f'{_MODULE}.qtw'),
            'select_cls': mocker.patch(f'{_MODULE}.ViPErLEEDSelectPlugin'),
            }
        mocks['app'] = mocks['qtw'].QApplication.return_value
        mocks['select'] = mocks['select_cls'].return_value
        return mocks

    def test_execution(self, mocks, capsys):
        """Check the result of calling gui_main."""
        qapp = mocks['app']
        selector = mocks['select']
        result = gui_main()
        assert result is qapp.exec_.return_value

        # Check important calls
        qapp.exec_.assert_called_once_with()
        selector.show.assert_called_once_with()

        # Check also printing
        stdout = capsys.readouterr().out
        assert stdout.rstrip() == 'Loading GUI...Done'

    def test_fonts_loaded(self, mocks, mocker):
        """Ensure the font database is updated."""
        gui_main()
        mocks['font_path'].assert_called_once_with('gui/fonts')
        base_path = Path(mocks['font_path'].return_value)
        font_db = mocks['qtg'].QFontDatabase.addApplicationFont
        assert font_db.mock_calls == [
            mocker.call(str(base_path/'DejaVuSans.ttf')),
            mocker.call(str(base_path/'cmunrm.otf')),
            ]

    def test_high_dpi_set(self, mocks, mocker):
        """Ensure that the application is set to support high-DPI screens."""
        gui_main()
        qtc = mocks['qtc']
        qapp = mocks['qtg'].QGuiApplication
        set_high_dpi = qapp.setAttribute
        assert set_high_dpi.mock_calls == [
            mocker.call(qtc.Qt.AA_EnableHighDpiScaling),
            mocker.call(qtc.Qt.AA_UseHighDpiPixmaps),
            ]

    def test_icon_set(self, mocks):
        """Ensure that the application icon was set."""
        gui_main()
        qapp = mocks['app']
        qapp.setWindowIcon.assert_called_once()

    def test_implementation(self, mocks):
        """Check inner implementation details."""
        gui_main()
        mocks['catch_crash'].assert_called_once_with()
        mocks['select_cls'].assert_called_once_with()
        mocks['qtw'].QApplication.assert_called_once_with(sys.argv)
