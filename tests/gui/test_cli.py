"""Tests for module cli of viperleed.gui."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-27'
__license__ = 'GPLv3+'

from contextlib import nullcontext
from pathlib import Path
import sys

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.cli import ViPErLEEDGUICLI
from viperleed.gui.cli import commandline_main
from viperleed.gui.cli import gui_main
from viperleed.gui.cli import import_graphics_modules

_MODULE = 'viperleed.gui.cli'


class TestCheckCanRunGui:
    """Tests for the check_can_run_gui method of ViPErLEEDGUICLI."""

    _err_msg = ('Cannot execute the ViPErLEED graphical user '
                'interface because ')
    _install = {
        'not found': '',
        'found': 'sudo apt install lib_one lib_two',
        }

    @parametrize(install=_install.values(), ids=_install)
    def test_missing_dependencies(self, install, mocker, capsys):
        """Check complaints when dependencies are missing."""
        mocker.patch(f'{_MODULE}.has_pyqt', return_value=True)
        missing_dict = {'base_one': ('dep_one', 'dep_two'),
                        'base_two': ('dep_three', 'dep_four')}
        mocker.patch(f'{_MODULE}.find_missing_qt_dependencies',
                     return_value=missing_dict)
        mock_find_install = mocker.patch(
            f'{_MODULE}.Qt5DependencyFinder.find_install_for_libs',
            return_value=install
            )

        cli = ViPErLEEDGUICLI()
        with pytest.raises(SystemExit):
            cli.check_can_run_gui()
        missing = ('dep_one', 'dep_two', 'dep_three', 'dep_four')
        mock_find_install.assert_called_once_with(missing)

        expect_err = self._err_msg + '''\
the following PyQt5 dependencies are missing on your system:
    dep_one
    dep_two
    dep_three
    dep_four
Try again after installing them'''
        expect_err += ('.' if not install else ''' via
    sudo apt install lib_one lib_two''')
        stderr = capsys.readouterr().err.strip()
        assert stderr.endswith(expect_err)

    _can_run = {
        'normal gui': (
            {'has_pyqt': True,
             'has_graphics': True,
             'find_missing_qt_dependencies': {}},
             '',
             ),
        'no qt': (
            {'has_pyqt': False,
             'has_graphics': True,
             'find_missing_qt_dependencies': {}},
            'PyQt5 is not installed.',
            ),
        'no graphics': (
            {'has_pyqt': True,
             'has_graphics': False,
             'find_missing_qt_dependencies': {}},
            ('the system appears to have no graphics capability (i.e., '
             'no monitor was detected). Try once more if this is the first '
             'time you execute the GUI, or have just updated viperleed.'),
            ),
        }

    @parametrize('mocks,err_msg', _can_run.values(), ids=_can_run)
    def test_result(self, mocks, err_msg, mocker, capsys):
        """Check the expected outcome."""
        for mock_name in ('has_pyqt', 'has_graphics'):
            mocker.patch(f'{_MODULE}.{mock_name}',
                         return_value=mocks[mock_name])
        cli = ViPErLEEDGUICLI()
        context = pytest.raises(SystemExit) if err_msg else nullcontext()
        with context:
            cli.check_can_run_gui()
        stderr = capsys.readouterr().err.strip()
        if not err_msg:
            assert not stderr
        else:
            expect = self._err_msg + err_msg
            assert stderr.endswith(expect)


# pylint: disable-next=too-few-public-methods  # It's parameterized
class TestCommandlineMain:
    """Tests for the commandline_main function."""

    _print = {
        'no pyqt': (
            False,
            'Running in command line version...'
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


class TestGuiMain:
    """Tests for the gui_main function."""

    @fixture(name='mocks')
    def fixture_mock_gui_components(self, mocker):
        """Replace implementation details with mocks."""
        mock_imports = {
            'qtc': mocker.MagicMock(),
            'qtg': mocker.MagicMock(),
            'qtw': mocker.MagicMock(),
            'widgets_lib': mocker.MagicMock(),
            'select_cls': mocker.MagicMock(name='ViPErLEEDSelectPlugin'),
            }
        mocker.patch(f'{_MODULE}.import_graphics_modules',
                     return_value=mock_imports.values())
        return {
            'font_path': mocker.patch(f'{_MODULE}.resources_path',
                                      return_value='fake/path/to/fonts'),
            **mock_imports,
            'catch_crash': mock_imports['widgets_lib'].catch_gui_crash,
            'raise_qt': mock_imports['widgets_lib'].raise_on_qt_messages,
            'select': mock_imports['select_cls'].return_value,
            'app': mock_imports['qtw'].QApplication.return_value,
            'suppress_warnings': mocker.patch(
                f'{_MODULE}.suppress_file_permission_warnings'
                ),
            }

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
        success_load_prints = 'Loading GUI...Done'
        stdout = capsys.readouterr().out
        assert stdout.rstrip() == success_load_prints

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
        mocks['catch_crash'].assert_called_once()
        mocks['raise_qt'].assert_called_once_with()
        mocks['select_cls'].assert_called_once_with()
        mocks['suppress_warnings'].assert_called_once_with()
        mocks['qtw'].QApplication.assert_called_once_with(sys.argv)


class TestImportGrpahicsModules:
    """Tests for the import_graphics_modules function."""

    def test_raises(self, mocker):
        """Check complaints without PyQt."""
        mocker.patch(f'{_MODULE}.has_pyqt', return_value=False)
        with pytest.raises(ImportError):
            import_graphics_modules()

    def test_success(self, mocker):
        """Check expected result of importing modules."""
        mock_import = mocker.patch(f'{_MODULE}.import_module')
        result = import_graphics_modules()
        imported_modules = (
            'viperleed.gui.selectplugin',
            'PyQt5.QtCore',
            'PyQt5.QtGui',
            'PyQt5.QtWidgets',
            'viperleed.gui.widgets.lib',
            )
        assert mock_import.mock_calls == [mocker.call(m)
                                          for m in imported_modules]
        n_return_values = 5
        assert len(result) == n_return_values


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
        mock_cmd = mocker.patch(f'{_MODULE}.commandline_main')
        mock_parse = mocker.patch.object(cli, 'parse_cli_args')
        mock_parse.return_value.nogui = True

        result = cli([])
        assert result is mock_cmd.return_value
        mock_cmd.assert_called_once_with()

    def test_call_gui_mode(self, mocker):
        """Check inner calls when running in GUI mode."""
        cli = ViPErLEEDGUICLI()
        mock_cmd = mocker.patch(f'{_MODULE}.gui_main')
        mock_check_gui = mocker.patch.object(cli, 'check_can_run_gui')
        mock_parse = mocker.patch.object(cli, 'parse_cli_args')
        mock_parse.return_value.nogui = False

        result = cli([])
        assert result is mock_cmd.return_value
        mock_cmd.assert_called_once_with()
        mock_check_gui.assert_called_once_with()
