"""Tests for module cli of viperleed.gui."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-05-27'
__license__ = 'GPLv3+'

from contextlib import nullcontext
from pathlib import Path
import sys

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.cli import _SANITY_TO_ERR_MSG
from viperleed.gui.cli import ViPErLEEDGUICLI
from viperleed.gui.cli import commandline_main
from viperleed.gui.cli import gui_main
from viperleed.gui.cli import import_graphics_modules
from viperleed.gui.detect_graphics import BindingCheck
from viperleed.gui.detect_graphics import PyQtSanity
from viperleed.gui.detect_graphics import has_pyqt
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6
from viperleed.gui.qt_compat.binding import QtApiError

_MODULE = 'viperleed.gui.cli'
_QT_COMPAT = 'viperleed.gui.qt_compat'


class TestCheckCanRunGui:
    """Tests for the check_can_run_gui method of ViPErLEEDGUICLI."""

    _err_msg = ('Cannot execute the ViPErLEED graphical user '
                'interface because ')
    _checking = 'Checking graphics capabilities...'
    _done = 'Done\n'
    _install = {
        'not found': '',
        'found': 'sudo apt install lib_one lib_two',
        }
    # Missing dependencies are not looked for with these
    _no_deps_check = {PyQtSanity.OK,
                      PyQtSanity.NOT_FOUND,
                      PyQtSanity.INVALID_QT_API,
                      PyQtSanity.TOO_OLD}

    _mock_abi_msg = '''
If you are executing viperleed in a conda environment, try one of the
following:
    1. Creating a new, clean environment without Qt by calling
       conda create with the --no-default-packages flag, then
            pip install "viperleed[GUI]"
       there.
    2. Deactivating the current environment first, if you have
       installed viperleed globally.
If none of the above works, or you're not in a conda environment, please
open an issue under https://github.com/viperleed/viperleed/issues.'''
    _mock_skip_msg = '''
If your system has a display, you can start the GUI anyway with
    viperleed gui --skip-graphics-check'''
    _can_run = {
        'normal gui': (PyQtSanity.OK, ''),
        'no qt': (PyQtSanity.NOT_FOUND,
                  'no Qt binding is installed (looked for PyQt6 or PyQt5).'),
        'no graphics': (
            PyQtSanity.NO_DISPLAY,
            ('the system appears to have no graphics capability (i.e., '
             'no monitor was detected).' + _mock_skip_msg),
            ),
        'unresponsive': (
            PyQtSanity.UNRESPONSIVE,
            ('loading PyQt6 did not complete in time. Try once more, or '
             'check whether your system is very busy.' + _mock_skip_msg),
            ),
        'import fail': (
            PyQtSanity.IMPORT_ERROR,
            'PyQt6 was found, but it could not be loaded.' + _mock_abi_msg,
            ),
        'runtime crash': (
            PyQtSanity.RUNTIME_CRASH,
            ('PyQt6 was found, but the GUI crashed at startup.'
             + _mock_abi_msg),
            ),
        }
    # Results of the checks of two bindings that both failed
    _both_failed = (BindingCheck(PYQT6, PyQtSanity.RUNTIME_CRASH),
                    BindingCheck(PYQT5, PyQtSanity.IMPORT_ERROR))
    _also_tried = '\nPyQt5 was also tried, but it could not be loaded.'

    @fixture(autouse=True)
    def fake_candidates(self, mocker):
        """Look for both bindings, independently of the environment."""
        mocker.patch(f'{_MODULE}.get_candidates',
                     return_value=(PYQT6, PYQT5))

    @fixture(name='mock_sanity')
    def fixture_mock_sanity(self, mocker):
        """Return a function that fakes the outcome of the Qt check.

        The function takes the status, the binding it is about, and
        the results of all bindings checked (by default, only that of
        the binding).
        """
        def _mock(sanity, binding=PYQT6, results=None):
            if results is None:
                results = (BindingCheck(binding, sanity),)
            return mocker.patch(f'{_MODULE}.check_pyqt_sanity',
                                return_value=sanity,
                                binding=binding,
                                results=results)
        return _mock

    def test_also_tried(self, mock_sanity, mocker, capsys):
        """Check that other bindings that failed are mentioned."""
        mock_sanity(PyQtSanity.RUNTIME_CRASH, results=self._both_failed)
        mocker.patch(f'{_MODULE}.find_missing_qt_dependencies',
                     return_value={})
        with pytest.raises(SystemExit):
            ViPErLEEDGUICLI().check_can_run_gui()
        stderr = capsys.readouterr().err.strip()
        expect = ('PyQt6 was found, but the GUI crashed at startup.'
                  + self._also_tried + self._mock_abi_msg)
        assert stderr.endswith(expect)

    def test_also_tried_dependencies(self, mock_sanity, mocker, capsys):
        """Check that other bindings follow missing libraries."""
        mock_sanity(PyQtSanity.RUNTIME_CRASH, results=self._both_failed)
        mock_find_deps = mocker.patch(
            f'{_MODULE}.find_missing_qt_dependencies',
            return_value={'libqxcb.so': ('libxcb-cursor',)},
            )
        mocker.patch(f'{_MODULE}.QtDependencyFinder.find_install_for_libs',
                     return_value='')
        with pytest.raises(SystemExit):
            ViPErLEEDGUICLI().check_can_run_gui()
        stderr = capsys.readouterr().err.strip()
        expect = ('the following Qt dependencies of PyQt6 are missing on '
                  'your system:\n    libxcb-cursor\n'
                  'Try again after installing them.' + self._also_tried)
        assert stderr.endswith(expect)
        mock_find_deps.assert_called_once_with()  # Only for PyQt6

    def test_also_tried_not_installed(self, mock_sanity, mocker, capsys):
        """Check that bindings not installed are not mentioned."""
        mock_sanity(PyQtSanity.IMPORT_ERROR,
                    binding=PYQT5,
                    results=(BindingCheck(PYQT6, PyQtSanity.NOT_FOUND),
                             BindingCheck(PYQT5, PyQtSanity.IMPORT_ERROR)))
        mocker.patch(f'{_MODULE}.find_missing_qt_dependencies',
                     return_value={})
        with pytest.raises(SystemExit):
            ViPErLEEDGUICLI().check_can_run_gui()
        stderr = capsys.readouterr().err.strip()
        assert PYQT6 not in stderr  # Not installed
        assert stderr.endswith('PyQt5 was found, but it could not be loaded.'
                               + self._mock_abi_msg)

    def test_invalid_qt_api(self, mocker, capsys):
        """Check complaints when QT_API prevents selecting a binding."""
        mocker.patch(f'{_MODULE}.get_candidates',
                     side_effect=QtApiError('Some reason.'))
        mocker.patch(f'{_MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.INVALID_QT_API)
        find_deps = mocker.patch(f'{_MODULE}.find_missing_qt_dependencies')
        with pytest.raises(SystemExit):
            ViPErLEEDGUICLI().check_can_run_gui()
        captured = capsys.readouterr()
        assert captured.out == self._checking + self._done
        expect = 'no Qt binding can be selected. Some reason.'
        assert captured.err.strip().endswith(self._err_msg + expect)
        find_deps.assert_not_called()

    def test_message_before_check(self, mocker, capsys):
        """Ensure users are told about the check before it starts."""
        def _check_sanity():
            # Only the start of the line is out when checking
            assert capsys.readouterr().out == self._checking
            return PyQtSanity.OK
        mocker.patch(f'{_MODULE}.check_pyqt_sanity', side_effect=_check_sanity)
        mocker.patch(f'{_MODULE}.find_missing_qt_dependencies',
                     return_value={})
        ViPErLEEDGUICLI().check_can_run_gui()
        assert capsys.readouterr().out == self._done

    @parametrize(install=_install.values(), ids=_install)
    @parametrize(sanity=PyQtSanity)
    def test_missing_dependencies(self, install, sanity, mock_sanity,
                                  mocker, capsys):
        """Check complaints when dependencies are missing."""
        mock_sanity(sanity)
        missing_dict = {'base_one': ('dep_one', 'dep_two'),
                        'base_two': ('dep_three', 'dep_four')}
        mock_find_deps = mocker.patch(
            f'{_MODULE}.find_missing_qt_dependencies',
            return_value=missing_dict,
            )
        mock_find_install = mocker.patch(
            f'{_MODULE}.QtDependencyFinder.find_install_for_libs',
            return_value=install,
            )

        cli = ViPErLEEDGUICLI()
        # Qt works, so missing dependencies do not matter
        can_run = sanity is PyQtSanity.OK
        with nullcontext() if can_run else pytest.raises(SystemExit):
            cli.check_can_run_gui()
        missing = ('dep_one', 'dep_two', 'dep_three', 'dep_four')
        if sanity in self._no_deps_check:
            mock_find_deps.assert_not_called()
            return
        mock_find_deps.assert_called_once()
        mock_find_install.assert_called_once_with(missing)

        expect_err = self._err_msg + '''\
the following Qt dependencies of PyQt6 are missing on your system:
    dep_one
    dep_two
    dep_three
    dep_four
Try again after installing them'''
        expect_err += ('.' if not install else ''' via
    sudo apt install lib_one lib_two''')
        stderr = capsys.readouterr().err.strip()
        assert stderr.endswith(expect_err)

    def test_requested_binding(self, mocker, capsys):
        """Check that only a binding requested via QT_API is named."""
        mocker.patch(f'{_MODULE}.get_candidates', return_value=(PYQT5,))
        mocker.patch(f'{_MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.NOT_FOUND)
        with pytest.raises(SystemExit):
            ViPErLEEDGUICLI().check_can_run_gui()
        stderr = capsys.readouterr().err.strip()
        assert stderr.endswith('no Qt binding is installed (looked for '
                               'PyQt5).')

    @parametrize('sanity,err_msg', _can_run.values(), ids=_can_run)
    def test_result(self, sanity, err_msg, mock_sanity, mocker, capsys):
        """Check the expected outcome."""
        mock_sanity(sanity)
        find_deps = mocker.patch(f'{_MODULE}.find_missing_qt_dependencies')
        cli = ViPErLEEDGUICLI()
        context = pytest.raises(SystemExit) if err_msg else nullcontext()
        with context:
            cli.check_can_run_gui()
        captured = capsys.readouterr()
        # The line ends before any error is printed
        assert captured.out == self._checking + self._done
        stderr = captured.err.strip()
        if sanity is PyQtSanity.OK:
            assert not stderr
        else:
            expect = self._err_msg + err_msg
            assert stderr.endswith(expect)
        if sanity in self._no_deps_check:
            find_deps.assert_not_called()
        else:
            find_deps.assert_called_once()

    def test_sanity_invalid(self, mocker):
        """Check that an invalid sanity fails with a KeyError."""
        mocker.patch(f'{_MODULE}.check_pyqt_sanity', return_value='invalid')
        cli = ViPErLEEDGUICLI()
        with pytest.raises(KeyError):
            cli.check_can_run_gui()

    @parametrize(sanity=PyQtSanity)
    def test_sanity_msg(self, sanity):
        """Ensure there is a message for each PyQtSanity."""
        assert sanity in _SANITY_TO_ERR_MSG

    def test_too_old(self, mock_sanity, mocker, capsys):
        """Check complaints when the binding found is too old."""
        mock_sanity(PyQtSanity.TOO_OLD)
        find_deps = mocker.patch(f'{_MODULE}.find_missing_qt_dependencies')
        with pytest.raises(SystemExit):
            ViPErLEEDGUICLI().check_can_run_gui()
        stderr = capsys.readouterr().err.strip()
        expect = ('the installed PyQt6 is too old. ViPErLEED needs '
                  'PyQt5>=5.15 or PyQt6>=6.1. Upgrade it, e.g., via')
        assert self._err_msg + expect in stderr
        assert stderr.endswith('python -m pip install --upgrade PyQt6')
        find_deps.assert_not_called()


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
            'select_cls': mocker.MagicMock(name='ViPErLEEDSelectPlugin'),
            }
        mocker.patch(f'{_MODULE}.import_graphics_modules',
                     return_value=mock_imports.values())
        return {
            'catch_crash': mocker.patch(f'{_MODULE}.catch_gui_crash'),
            'font_path': mocker.patch(f'{_MODULE}.resources_path',
                                      return_value='fake/path/to/fonts'),
            **mock_imports,
            'app': mock_imports['qtw'].QApplication.return_value,
            'select': mock_imports['select_cls'].return_value,
            'suppress_warnings': mocker.patch(
                f'{_MODULE}.suppress_file_permission_warnings'
                ),
            }

    def test_execution(self, mocks, capsys):
        """Check the result of calling gui_main."""
        qapp = mocks['app']
        selector = mocks['select']
        result = gui_main()
        assert result is qapp.exec.return_value

        # Check important calls
        qapp.exec.assert_called_once_with()
        selector.show.assert_called_once_with()

        # Check also printing
        success_load_prints = 'Loading GUI...Done'
        stdout = capsys.readouterr().out
        assert stdout.rstrip() == success_load_prints

    def test_fonts_loaded(self, mocks, mocker):
        """Ensure the font database is updated."""
        gui_main()
        mocks['font_path'].assert_called_once_with('gui/fonts')
        base_path = Path(mocks['font_path'].return_value).resolve()
        font_db = mocks['qtg'].QFontDatabase.addApplicationFont
        assert font_db.mock_calls == [
            mocker.call(str(base_path/'DejaVuSans.ttf')),
            mocker.call(str(base_path/'cmunrm.otf')),
            ]

    @parametrize(binding=(PYQT5, PYQT6))
    def test_high_dpi_set(self, binding, mocks, mocker):
        """Ensure that the application supports high-DPI screens.

        Only with PyQt5. Qt 6 always supports them, and does not have
        the corresponding application attributes any longer.
        """
        mocker.patch(f'{_MODULE}.get_binding', return_value=binding)
        gui_main()
        attribute = mocks['qtc'].Qt.ApplicationAttribute
        qapp = mocks['qtg'].QGuiApplication
        set_high_dpi = qapp.setAttribute
        expect_calls = [
            mocker.call(attribute.AA_EnableHighDpiScaling),
            mocker.call(attribute.AA_UseHighDpiPixmaps),
            ] if binding == PYQT5 else []
        assert set_high_dpi.mock_calls == expect_calls

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
        mocker.patch(f'{_MODULE}.has_pyqt', return_value=True)
        mock_import = mocker.patch(f'{_MODULE}.import_module')
        result = import_graphics_modules()
        imported_modules = (
            'viperleed.gui.selectplugin',
            f'{_QT_COMPAT}.QtCore',
            f'{_QT_COMPAT}.QtGui',
            f'{_QT_COMPAT}.QtWidgets',
            )
        assert mock_import.mock_calls == [mocker.call(m)
                                          for m in imported_modules]
        n_return_values = 4
        assert len(result) == n_return_values


class TestViPErLEEDGUICLI:
    """Tests for the ViPErLEEDGUICLI class."""

    _valid_args = (
        ([], {'nogui': False, 'skip_graphics_check': False}),
        (['--nogui'], {'nogui': True, 'skip_graphics_check': False}),
        (['--skip-graphics-check'],
         {'nogui': False, 'skip_graphics_check': True}),
        (['--nogui', '--skip-graphics-check'],
         {'nogui': True, 'skip_graphics_check': True}),
        )

    @parametrize('args,attrs', _valid_args)
    def test_parse_args_valid(self, args, attrs):
        """Check result of successfully parsing command-line arguments."""
        cli = ViPErLEEDGUICLI()
        parsed = cli.parse_cli_args(args)
        for attr_name, attr_value in attrs.items():
            assert getattr(parsed, attr_name) == attr_value

    @fixture(name='mocks')
    def fixture_mock_calls(self, mocker):
        """Replace the functions called by ViPErLEEDGUICLI with mocks."""
        cli = ViPErLEEDGUICLI()
        mock_check_sanity = mocker.patch(f'{_MODULE}.check_pyqt_sanity')
        mocks = mocker.MagicMock()  # Also records the order of calls
        mocks.attach_mock(mocker.patch(f'{_MODULE}.commandline_main'),
                          'commandline_main')
        mocks.attach_mock(mocker.patch(f'{_MODULE}.gui_main'), 'gui_main')
        mocks.attach_mock(mocker.patch.object(cli, 'check_can_run_gui'),
                          'check_can_run_gui')
        mocks.attach_mock(mock_check_sanity.assume_ok, 'assume_ok')
        return cli, mocks

    @parametrize(skip_check=(False, True))
    def test_call_commandline_mode(self, skip_check, mocks, mocker):
        """Check inner calls when running in command-line mode."""
        cli, mocks = mocks
        args = ['--nogui', *(['--skip-graphics-check'] if skip_check else [])]
        result = cli(args)
        assert result is mocks.commandline_main.return_value
        assert mocks.mock_calls == [mocker.call.commandline_main()]

    def test_call_gui_mode(self, mocks, mocker):
        """Check inner calls when running in GUI mode."""
        cli, mocks = mocks
        result = cli([])
        assert result is mocks.gui_main.return_value
        assert mocks.mock_calls == [
            mocker.call.check_can_run_gui(),
            mocker.call.gui_main(),
            ]

    def test_call_skip_graphics_check(self, mocks, mocker):
        """Check inner calls when skipping the check of graphics."""
        cli, mocks = mocks
        mock_process = mocker.patch('viperleed.gui.detect_graphics.Process')
        result = cli(['--skip-graphics-check'])
        assert result is mocks.gui_main.return_value
        # Sanity must be assumed before starting the GUI, so
        # that the GUI does not check it again via has_pyqt
        assert mocks.mock_calls == [
            mocker.call.assume_ok(),
            mocker.call.gui_main(),
            ]
        mock_process.assert_not_called()

    def test_call_skip_graphics_check_qt_api(self, mocks, mocker, capsys):
        """Check complaints about QT_API also when skipping checks."""
        cli, mocks = mocks
        mocker.patch(f'{_MODULE}.get_candidates',
                     side_effect=QtApiError('Some reason.'))
        with pytest.raises(SystemExit):
            cli(['--skip-graphics-check'])
        assert not mocks.mock_calls
        stderr = capsys.readouterr().err.strip()
        assert stderr.endswith('no Qt binding can be selected. Some reason.')

    def test_skip_check_prevents_later_checks(self, mocker):
        """Ensure that has_pyqt does not check after skipping the check."""
        detect_graphics = 'viperleed.gui.detect_graphics'
        mocker.patch(f'{detect_graphics}.PyQtSanityChecker._summary', None)
        mock_process = mocker.patch(f'{detect_graphics}.Process')
        mock_gui_main = mocker.patch(f'{_MODULE}.gui_main',
                                     side_effect=has_pyqt)
        cli = ViPErLEEDGUICLI()
        assert cli(['--skip-graphics-check'])
        mock_gui_main.assert_called_once_with()
        mock_process.assert_not_called()
