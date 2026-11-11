"""Tests for module detect_graphics of viperleed.gui."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-05-22'
__license__ = 'GPLv3+'

import os
from pathlib import Path
from subprocess import SubprocessError
from types import SimpleNamespace

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.detect_graphics import _DISPLAY_TIMEOUT
from viperleed.gui.detect_graphics import _EXIT_CODE_OFFSET
from viperleed.gui.detect_graphics import _QtUnixPlatformsGetter
from viperleed.gui.detect_graphics import _STARTUP_TIMEOUT
from viperleed.gui.detect_graphics import PyQtSanity
from viperleed.gui.detect_graphics import PyQtSanityChecker
from viperleed.gui.detect_graphics import QtDependencyFinder
from viperleed.gui.detect_graphics import _import_dummy_pyqt_and_run_app
from viperleed.gui.detect_graphics import _init_dummy_qapp
from viperleed.gui.detect_graphics import find_missing_qt_dependencies
from viperleed.gui.detect_graphics import has_graphics
from viperleed.gui.detect_graphics import has_pyqt
from viperleed.gui.detect_graphics import suppress_file_permission_warnings
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6

_MODULE = 'viperleed.gui.detect_graphics'
BOTH = (PYQT6, PYQT5)  # In the default order of preference


# TODO: would be better to find a way to mock at least
# _import_dummy_pyqt_and_run_app instead of mocking
# multiprocessing.Process. My tests seem to reveal
# that, while _check_pyqt_sanity gets the correct
# target, that target is never actually called.


@fixture(name='mock_multiprocessing')
def fixture_mock_multiprocessing(mocker):
    """Replace the multiprocessing tools of the checker with mocks.

    Notice that it is necessary to do so in order to correctly
    test the behavior of has_display. In fact, we can't mock
    detect_graphics.QApplication as the module is re-imported
    by multiprocessing. We could, in principle, mock the target
    of multiprocessing.Process, i.e., _init_dummy_qapp. However,
    I tried without success to do so.

    The returned factory takes the following arguments:
    ready : {'imports', 'exited', 'none'}
        Which of the waited-for objects is ready within the start-up
        time: the pipe (a binding was imported), the process sentinel
        (the subprocess exited before), or none of them.
    alive : bool
        Whether the subprocess is still alive after the QApplication
        time limit, i.e., whether creating the QApplication hangs.
    exitcode : object
        The exit code of the subprocess.
    """
    def _mock(ready='imports', alive=False, exitcode=0):
        imports_done = mocker.MagicMock(name='imports_done')
        imports_done.recv.return_value = PYQT5
        child_end = mocker.MagicMock(name='child_end')
        process = mocker.MagicMock(name='process', exitcode=exitcode)
        process.is_alive.return_value = alive
        ready_objects = {
            'imports': [imports_done],
            'exited': [process.sentinel],
            'none': [],
            }
        order = []  # Names of selected calls, in the order they're made
        process.start.side_effect = lambda: order.append('start')
        child_end.close.side_effect = lambda: order.append('close_child')

        def _wait(*_):
            order.append('wait')
            return ready_objects[ready]
        return SimpleNamespace(
            order=order,
            imports_done=imports_done,
            child_end=child_end,
            process=process,
            pipe=mocker.patch(f'{_MODULE}.Pipe',
                              return_value=(imports_done, child_end)),
            process_cls=mocker.patch(f'{_MODULE}.Process',
                                     return_value=process),
            wait=mocker.patch(f'{_MODULE}.wait', side_effect=_wait),
            candidates=mocker.patch(f'{_MODULE}.get_candidates',
                                    return_value=BOTH),
            )
    return _mock


class TestInitDummyQapp:
    """Tests for the _init_dummy_qapp subprocess worker."""

    @fixture(name='binding')
    @parametrize(binding=BOTH)
    def fixture_binding(self, binding):
        """Return the name of a Qt binding."""
        return binding

    @fixture(name='mock_qapp')
    def fixture_mock_qapp(self, binding, mocker, mock_import_module):
        """Replace QApplication of `binding` with a mock."""
        mocker.patch(f'{_MODULE}.suppress_file_permission_warnings')
        # Fake only QtWidgets: newer Pythons (e.g., 3.13) find the
        # targets of mock.patch via importlib.import_module, and later
        # patches (e.g., in test_calls) must reach the real modules.
        mock_qtwidgets = mocker.MagicMock()
        mock_import_module(**{f'{binding}.QtWidgets': mock_qtwidgets})
        return mock_qtwidgets.QApplication

    def test_calls(self, binding, mock_qapp, mocker):
        """Check calls in _init_dummy_qapp."""
        mock_suppress = mocker.patch(
            f'{_MODULE}.suppress_file_permission_warnings'
            )
        with pytest.raises(SystemExit):
            _init_dummy_qapp(binding)
        mock_suppress.assert_called_once_with(binding)
        # QApplication is created only if the QtWidgets of binding
        # is imported, as the QtWidgets of others are not mocked
        mock_qapp.assert_called_once_with([])

    _screen = {
        'with screen': (object(), PyQtSanity.OK),
        'no screen': (None, PyQtSanity.NO_DISPLAY),
        }

    @parametrize('screen,expect', _screen.values(), ids=_screen)
    def test_exit_code(self, binding, screen, expect, mock_qapp):
        """Check the exit code with and without a primary screen."""
        mock_qapp.return_value.primaryScreen.return_value = screen
        with pytest.raises(SystemExit) as excinfo:
            _init_dummy_qapp(binding)
        assert excinfo.value.code == _EXIT_CODE_OFFSET + expect.value


@parametrize(pyqt_present=(True, False))
@parametrize(has_display=(True, False))
def test_has_graphics(pyqt_present, has_display, mocker):
    """Ensure correct identification of graphics capabilities."""
    mocker.patch(f'{_MODULE}.has_pyqt', return_value=pyqt_present)
    check_sanity = mocker.patch(
        f'{_MODULE}.check_pyqt_sanity',
        return_value=PyQtSanity.OK if has_display else PyQtSanity.NO_DISPLAY
        )
    result = has_graphics()
    assert result == (pyqt_present and has_display)

    # Check also internal calls
    if pyqt_present:
        check_sanity.assert_called_once_with()
    else:
        check_sanity.assert_not_called()


@parametrize(pyqt_sanity=PyQtSanity)
def test_has_pyqt(pyqt_sanity, mocker):
    """Ensure correct detection of the presence of a Qt binding."""
    mocker.patch(f'{_MODULE}.check_pyqt_sanity', return_value=pyqt_sanity)
    pyqt_present = pyqt_sanity in {PyQtSanity.OK, PyQtSanity.NO_DISPLAY}
    assert has_pyqt() == pyqt_present


class TestImportDummyPyqtAndRunApp:
    """Tests for the _import_dummy_pyqt_and_run_app subprocess worker."""

    # The modules that we attempt to import, besides the package
    _modules = (
        'QtCore',
        'QtWidgets',
        )

    @fixture(name='mock_dummy_qapp', autouse=True)
    def fixture_mock_dummy_qapp(self, mocker):
        """Replace _init_dummy_qapp with a mock."""
        # This is to ensure tests don't hang forever on a system
        # without graphics.
        return mocker.patch(f'{_MODULE}._init_dummy_qapp')

    @fixture(name='imports_done')
    def fixture_imports_done(self, mocker):
        """Return a fake sending end of a pipe."""
        return mocker.MagicMock(name='imports_done')

    @staticmethod
    def check_exit_code(excinfo, expect):
        """Check that the exit code of `excinfo` is for `expect`."""
        assert excinfo.value.code == _EXIT_CODE_OFFSET + expect.value

    def test_cant_import_pyqt(self, imports_done, mocker):
        """Check the expected exit code when importing any binding fails."""
        mock_import = mocker.patch('importlib.import_module',
                                   side_effect=ImportError)
        with pytest.raises(SystemExit) as excinfo:
            _import_dummy_pyqt_and_run_app(imports_done, BOTH)
        self.check_exit_code(excinfo, PyQtSanity.IMPORT_ERROR)
        imports_done.send.assert_not_called()
        assert mock_import.mock_calls == [mocker.call(b) for b in BOTH]

    # Exceptions raised when importing a submodule of a binding, as a
    # function of the submodule's name, and the expected exit status.
    _submodule_errors = {
        'import error': (lambda _: ImportError(),
                         PyQtSanity.IMPORT_ERROR),
        # E.g., only PyQt5-sip and PyQt5-Qt5 are left after
        # uninstalling PyQt5: 'PyQt5' is a namespace package
        'submodule missing': (lambda mod: ModuleNotFoundError(name=mod),
                              PyQtSanity.NOT_FOUND),
        # A dependency of the submodule is missing: broken PyQt5
        'dependency missing': (
            lambda mod: ModuleNotFoundError(name=f'{mod.split(".")[0]}.sip'),
            PyQtSanity.IMPORT_ERROR,
            ),
        }

    @parametrize(binding=BOTH)
    @parametrize(module=_modules)
    @parametrize('make_exc,expect', _submodule_errors.values(),
                 ids=_submodule_errors)
    # pylint: disable-next=too-many-arguments  # 2/7 fixtures
    def test_cant_import_pyqt_submodule(self, binding, module, make_exc,
                                        expect, imports_done, mocker):
        """Check the expected exit code when importing a submodule fails."""
        failing = f'{binding}.{module}'
        def _fail_module_import(mod_name):
            if mod_name == failing:
                raise make_exc(mod_name)
        mocker.patch('importlib.import_module',
                     side_effect=_fail_module_import)
        with pytest.raises(SystemExit) as excinfo:
            _import_dummy_pyqt_and_run_app(imports_done, (binding,))
        self.check_exit_code(excinfo, expect)
        imports_done.send.assert_not_called()

    def test_pyqt_not_found(self, imports_done, mocker):
        """Check the expected exit code when no binding is found."""
        mocker.patch('importlib.import_module',
                     side_effect=ModuleNotFoundError)
        with pytest.raises(SystemExit) as excinfo:
            _import_dummy_pyqt_and_run_app(imports_done, BOTH)
        self.check_exit_code(excinfo, PyQtSanity.NOT_FOUND)
        imports_done.send.assert_not_called()

    def test_import_error_wins(self, imports_done, mocker):
        """Check that one broken binding gives IMPORT_ERROR."""
        def _import(name):
            if name.startswith(PYQT6):
                raise ImportError
            raise ModuleNotFoundError
        mocker.patch('importlib.import_module', side_effect=_import)
        with pytest.raises(SystemExit) as excinfo:
            _import_dummy_pyqt_and_run_app(imports_done, BOTH)
        self.check_exit_code(excinfo, PyQtSanity.IMPORT_ERROR)

    @parametrize(binding=BOTH)
    def test_success(self, binding, imports_done, mocker, mock_dummy_qapp):
        """Test an execution that does not cause the interpreter to exit."""
        # NB: normally, the interpreter always exits when
        # _import_dummy_pyqt_and_run_app is called. The only
        # reason why we don't exit here is the patching of
        # _init_dummy_qapp via mock_dummy_qapp.
        mock_import = mocker.patch('importlib.import_module')
        calls = mocker.MagicMock()
        calls.attach_mock(imports_done, 'imports_done')
        calls.attach_mock(mock_dummy_qapp, 'init_dummy_qapp')
        _import_dummy_pyqt_and_run_app(imports_done, (binding,))

        expect_calls = [
            mocker.call(binding),
            *(mocker.call(f'{binding}.{m}') for m in self._modules)
            ]
        assert mock_import.mock_calls == expect_calls
        # The parent must be told before creating the QApplication.
        # Build the expected calls outside the assert: code-quality
        # checkers take mocker.call.<...>.send() for a real call.
        expect_order = [
            mocker.call.imports_done.send(binding),
            mocker.call.imports_done.close(),
            mocker.call.init_dummy_qapp(binding),
            ]
        assert calls.mock_calls == expect_order

    def test_falls_back(self, imports_done, mocker, mock_dummy_qapp):
        """Check that the next candidate is used if the first fails."""
        def _import(name):
            if name.startswith(PYQT6):
                raise ModuleNotFoundError
        mocker.patch('importlib.import_module', side_effect=_import)
        _import_dummy_pyqt_and_run_app(imports_done, BOTH)
        imports_done.send.assert_called_once_with(PYQT5)
        mock_dummy_qapp.assert_called_once_with(PYQT5)


class TestFindMissingDependencies:
    """Tests for the find_missing_qt_dependencies function."""

    @fixture(autouse=True)
    def mock_sane_pyqt(self, mocker):
        """Do not start a real QApplication to check Qt."""
        # Unmocked, the first call would run a subprocess that starts
        # a QApplication, which aborts (and makes pytest's faulthandler
        # print "Fatal Python error") on systems without a display.
        mocker.patch(f'{_MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.OK,
                     binding=PYQT5)

    def test_raises(self, mocker):
        """Check complaints when calling the function without Qt."""
        mocker.patch(f'{_MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.NOT_FOUND)
        with pytest.raises(NotImplementedError):
            find_missing_qt_dependencies()

    def test_windows(self, mocker):
        """Check the result of finding dependencies on Windows."""
        mocker.patch('sys.platform', 'win32')
        missing = find_missing_qt_dependencies()
        assert not missing

    def test_not_implemented(self, mocker):
        """Check the result when dependency finding is not implemented."""
        mocker.patch('sys.platform', 'some-unsupported-platform')
        mock_finder = mocker.MagicMock()
        mocker.patch(f'{_MODULE}.QtDependencyFinder',
                     return_value=mock_finder)
        mock_finder.find_missing_dependencies.side_effect = NotImplementedError
        missing = find_missing_qt_dependencies()
        assert not missing

    def test_missing_deps(self, mocker):
        """Check the expected result when missing dependencies are found."""
        mocker.patch('sys.platform', 'some-supported-platform')
        mock_finder = mocker.MagicMock()
        mock_finder_cls = mocker.patch(f'{_MODULE}.QtDependencyFinder',
                                       return_value=mock_finder)
        missing = find_missing_qt_dependencies()
        assert missing is mock_finder.find_missing_dependencies.return_value
        mock_finder_cls.assert_called_once_with(PYQT5)

    @parametrize(installed=BOTH)
    def test_binding_not_imported(self, installed, mocker):
        """Check that the first installed binding is used if none worked."""
        mocker.patch('sys.platform', 'some-supported-platform')
        mocker.patch(f'{_MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.IMPORT_ERROR,
                     binding=None)
        mocker.patch(f'{_MODULE}.get_candidates', return_value=BOTH)
        mocker.patch('importlib.util.find_spec',
                     side_effect=lambda name: name == installed or None)
        mock_finder_cls = mocker.patch(f'{_MODULE}.QtDependencyFinder')
        find_missing_qt_dependencies()
        mock_finder_cls.assert_called_once_with(installed)

    def test_no_binding_installed(self, mocker):
        """Check the result when no binding is installed at all."""
        mocker.patch('sys.platform', 'some-supported-platform')
        mocker.patch(f'{_MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.IMPORT_ERROR,
                     binding=None)
        mocker.patch('importlib.util.find_spec', return_value=None)
        mock_finder_cls = mocker.patch(f'{_MODULE}.QtDependencyFinder')
        assert not find_missing_qt_dependencies()
        mock_finder_cls.assert_not_called()


class TestSuppressFilePermissionWarnings:
    """Tests for the suppress_file_permission_warnings function."""

    def test_no_xdg_runtime_dir(self, mocker):
        """Check no changes if XDG_RUNTIME_DIR is not set (e.g., not UNIX)."""
        mocker.patch(f'{_MODULE}._UNIX_RUNTIME_DIR', '')
        mock_mkdir = mocker.patch('pathlib.Path.mkdir')
        mocker.patch.dict('os.environ', clear=True)

        suppress_file_permission_warnings()
        mock_mkdir.assert_not_called()
        # pylint: disable-next=magic-value-comparison
        assert 'XDG_RUNTIME_DIR' not in os.environ

    def test_xdg_set(self, tmp_path, mocker):
        """Check implementation when XDG_RUNTIME_DIR is set."""
        expect_kwargs = {'mode': 0o0700, 'parents': True, 'exist_ok': True}
        mkdir = Path.mkdir
        def mock_mkdir(*args, **kwargs):
            assert kwargs == expect_kwargs
            mkdir(*args, **kwargs)

        mocker.patch(f'{_MODULE}._UNIX_RUNTIME_DIR', str(tmp_path))
        mocker.patch('pathlib.Path.mkdir', mock_mkdir)
        expect_path = (tmp_path / 'viperleed-gui').resolve()
        assert not expect_path.exists()

        suppress_file_permission_warnings()
        assert expect_path.is_dir()
        assert os.environ['XDG_RUNTIME_DIR'] == str(expect_path)

    def test_mkdir_fails(self, tmp_path, mocker):
        """Check setting of fallback directory if the default fails."""
        fails = tmp_path/'fails'
        fallback = tmp_path/'fallback'
        fallback.mkdir()  # Otherwise kwargs['parents'] is mismatched
        expect_kwargs = {'mode': 0o0700, 'parents': True, 'exist_ok': True}
        not_created = fails/'viperleed-gui'
        expect_path = fallback/'viperleed-gui'

        mocker.patch(f'{_MODULE}._UNIX_RUNTIME_DIR', str(fails))

        mkdir = Path.mkdir
        def mkdir_raises(path, *args, **kwargs):
            check_kwargs = kwargs.copy()
            if path in (not_created, expect_path):
                if args:
                    check_kwargs['mode'] = args[0]
                assert check_kwargs == expect_kwargs
            if path == not_created:
                raise OSError
            mkdir(path, *args, **kwargs)

        mocker.patch('pathlib.Path.mkdir', mkdir_raises)
        mock_import = mocker.patch('importlib.import_module')
        mock_qtc = mock_import.return_value
        mock_qtc.QStandardPaths.writableLocation.return_value = str(fallback)

        suppress_file_permission_warnings(PYQT6)
        mock_import.assert_called_once_with('PyQt6.QtCore')
        assert not not_created.exists()
        assert expect_path.is_dir()
        assert os.environ['XDG_RUNTIME_DIR'] == str(expect_path)


class TestPyQtSanityChecker:
    """Tests for the PyQtSanityChecker class."""

    @fixture(autouse=True)
    def fresh_checker(self, mocker):
        """Forget previous checks, and do not select bindings globally."""
        mocker.patch.object(PyQtSanityChecker, '_status', None)
        mocker.patch.object(PyQtSanityChecker, '_binding', None)
        return mocker.patch(f'{_MODULE}.set_binding')

    def test_checks_only_once(self, mocker):
        """Ensure that checking is skipped when already done."""
        mock_sanity = mocker.MagicMock()
        mock_do_check = mocker.patch.object(PyQtSanityChecker,
                                            '_check_pyqt_sanity',
                                            return_value=(mock_sanity, None))
        checker = PyQtSanityChecker()

        # Call it once: performs the check.
        assert checker() is mock_sanity
        # pylint: disable-next=protected-access
        assert PyQtSanityChecker._status is mock_sanity

        # Call it again: should skip the check
        assert checker() is mock_sanity
        mock_do_check.assert_called_once_with()

    @parametrize(binding=BOTH)
    def test_selects_binding(self, binding, fresh_checker, mocker):
        """Check that the binding found is selected for qt_compat."""
        mocker.patch.object(PyQtSanityChecker, '_check_pyqt_sanity',
                            return_value=(PyQtSanity.OK, binding))
        checker = PyQtSanityChecker()
        assert checker() is PyQtSanity.OK
        assert checker.binding == binding
        fresh_checker.assert_called_once_with(binding)

    def test_no_binding(self, fresh_checker, mocker):
        """Check that no binding is selected if none was found."""
        mocker.patch.object(PyQtSanityChecker, '_check_pyqt_sanity',
                            return_value=(PyQtSanity.NOT_FOUND, None))
        checker = PyQtSanityChecker()
        assert checker() is PyQtSanity.NOT_FOUND
        assert checker.binding is None
        fresh_checker.assert_not_called()

    def test_forget_status(self, mocker):
        """Ensure resetting re-checks."""
        mock_do_check = mocker.patch.object(PyQtSanityChecker,
                                            '_check_pyqt_sanity',
                                            return_value=(None, PYQT5))
        checker = PyQtSanityChecker()
        # Check once to ensure forget_status does the right thing
        checker()

        checker.forget_status()
        assert checker.binding is None
        checker()

        # pylint: disable-next=magic-value-comparison
        assert mock_do_check.call_count == 2

    def test_assume_ok(self, mocker):
        """Ensure that assuming a sane Qt binding skips checking."""
        mock_sanity = mocker.MagicMock()
        mock_do_check = mocker.patch.object(PyQtSanityChecker,
                                            '_check_pyqt_sanity',
                                            return_value=(mock_sanity, None))
        checker = PyQtSanityChecker()
        checker.assume_ok()
        assert checker() is PyQtSanity.OK
        assert checker.binding is None  # qt_compat selects one
        mock_do_check.assert_not_called()

        # forget_status undoes it
        checker.forget_status()
        assert checker() is mock_sanity
        mock_do_check.assert_called_once_with()

    _exit_codes = {
        # Offset exit codes are mapped to the PyQtSanity members
        **{_EXIT_CODE_OFFSET + s.value: s for s in PyQtSanity},
        # Any other exit code is a crash. In particular, the ones that
        # coincide with the values of PyQtSanity members, like 1 for
        # an unhandled exception, must not be taken as such members.
        **{s.value: PyQtSanity.RUNTIME_CRASH for s in PyQtSanity},
        -6: PyQtSanity.RUNTIME_CRASH,   # SIGABRT
        134: PyQtSanity.RUNTIME_CRASH,  # SIGABRT on Unix
        _EXIT_CODE_OFFSET - 1: PyQtSanity.RUNTIME_CRASH,
        _EXIT_CODE_OFFSET + len(PyQtSanity): PyQtSanity.RUNTIME_CRASH,
        }

    @staticmethod
    def check_common_calls(mocks):
        """Check the calls made whatever the outcome of the check."""
        mocks.pipe.assert_called_once_with(duplex=False)
        mocks.process_cls.assert_called_once_with(
            target=_import_dummy_pyqt_and_run_app,
            args=(mocks.child_end, BOTH),
            )
        mocks.process.start.assert_called_once_with()
        mocks.child_end.close.assert_called_once_with()
        mocks.wait.assert_called_once_with(
            (mocks.imports_done, mocks.process.sentinel),
            _STARTUP_TIMEOUT,
            )
        mocks.imports_done.close.assert_called_once_with()
        # The parent closes its copy of the child's end of the
        # pipe right after starting, i.e., before waiting.
        assert mocks.order == ['start', 'close_child', 'wait']

    @parametrize(ready=('imports', 'exited'))
    @parametrize('exit_code,expect', _exit_codes.items())
    def test_check_sanity_exits(self, ready, exit_code, expect,
                                mock_multiprocessing):
        """Check the result when the subprocess exits by itself."""
        mocks = mock_multiprocessing(ready=ready, exitcode=exit_code)
        # pylint: disable-next=protected-access
        sanity, binding = PyQtSanityChecker._check_pyqt_sanity()
        assert sanity is expect
        # The binding is known only if it was sent before exiting,
        # and it is reported only if it can be used
        usable = expect in {PyQtSanity.OK, PyQtSanity.NO_DISPLAY}
        expect_binding = PYQT5 if ready == 'imports' and usable else None
        assert binding == expect_binding
        self.check_common_calls(mocks)
        mocks.process.join.assert_called_once_with(_DISPLAY_TIMEOUT)
        mocks.process.terminate.assert_not_called()

    def test_check_sanity_no_display(self, mock_multiprocessing):
        """Check the result when creating a QApplication hangs."""
        mocks = mock_multiprocessing(ready='imports', alive=True)
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_pyqt_sanity()
        assert result == (PyQtSanity.NO_DISPLAY, PYQT5)
        self.check_common_calls(mocks)
        mocks.process.join.assert_called_once_with(_DISPLAY_TIMEOUT)
        mocks.process.terminate.assert_called_once_with()

    def test_check_sanity_no_message(self, mock_multiprocessing):
        """Check the result when the subprocess dies without sending."""
        mocks = mock_multiprocessing(ready='imports',
                                     exitcode=-6)  # SIGABRT
        mocks.imports_done.recv.side_effect = EOFError
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_pyqt_sanity()
        assert result == (PyQtSanity.RUNTIME_CRASH, None)

    def test_check_sanity_unresponsive(self, mock_multiprocessing):
        """Check the result when importing does not finish in time."""
        mocks = mock_multiprocessing(ready='none', alive=True)
        # pylint: disable-next=protected-access
        result = PyQtSanityChecker._check_pyqt_sanity()
        assert result == (PyQtSanity.UNRESPONSIVE, None)
        self.check_common_calls(mocks)
        mocks.process.join.assert_not_called()
        mocks.process.terminate.assert_called_once_with()


@fixture(name='mock_pyqt_root')
def fixture_mock_pyqt_root(mocker, mock_import_module):
    """Fake the existence of the root path of PyQt5 and PyQt6."""
    mock_root = '/path/to/pyqt'
    mock_pyqt = mocker.MagicMock()
    mock_pyqt.__file__ = f'{mock_root}/__init__.py'
    mock_import_module(PyQt5=mock_pyqt, PyQt6=mock_pyqt)
    mocker.patch('pathlib.Path.exists', return_value=True)
    return Path(mock_root).resolve()


class TestQtDependencyFinder:
    """Tests for the QtDependencyFinder helper class."""

    @parametrize('binding,qt_dir', ((PYQT5, 'Qt5'), (PYQT6, 'Qt6')))
    def test_init(self, binding, qt_dir, mock_pyqt_root):
        """Check initialization of a finder."""
        expect_root = mock_pyqt_root / f'{qt_dir}/plugins/platforms'
        finder = QtDependencyFinder(binding)
        # pylint: disable-next=protected-access
        assert finder._root == expect_root

    _install_linux = {
        'found': (('lib_one', 'lib_two'), 'sudo apt install lib_one lib_two'),
        'not found': ((), ''),
        }

    @parametrize('libs,expect', _install_linux.values(), ids=_install_linux)
    def test_install_libs_linux(self, libs, expect, mocker):
        """Check result of _find_install_for_libs_linux."""
        mock_list = mocker.patch(
            f'{_MODULE}.QtDependencyFinder._list_install_for_libs_linux',
            return_value=libs,
            )
        missing = mocker.MagicMock()
        # pylint: disable-next=protected-access
        result = QtDependencyFinder._find_install_for_libs_linux(missing)
        mock_list.assert_called_once_with(missing)
        assert result == expect
        mocker.patch('sys.platform', 'linux')
        assert QtDependencyFinder.find_install_for_libs(missing) == result

    _unsupported = (
        'aix',     # Python >= 3.8
        'aix5',    # Python < 3.8
        'aix7',    # Python < 3.8
        'atheos',
        'cygwin',
        'darwin',  # macOS
        'win32',
        'msys',
        'os2',
        'os2emx',
        'riscos',
        'freebsd7',
        'freebsd8',
        'freebsdN',
        'openbsd6',
        'invalidplatform',
        )
    _supported = (
        'linux',
        'validplatform',
        )

    @pytest.mark.usefixtures('mock_pyqt_root')
    @parametrize(platform=_supported)
    def test_find_install_for_libs_supported(self, platform, mocker):
        """Check result of find_install_for_libs on supported platforms."""
        mocker.patch('sys.platform', platform)
        mock_impl = mocker.patch(
            # Patch the class, as find_install_for_libs is @classmethod
            f'{_MODULE}.QtDependencyFinder._find_install_for_libs_{platform}',
            create=True,
            )
        missing = mocker.MagicMock()
        finder = QtDependencyFinder(PYQT5)
        result = finder.find_install_for_libs(missing)
        assert result is mock_impl.return_value
        mock_impl.assert_called_once_with(missing)

    @parametrize(platform=_unsupported)
    def test_find_install_for_libs_unsupported(self, platform, mocker):
        """Check result of find_install_for_libs on unsupported platforms."""
        mocker.patch('sys.platform', platform)
        result = QtDependencyFinder.find_install_for_libs(mocker.MagicMock())
        no_install_suggestion = ''
        assert result == no_install_suggestion

    @parametrize(platform=(*_supported, *_unsupported))
    def test_find_install_for_libs_no_missing(self, platform, mocker):
        """Check shortcut if nothing is missing."""
        mocker.patch('sys.platform', platform)
        mock_delegate = mocker.patch.object(
            QtDependencyFinder,
            f'_find_install_for_libs_{platform}',
            create=True,
            )
        none_missing = {}
        assert not QtDependencyFinder.find_install_for_libs(none_missing)
        mock_delegate.assert_not_called()

    @pytest.mark.usefixtures('mock_pyqt_root')
    @parametrize(platform=_supported)
    def test_find_missing_dependencies_supported(self, platform, mocker):
        """Check result of find_missing_dependencies on supported platforms."""
        mocker.patch('sys.platform', platform)
        finder = QtDependencyFinder(PYQT5)
        mock_impl = mocker.patch.object(
            finder,
            f'_find_missing_deps_{platform}',
            create=True,
            )
        result = finder.find_missing_dependencies()
        assert result is mock_impl.return_value
        mock_impl.assert_called_once_with()

    @pytest.mark.usefixtures('mock_pyqt_root')
    @parametrize(platform=_unsupported)
    def test_find_missing_dependencies_unsupported(self, platform, mocker):
        """Check find_missing_dependencies raises on unsupported platforms."""
        mocker.patch('sys.platform', platform)
        finder = QtDependencyFinder(PYQT5)
        with pytest.raises(NotImplementedError):
            finder.find_missing_dependencies()

    @pytest.mark.usefixtures('mock_pyqt_root')
    @parametrize(exc=(FileNotFoundError, SubprocessError))
    def test_find_missing_dependencies_implementation_fails(self, exc, mocker):
        """Check find_missing_dependencies when implementation raises."""
        mocker.patch('sys.platform', 'linux')
        finder = QtDependencyFinder(PYQT5)
        mocker.patch.object(finder,
                            '_find_missing_deps_linux',
                            side_effect=exc)
        assert not finder.find_missing_dependencies()


class TestQtDependencyFinderLinux:
    """Tests for finding missing dependencies on Linux."""

    @pytest.mark.usefixtures('mock_pyqt_root')
    def test_find_missing_deps_linux(self, mocker):
        """Check the result of _find_missing_deps_linux."""
        finder = QtDependencyFinder(PYQT5)
        mock_get = mocker.patch(f'{_MODULE}._QtUnixPlatformsGetter.get',
                                return_value=('missing_deps',
                                              'no_missing_deps'))
        root = finder._root  # pylint: disable=protected-access
        so_files = (
            root/'libmissing_deps.so',
            root/'anothermissing_deps.so.999',
            root/'nested/missing_deps.so.12',
            root/'no_missing_deps.so.9',
            )
        dependencies = {
            str(root/'libmissing_deps.so'): '''\
dep_one.so => not found
dep_one.so => not found, it's a duplicate
dep_two.so.3 => dep was found, here the info\n''',
            str(root/'anothermissing_deps.so.999'): '''\
dep_three.so.9 (0x00007ffd6f7d7000) => not found, hex code at left\n''',
            str(root/'nested/missing_deps.so.12'): '''\
nested/deeper/dep_four => not found\n''',
            str(root/'no_missing_deps.so.9'): '''\
dep_five => was found
dep_six.so => this one is also already installed''',
            }
        expect = {
            'libmissing_deps.so': ['dep_one'],
            'anothermissing_deps.so.999': ['dep_three'],
            'nested/missing_deps.so.12': ['dep_four'],
            }

        def mock_glob(_, pattern):
            pattern = pattern.replace('*', '')
            for file in so_files:
                if pattern in file.name:
                    yield file

        def mock_ldd(cmd, **_):
            """Fake the result of subprocess.run(['ldd' ...])."""
            _, so_name = cmd
            mock_process = mocker.MagicMock()
            mock_process.stdout = dependencies[so_name].encode()
            return mock_process

        mocker.patch('pathlib.Path.rglob', mock_glob)
        mocker.patch('subprocess.run', mock_ldd)

        # pylint: disable-next=protected-access
        missing = finder._find_missing_deps_linux()
        assert missing == expect
        mock_get.assert_called_once_with(PYQT5)

    def test_list_install_libs_success(self, mocker):
        """Check the successful listing of apt modules."""
        missing = ['multiple', 'devonly', 'another']
        expect_rgx = 'multiple|devonly|another'
        def mock_apt(cmd, **__):
            """Fake the result of calling 'apt-cache search'."""
            *_, rgx = cmd
            assert rgx == expect_rgx
            mock_proc = mocker.MagicMock()
            mock_proc.stdout = '''\
multiplelib-one - description -- some more text
multiplelib-one4 - more recent than the previous - should be selected
multiplelib-one5-dev - this one looks more recent but it is a development one
devonly-three-dev - this comes only with developer options
another-lib - this has no number'''.encode()
            return mock_proc

        expect = ['multiplelib-one4',
                  'devonly-three-dev',
                  'another-lib']
        mocker.patch('subprocess.run', mock_apt)
        # pylint: disable-next=protected-access
        libs = QtDependencyFinder._list_install_for_libs_linux(missing)
        assert libs == expect

    def test_list_install_libs_one_missing(self, mocker):
        """Check that no info is returned if any install is missing."""
        missing = ['one', 'cant-find-install']
        def mock_apt(_, **__):
            """Fake the result of calling 'apt-cache search'."""
            mock_proc = mocker.MagicMock()
            mock_proc.stdout = 'onelib - this is found'.encode()
            return mock_proc

        mocker.patch('subprocess.run', mock_apt)
        # pylint: disable-next=protected-access
        assert not QtDependencyFinder._list_install_for_libs_linux(missing)

    @parametrize(exc=(FileNotFoundError, SubprocessError))
    def test_list_install_libs_fails(self, exc, mocker):
        """Check result when subprocess raises."""
        missing = 'a', 'b'
        mocker.patch('subprocess.run', side_effect=exc)
        # pylint: disable-next=protected-access
        assert not QtDependencyFinder._list_install_for_libs_linux(missing)


class TestQtUnixPlatformGetter:
    """Tests for the _QtUnixPlatformsGetter helper class."""

    _supported = (
        '5.15',
        *(f'5.15.{v}' for v in range(18)),
        )
    _unsupported = (
        '4.12.99',
        '5.15.18',
        '5.16',
        '6.0.0',
        )

    @fixture(name='mock_qtc')
    def fixture_mock_qtc(self, mocker, mock_import_module):
        """Return a fake QtCore module."""
        mock_qtc = mocker.MagicMock()
        fake_modules = {'PyQt5.QtCore': mock_qtc}
        mock_import_module(**fake_modules)
        return mock_qtc

    @fixture
    @parametrize(version=_supported)
    def mock_supported_version(self, version, mock_qtc):
        """Fake a supported Qt version."""
        mock_qtc.QT_VERSION_STR = version

    def test_cant_import_qtcore(self, mocker):
        """Check outcome when QtCore fails to import."""
        mocker.patch('importlib.import_module', side_effect=ImportError)
        assert not _QtUnixPlatformsGetter.get(PYQT5)

    @parametrize(version=_unsupported)
    def test_unsupported_version(self, version, mock_qtc):
        """Check complaints for unsupported Qt versions."""
        mock_qtc.QT_VERSION_STR = version
        with pytest.raises(NotImplementedError):
            _QtUnixPlatformsGetter.get(PYQT5)

    def test_qt6(self, mocker, mock_import_module):
        """Check that Qt 6 is not supported yet, via PyQt6.QtCore."""
        mock_qtc = mocker.MagicMock(QT_VERSION_STR='6.8.0')
        mock_import_module(**{'PyQt6.QtCore': mock_qtc})
        with pytest.raises(NotImplementedError):
            _QtUnixPlatformsGetter.get(PYQT6)

    @pytest.mark.usefixtures('mock_supported_version')
    def test_qt_qpa_platform(self, mocker):
        """Check result when QT_QPA_PLATFORM is set."""
        mock_qpa = str(mocker.MagicMock())
        mocker.patch.dict('os.environ',
                          QT_QPA_PLATFORM=mock_qpa,
                          clear=True)
        result = _QtUnixPlatformsGetter.get(PYQT5)
        assert result == {mock_qpa}

    @pytest.mark.usefixtures('mock_supported_version')
    @parametrize(default=(None, 'some-default'))
    def test_no_xdg_session(self, default, mocker):
        """Check result when no XDG_SESSION_TYPE is set."""
        kwargs = {'clear': True}
        if default:
            kwargs['QT_QPA_DEFAULT_PLATFORM'] = default
        else:
            default = 'xcb'
        mocker.patch.dict('os.environ', **kwargs)
        result = _QtUnixPlatformsGetter.get(PYQT5)
        assert result == {default}

    @parametrize(default=(None, 'some-default'))
    @parametrize(gnome=({},
                        {'XDG_CURRENT_DESKTOP': 'some-gnome-value'},
                        {'XDG_SESSION_DESKTOP': 'some-gnome-value'}))
    @pytest.mark.usefixtures('mock_supported_version')
    def test_wayland(self, default, gnome, mocker):
        """Check result when XDG_SESSION_TYPE is set to wayland."""
        kwargs = {'clear': True,
                  'XDG_SESSION_TYPE': 'wayland',
                  **gnome}
        if default:
            kwargs['QT_QPA_DEFAULT_PLATFORM'] = default
        else:
            default = 'xcb'
        mocker.patch.dict('os.environ', **kwargs)
        result = _QtUnixPlatformsGetter.get(PYQT5)
        expect = {default} if gnome else {default, 'wayland'}
        assert result == expect

    @parametrize(default=(None, 'some-default'))
    @parametrize(xdg=('valid', 'invalid', 'xcb'))
    @pytest.mark.usefixtures('mock_supported_version')
    def test_xdg_other(self, default, xdg, mocker):
        """Check result when XDG_SESSION_TYPE is set to non-wayland/x11."""
        kwargs = {'clear': True,
                  'XDG_SESSION_TYPE': xdg}
        if default:
            kwargs['QT_QPA_DEFAULT_PLATFORM'] = default
        else:
            default = 'xcb'
        mocker.patch.dict('os.environ', **kwargs)
        result = _QtUnixPlatformsGetter.get(PYQT5)
        assert result == {xdg, default}

    @parametrize(default=(None, 'some-default'))
    @pytest.mark.usefixtures('mock_supported_version')
    def test_x_eleven(self, default, mocker):
        """Check result when XDG_SESSION_TYPE is set to x11."""
        kwargs = {'clear': True,
                  'XDG_SESSION_TYPE': 'x11'}
        if default:
            kwargs['QT_QPA_DEFAULT_PLATFORM'] = default
        else:
            default = 'xcb'
        mocker.patch.dict('os.environ', **kwargs)
        result = _QtUnixPlatformsGetter.get(PYQT5)
        assert result == {default, 'xcb'}
