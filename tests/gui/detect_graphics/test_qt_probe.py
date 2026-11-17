"""Tests for the _QtProbe class of detect_graphics."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'

import pickle

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.detect_graphics import _EXIT_CODE_OFFSET
from viperleed.gui.detect_graphics import _QtProbe
from viperleed.gui.detect_graphics import PyQtSanity
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6

from .conftest import BOTH
from .conftest import MODULE


class TestQtProbe:
    """Tests for calling a _QtProbe, as in the checker subprocess."""

    # The modules that we attempt to import, besides the package
    _modules = (
        'QtCore',
        'QtWidgets',
        )

    @fixture(name='mock_dummy_qapp', autouse=True)
    def fixture_mock_dummy_qapp(self, mocker):
        """Replace _QtProbe._init_dummy_qapp with a mock."""
        # This is to ensure tests don't hang forever on a system
        # without graphics.
        return mocker.patch(f'{MODULE}._QtProbe._init_dummy_qapp')

    @fixture(name='imports_done')
    def fixture_imports_done(self, mocker):
        """Return a fake sending end of a pipe."""
        return mocker.MagicMock(name='imports_done')

    @staticmethod
    def check_exit_code(excinfo, expect):
        """Check that the exit code of `excinfo` is for `expect`."""
        assert excinfo.value.code == _EXIT_CODE_OFFSET + expect.value

    def test_candidates_pickled(self):
        """Check that a probe survives pickling, as when spawned."""
        probe = pickle.loads(pickle.dumps(_QtProbe(list(BOTH))))
        # pylint: disable-next=protected-access
        assert probe._candidates == BOTH

    def test_cant_import_pyqt(self, imports_done, mocker):
        """Check the exit code when importing any binding fails."""
        mock_import = mocker.patch('importlib.import_module',
                                   side_effect=ImportError)
        with pytest.raises(SystemExit) as excinfo:
            _QtProbe(BOTH)(imports_done)
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
            lambda mod: ModuleNotFoundError(name=mod.split('.')[0] + '.sip'),
            PyQtSanity.IMPORT_ERROR,
            ),
        }

    @parametrize(binding=BOTH)
    @parametrize(module=_modules)
    @parametrize('make_exc,expect', _submodule_errors.values(),
                 ids=_submodule_errors)
    # pylint: disable-next=too-many-arguments,too-many-positional-arguments
    def test_cant_import_pyqt_submodule(self, binding, module, make_exc,
                                        expect, imports_done, mocker):
        """Check the exit code when importing a submodule fails."""
        failing = f'{binding}.{module}'
        def _fail_module_import(mod_name):
            if mod_name == failing:
                raise make_exc(mod_name)
        mocker.patch('importlib.import_module',
                     side_effect=_fail_module_import)
        with pytest.raises(SystemExit) as excinfo:
            _QtProbe((binding,))(imports_done)
        self.check_exit_code(excinfo, expect)
        imports_done.send.assert_not_called()

    def test_pyqt_not_found(self, imports_done, mocker):
        """Check the expected exit code when no binding is found."""
        mocker.patch('importlib.import_module',
                     side_effect=ModuleNotFoundError)
        with pytest.raises(SystemExit) as excinfo:
            _QtProbe(BOTH)(imports_done)
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
            _QtProbe(BOTH)(imports_done)
        self.check_exit_code(excinfo, PyQtSanity.IMPORT_ERROR)

    @parametrize(binding=BOTH)
    def test_success(self, binding, imports_done, mocker, mock_dummy_qapp):
        """Test an execution that does not make the interpreter exit."""
        # NB: normally, the interpreter always exits when
        # a _QtProbe is called. The only
        # reason why we don't exit here is the patching of
        # _init_dummy_qapp via mock_dummy_qapp.
        mock_import = mocker.patch('importlib.import_module')
        calls = mocker.MagicMock()
        calls.attach_mock(imports_done, 'imports_done')
        calls.attach_mock(mock_dummy_qapp, 'init_dummy_qapp')
        _QtProbe((binding,))(imports_done)

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
        _QtProbe(BOTH)(imports_done)
        imports_done.send.assert_called_once_with(PYQT5)
        mock_dummy_qapp.assert_called_once_with(PYQT5)


class TestQtProbeDummyQapp:
    """Tests for the _init_dummy_qapp method of _QtProbe."""

    @fixture(name='binding')
    @parametrize(binding=BOTH)
    def fixture_binding(self, binding):
        """Return the name of a Qt binding."""
        return binding

    @fixture(name='mock_qapp')
    def fixture_mock_qapp(self, binding, mocker, mock_import_module):
        """Replace QApplication of `binding` with a mock."""
        mocker.patch(f'{MODULE}.suppress_file_permission_warnings')
        # Fake only QtWidgets: newer Pythons (e.g., 3.13) find the
        # targets of mock.patch via importlib.import_module, and later
        # patches (e.g., in test_calls) must reach the real modules.
        mock_qtwidgets = mocker.MagicMock()
        mock_import_module(**{f'{binding}.QtWidgets': mock_qtwidgets})
        return mock_qtwidgets.QApplication

    def test_calls(self, binding, mock_qapp, mocker):
        """Check calls in _init_dummy_qapp."""
        mock_suppress = mocker.patch(
            f'{MODULE}.suppress_file_permission_warnings'
            )
        with pytest.raises(SystemExit):
            # pylint: disable-next=protected-access
            _QtProbe._init_dummy_qapp(binding)
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
            # pylint: disable-next=protected-access
            _QtProbe._init_dummy_qapp(binding)
        assert excinfo.value.code == _EXIT_CODE_OFFSET + expect.value
