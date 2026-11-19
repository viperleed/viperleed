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
from viperleed.gui.detect_graphics import _QT_MODULES
from viperleed.gui.detect_graphics import _QtProbe
from viperleed.gui.detect_graphics import PyQtSanity
from viperleed.gui.qt_compat.binding import QtBindingTooOldError

from .conftest import BOTH
from .conftest import MODULE


@parametrize(binding=BOTH)
class TestQtProbe:
    """Tests for calling a _QtProbe, as in the checker subprocess."""

    # The modules that we attempt to import, besides the package
    _modules = _QT_MODULES
    # Exceptions raised when importing a submodule of a binding, as a
    # function of the submodule's name. All mean a broken binding.
    _submodule_errors = {
        'import error': lambda _: ImportError(),
        'sip mismatch': lambda _: RuntimeError('the sip module implements '
                                               'API v12.8 but ...'),
        # A dependency of the submodule is missing: broken PyQt5
        'dependency missing': (
            lambda mod: ModuleNotFoundError(name=mod.split('.')[0] + '.sip')
            ),
        }
    _qt_version = '6.12.0'

    @fixture(name='mock_dummy_qapp', autouse=True)
    def fixture_mock_dummy_qapp(self, mocker):
        """Replace _QtProbe._init_dummy_qapp with a mock."""
        # This is to ensure tests don't hang forever on a system
        # without graphics.
        return mocker.patch(f'{MODULE}._QtProbe._init_dummy_qapp')

    @fixture(name='check_version', autouse=True)
    def fixture_check_version(self, mocker):
        """Take any binding as new enough, unless told otherwise."""
        return mocker.patch(
            'viperleed.gui.qt_compat.binding.check_version'
            )

    @fixture(name='imports_done')
    def fixture_imports_done(self, mocker):
        """Return a fake sending end of a pipe."""
        return mocker.MagicMock(name='imports_done')

    @fixture(name='mock_import')
    def fixture_mock_import(self, mocker):
        """Return a function that fakes importing modules.

        The function takes the name of a module that fails to import,
        and a function that makes the exception it raises from the
        name of the module. It returns the fake import_module.
        """
        def _mock(failing=None, make_exc=None):
            def _import(name):
                if name == failing:
                    raise make_exc(name)
                if name.endswith('.QtCore'):
                    return mocker.MagicMock(QT_VERSION_STR=self._qt_version)
                return mocker.MagicMock()
            return mocker.patch('importlib.import_module',
                                side_effect=_import)
        return _mock

    @staticmethod
    def check_failed(binding, imports_done, expect, qt_version=None):
        """Check that probing `binding` exits with `expect`.

        Parameters
        ----------
        binding : str
            The binding to probe.
        imports_done : MagicMock
            The fake sending end of the pipe.
        expect : PyQtSanity
            The expected exit status.
        qt_version : str or None, optional
            The Qt version expected to be sent. Default is None, i.e.,
            if QtCore could not be imported.
        """
        with pytest.raises(SystemExit) as excinfo:
            _QtProbe(binding)(imports_done)
        assert excinfo.value.code == _EXIT_CODE_OFFSET + expect.value
        # The parent is told when importing is over, also on failure
        imports_done.send.assert_called_once_with(qt_version)
        imports_done.close.assert_called_once_with()

    def test_binding_pickled(self, binding):
        """Check that a probe survives pickling, as when spawned."""
        probe = pickle.loads(pickle.dumps(_QtProbe(binding)))
        # pylint: disable-next=protected-access
        assert probe._binding == binding

    def test_cant_import_pyqt(self, binding, imports_done, mock_import):
        """Check the exit code when importing the binding fails."""
        mock_import(binding, lambda _: ImportError())
        self.check_failed(binding, imports_done, PyQtSanity.IMPORT_ERROR)

    @parametrize(module=_modules)
    @parametrize(make_exc=_submodule_errors.values(), ids=_submodule_errors)
    def test_cant_import_pyqt_submodule(self, binding, module, make_exc,
                                        imports_done, mock_import):
        """Check the exit code when importing a submodule fails."""
        mock_import(f'{binding}.{module}', make_exc)
        # The Qt version is known if QtCore could be imported
        # pylint: disable-next=magic-value-comparison
        qt_version = None if module == 'QtCore' else self._qt_version
        self.check_failed(binding, imports_done, PyQtSanity.IMPORT_ERROR,
                          qt_version)

    @parametrize(module=_modules)
    def test_missing_pyqt_submodule(self, binding, module,
                                    imports_done, mock_import):
        """Check that only a missing QtCore means 'not installed'."""
        # E.g., only PyQt5-sip and PyQt5-Qt5 are left after
        # uninstalling PyQt5: 'PyQt5' is a namespace package
        mock_import(f'{binding}.{module}',
                    lambda name: ModuleNotFoundError(name=name))
        # pylint: disable-next=magic-value-comparison
        if module == 'QtCore':
            self.check_failed(binding, imports_done, PyQtSanity.NOT_FOUND)
        else:
            self.check_failed(binding, imports_done,
                              PyQtSanity.IMPORT_ERROR, self._qt_version)

    def test_pyqt_dependency_missing(self, binding, imports_done,
                                     mock_import):
        """Check that a package missing a dependency is broken."""
        mock_import(binding,
                    lambda _: ModuleNotFoundError(name=f'{binding}.sip'))
        self.check_failed(binding, imports_done, PyQtSanity.IMPORT_ERROR)

    def test_pyqt_not_found(self, binding, imports_done, mock_import):
        """Check the expected exit code when the binding is missing."""
        mock_import(binding, lambda name: ModuleNotFoundError(name=name))
        self.check_failed(binding, imports_done, PyQtSanity.NOT_FOUND)

    def test_success(self, binding, imports_done, mock_import, mocker,
                     mock_dummy_qapp):
        """Test an execution that does not make the interpreter exit."""
        # NB: normally, the interpreter always exits when
        # a _QtProbe is called. The only
        # reason why we don't exit here is the patching of
        # _init_dummy_qapp via mock_dummy_qapp.
        fake_import = mock_import()
        calls = mocker.MagicMock()
        calls.attach_mock(imports_done, 'imports_done')
        calls.attach_mock(mock_dummy_qapp, 'init_dummy_qapp')
        _QtProbe(binding)(imports_done)

        imported = {c[0][0] for c in fake_import.call_args_list}
        assert imported == {binding,
                            *(f'{binding}.{m}' for m in self._modules)}
        # The parent must be told before creating the QApplication.
        # Build the expected calls outside the assert: code-quality
        # checkers take mocker.call.<...>.send() for a real call.
        expect_order = [
            mocker.call.imports_done.send(self._qt_version),
            mocker.call.imports_done.close(),
            mocker.call.init_dummy_qapp(binding),
            ]
        assert calls.mock_calls == expect_order

    def test_too_old(self, binding, check_version, imports_done,
                     mock_import):
        """Check the exit code when the binding is too old."""
        fake_import = mock_import()
        check_version.side_effect = QtBindingTooOldError(
            'too old',
            qt_version=self._qt_version,
            )
        # The Qt version is known also for an old binding
        self.check_failed(binding, imports_done, PyQtSanity.TOO_OLD,
                          self._qt_version)
        # Nothing is imported after QtCore
        assert fake_import.call_args_list[-1][0] == (f'{binding}.QtCore',)


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
