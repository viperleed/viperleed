"""Tests for module binding of viperleed.gui.qt_compat."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-09'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.qt_compat import binding
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6
from viperleed.gui.qt_compat.binding import QT_API_ENV
from viperleed.gui.qt_compat.binding import QtApiError
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError
from viperleed.gui.qt_compat.binding import QtBindingTooOldError
from viperleed.gui.qt_compat.binding import check_version
from viperleed.gui.qt_compat.binding import get_binding

BOTH = (PYQT5, PYQT6)
_MODULE = 'viperleed.gui.qt_compat.binding'
_BROKEN = {
    'missing sip': ModuleNotFoundError("No module named 'PyQt6.sip'",
                                       name='PyQt6.sip'),
    'DLL load failed': ImportError('DLL load failed'),
    }
# PYQT_VERSION of usable bindings, and of too old ones
_NEW_ENOUGH = {PYQT5: 0x050F0A, PYQT6: 0x060C00}  # 5.15.10, 6.12.0
_OLDEST = {PYQT5: 0x050F00, PYQT6: 0x060100}      # 5.15.0, 6.1.0
_TOO_OLD = {PYQT5: 0x050E02, PYQT6: 0x060003}     # 5.14.2, 6.0.3


def _fake_qt_core(mocker, qt_binding, version):
    """Return a fake QtCore with PYQT_VERSION `version`."""
    version_str = '.'.join(str((version >> s) & 0xFF) for s in (16, 8, 0))
    return mocker.MagicMock(name=f'{qt_binding}.QtCore',
                            PYQT_VERSION=version,
                            PYQT_VERSION_STR=version_str)


@fixture(autouse=True)
def fixture_fresh_selection(monkeypatch):
    """Forget any previously selected binding, and ignore QT_API."""
    monkeypatch.setattr(binding, '_SELECTED_BINDING', None)
    monkeypatch.delenv(QT_API_ENV, raising=False)


@fixture(name='fake_bindings')
def fixture_fake_bindings(mocker):
    """Fake which Qt bindings can be imported/were already imported.

    The returned function takes the bindings that can be imported,
    those that were already imported, a dict {binding: exception}
    of bindings that are installed but raise exception on import,
    and the bindings that can be imported but are too old.
    """
    def _fake(available, imported=(), broken=None, too_old=()):
        broken = broken or {}
        def _import(name):
            package = name.split('.')[0]
            if package in too_old:
                return _fake_qt_core(mocker, package, _TOO_OLD[package])
            if package in available:
                return _fake_qt_core(mocker, package, _NEW_ENOUGH[package])
            if package in broken:
                raise broken[package]
            raise ModuleNotFoundError(f'No module named {package!r}',
                                      name=package)
        mocker.patch(f'{_MODULE}.sys',
                     modules={f'{b}.QtCore': None for b in imported})
        mock_importlib = mocker.patch(f'{_MODULE}.importlib')
        mock_importlib.import_module.side_effect = _import
        return mock_importlib.import_module
    return _fake


class TestCheckVersion:
    """Tests for the check_version function."""

    @parametrize(qt_binding=BOTH)
    def test_new_enough(self, qt_binding, mocker):
        """Check that recent bindings pass."""
        oldest = _fake_qt_core(mocker, qt_binding, _OLDEST[qt_binding])
        recent = _fake_qt_core(mocker, qt_binding, _NEW_ENOUGH[qt_binding])
        check_version(qt_binding, oldest)
        check_version(qt_binding, recent)

    @parametrize(qt_binding=BOTH)
    def test_too_old(self, qt_binding, mocker):
        """Check complaints for a binding that is too old."""
        qt_core = _fake_qt_core(mocker, qt_binding, _TOO_OLD[qt_binding])
        with pytest.raises(QtBindingTooOldError, match='too old') as exc:
            check_version(qt_binding, qt_core)
        assert qt_core.PYQT_VERSION_STR in str(exc.value)
        assert exc.value.qt_version is qt_core.QT_VERSION_STR
        assert not isinstance(exc.value, QtBindingNotFoundError)


class TestDefaultSelection:
    """Tests for selecting a binding without QT_API."""

    _available = {
        'both': (BOTH, PYQT6),
        'only PyQt5': ((PYQT5,), PYQT5),
        'only PyQt6': ((PYQT6,), PYQT6),
        }

    @parametrize(imported=BOTH)
    def test_already_imported(self, imported, fake_bindings):
        """Check that an already-imported binding is preferred."""
        fake_bindings(BOTH, imported=(imported,))
        assert get_binding() == imported

    @parametrize('available,expect', _available.values(), ids=_available)
    def test_available(self, available, expect, fake_bindings):
        """Check that the preferred available binding is selected."""
        fake_bindings(available)
        assert get_binding() == expect

    @parametrize(error=_BROKEN.values(), ids=_BROKEN)
    def test_broken_none_available(self, error, fake_bindings):
        """Check that a broken binding is not reported as missing."""
        fake_bindings((), broken={PYQT6: error})
        with pytest.raises(ImportError, match=str(error)) as exc_info:
            get_binding()
        assert not isinstance(exc_info.value, QtBindingNotFoundError)

    @parametrize(error=_BROKEN.values(), ids=_BROKEN)
    def test_broken_skipped(self, error, fake_bindings):
        """Check that another binding is used if one is broken."""
        fake_bindings((PYQT5,), broken={PYQT6: error})
        assert get_binding() == PYQT5

    def test_fall_back_to_pyqt5(self, fake_bindings, mocker):
        """Check that PyQt5 is tried if PyQt6 cannot be imported."""
        mock_import = fake_bindings((PYQT5,))
        get_binding()
        assert mock_import.mock_calls == [mocker.call('PyQt6'),
                                          mocker.call('PyQt5'),
                                          mocker.call('PyQt5.QtCore')]

    def test_namespace_package_left(self, fake_bindings):
        """Check that a binding without QtCore counts as missing."""
        missing = ModuleNotFoundError("No module named 'PyQt6.QtCore'",
                                      name='PyQt6.QtCore')
        fake_bindings((), broken={PYQT6: missing})
        with pytest.raises(QtBindingNotFoundError):
            get_binding()

    def test_none_available(self, fake_bindings):
        """Check complaints when no binding can be imported."""
        fake_bindings(())
        with pytest.raises(QtBindingNotFoundError) as exc_info:
            get_binding()
        assert isinstance(exc_info.value, ImportError)
        msg = str(exc_info.value)
        assert all(b in msg for b in BOTH)
        assert QT_API_ENV not in msg

    def test_selection_cached(self, fake_bindings):
        """Check that a binding is selected only once."""
        mock_import = fake_bindings(BOTH)
        first = get_binding()
        mock_import.reset_mock()
        assert get_binding() == first
        mock_import.assert_not_called()

    def test_too_old_none_available(self, fake_bindings):
        """Check that a too-old binding is not reported as missing."""
        fake_bindings((), too_old=(PYQT6,))
        with pytest.raises(ImportError, match='too old') as exc_info:
            get_binding()
        assert not isinstance(exc_info.value, QtBindingNotFoundError)

    def test_too_old_skipped(self, fake_bindings):
        """Check that another binding is used if one is too old."""
        fake_bindings((PYQT5,), too_old=(PYQT6,))
        assert get_binding() == PYQT5


class TestQtApi:
    """Tests for selecting a binding via the QT_API variable."""

    _requested = {
        'pyqt5': PYQT5,
        'PyQt5': PYQT5,
        ' PYQT5 ': PYQT5,
        'pyqt6': PYQT6,
        'PyQt6': PYQT6,
        }

    @parametrize(qt_api=('', '   '))
    def test_empty_ignored(self, qt_api, fake_bindings, monkeypatch):
        """Check that an empty QT_API is treated as if it were unset."""
        monkeypatch.setenv(QT_API_ENV, qt_api)
        fake_bindings(BOTH)
        assert get_binding() == PYQT6

    def test_other_imported(self, fake_bindings, monkeypatch):
        """Check complaints if QT_API requests another binding."""
        monkeypatch.setenv(QT_API_ENV, 'pyqt5')
        mock_import = fake_bindings(BOTH, imported=(PYQT6,))
        with pytest.raises(QtApiError, match='already imported'):
            get_binding()
        mock_import.assert_not_called()

    @parametrize('qt_api,expect', _requested.items(), ids=_requested)
    def test_requested(self, qt_api, expect, fake_bindings, monkeypatch):
        """Check that the binding in QT_API is used."""
        monkeypatch.setenv(QT_API_ENV, qt_api)
        fake_bindings(BOTH)
        assert get_binding() == expect

    def test_requested_missing(self, fake_bindings, monkeypatch):
        """Check complaints when the requested binding is missing."""
        monkeypatch.setenv(QT_API_ENV, 'pyqt5')
        mock_import = fake_bindings((PYQT6,))
        with pytest.raises(QtApiError, match=QT_API_ENV):
            get_binding()
        mock_import.assert_called_once_with('PyQt5')

    def test_requested_overrides_preference(self, fake_bindings,
                                            monkeypatch, mocker):
        """Check that QT_API wins over the preferred binding."""
        monkeypatch.setenv(QT_API_ENV, 'pyqt5')
        mock_import = fake_bindings(BOTH)
        assert get_binding() == PYQT5
        assert mock_import.mock_calls == [mocker.call('PyQt5'),
                                          mocker.call('PyQt5.QtCore')]

    @parametrize(qt_api=('pyside6', 'pyqt4', 'qt5'))
    def test_unsupported(self, qt_api, fake_bindings, monkeypatch):
        """Check complaints for an unsupported QT_API value."""
        monkeypatch.setenv(QT_API_ENV, qt_api)
        mock_import = fake_bindings(BOTH)
        with pytest.raises(QtApiError, match=qt_api):
            get_binding()
        mock_import.assert_not_called()


class TestGetCandidates:
    """Tests for the get_candidates function."""

    @parametrize(imported=BOTH)
    def test_already_imported(self, imported, fake_bindings):
        """Check that an already imported binding is the only one."""
        mock_import = fake_bindings((), imported=(imported,))
        assert binding.get_candidates() == (imported,)
        mock_import.assert_not_called()

    def test_default(self, fake_bindings):
        """Check that both bindings are candidates, in order."""
        mock_import = fake_bindings(())
        assert binding.get_candidates() == (PYQT6, PYQT5)
        mock_import.assert_not_called()

    @parametrize(imported=((), (PYQT5,)))
    def test_requested(self, imported, fake_bindings, monkeypatch):
        """Check that a binding requested via QT_API is the only one."""
        monkeypatch.setenv(QT_API_ENV, 'pyqt5')
        fake_bindings((), imported=imported)
        assert binding.get_candidates() == (PYQT5,)

    def test_requested_other_imported(self, fake_bindings, monkeypatch):
        """Check complaints if QT_API requests another binding."""
        monkeypatch.setenv(QT_API_ENV, 'pyqt5')
        fake_bindings((), imported=(PYQT6,))
        with pytest.raises(QtApiError, match='already imported'):
            binding.get_candidates()

    def test_unsupported(self, monkeypatch):
        """Check complaints for an unsupported QT_API value."""
        monkeypatch.setenv(QT_API_ENV, 'pyside6')
        with pytest.raises(QtApiError):
            binding.get_candidates()


class TestImportBinding:
    """Tests for the import_binding function."""

    _widgets = f'{PYQT6}.QtWidgets'

    @fixture(name='mock_import')
    def fixture_mock_import(self, mocker):
        """Return a fake import_module, for a recent binding."""
        def _import(name):
            package = name.split('.')[0]
            return _fake_qt_core(mocker, package, _NEW_ENOUGH[package])
        mock_importlib = mocker.patch(f'{_MODULE}.importlib')
        mock_importlib.import_module.side_effect = _import
        return mock_importlib.import_module

    @parametrize(error=(ModuleNotFoundError(name=_widgets),
                        ModuleNotFoundError(name=f'{PYQT6}.sip'),
                        ImportError('DLL load failed')))
    def test_broken(self, error, mock_import):
        """Check that other import errors mean a broken binding."""
        default = mock_import.side_effect
        def _import(name):
            if name == self._widgets:
                raise error
            return default(name)
        mock_import.side_effect = _import
        with pytest.raises(ImportError) as exc_info:
            binding.import_binding(PYQT6, ('QtWidgets',))
        assert exc_info.value is error
        assert not isinstance(exc_info.value, QtBindingNotFoundError)

    @parametrize(qt_binding=BOTH)
    def test_imports_once_in_order(self, qt_binding, mock_import, mocker):
        """Check that each module is imported once, QtCore first."""
        binding.import_binding(qt_binding, ('QtWidgets', 'QtCore', 'QtGui'))
        assert mock_import.mock_calls == [
            mocker.call(qt_binding),
            mocker.call(f'{qt_binding}.QtCore'),
            mocker.call(f'{qt_binding}.QtWidgets'),
            mocker.call(f'{qt_binding}.QtGui'),
            ]

    @parametrize(missing=('', '.QtCore'))
    def test_not_installed(self, missing, mock_import):
        """Check that a missing package or QtCore means 'missing'."""
        missing = PYQT6 + missing
        mock_import.side_effect = ModuleNotFoundError(name=missing)
        with pytest.raises(QtBindingNotFoundError):
            binding.import_binding(PYQT6, ('QtWidgets',))

    @parametrize(qt_binding=BOTH)
    def test_returns_qt_core(self, qt_binding, mock_import):
        """Check that the QtCore module imported is returned."""
        imported = {}
        default = mock_import.side_effect
        def _import(name):
            imported[name] = default(name)
            return imported[name]
        mock_import.side_effect = _import
        qt_core = binding.import_binding(qt_binding, ('QtWidgets',))
        assert qt_core is imported[f'{qt_binding}.QtCore']

    def test_too_old(self, mock_import, mocker):
        """Check that nothing else is imported for an old binding."""
        mock_import.side_effect = (
            lambda name: _fake_qt_core(mocker, PYQT6, _TOO_OLD[PYQT6])
            )
        with pytest.raises(QtBindingTooOldError):
            binding.import_binding(PYQT6, ('QtWidgets',))
        assert mock_import.mock_calls[-1] == mocker.call(f'{PYQT6}.QtCore')


class TestSetBinding:
    """Tests for the set_binding function."""

    def test_different(self):
        """Check complaints when another binding was selected."""
        binding.set_binding(PYQT5)
        with pytest.raises(ValueError, match='already in use'):
            binding.set_binding(PYQT6)

    def test_same_again(self):
        """Check that selecting the same binding again is fine."""
        binding.set_binding(PYQT5)
        binding.set_binding(PYQT5)
        assert get_binding() == PYQT5

    @parametrize(selected=BOTH)
    def test_selects(self, selected, fake_bindings):
        """Check that the binding is selected without importing it."""
        mock_import = fake_bindings(())
        binding.set_binding(selected)
        assert get_binding() == selected
        mock_import.assert_not_called()

    @parametrize(unsupported=('PySide6', 'pyqt5', None))
    def test_unsupported(self, unsupported):
        """Check complaints for an unsupported binding."""
        with pytest.raises(ValueError, match='Unsupported'):
            binding.set_binding(unsupported)
