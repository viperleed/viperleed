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
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError
from viperleed.gui.qt_compat.binding import get_binding

_MODULE = 'viperleed.gui.qt_compat.binding'
BOTH = (PYQT5, PYQT6)


@fixture(autouse=True)
def fixture_fresh_selection(monkeypatch):
    """Forget any previously selected binding, and ignore QT_API."""
    monkeypatch.setattr(binding, '_SELECTED_BINDING', None)
    monkeypatch.delenv(QT_API_ENV, raising=False)


@fixture(name='fake_bindings')
def fixture_fake_bindings(mocker):
    """Fake which Qt bindings can be imported/were already imported."""
    def _fake(available, imported=()):
        def _import(name):
            if name.split('.')[0] in available:
                return mocker.MagicMock()
            raise ImportError(f'No module named {name!r}')
        mocker.patch(f'{_MODULE}.sys',
                     modules={f'{b}.QtCore': None for b in imported})
        mock_importlib = mocker.patch(f'{_MODULE}.importlib')
        mock_importlib.import_module.side_effect = _import
        return mock_importlib.import_module
    return _fake


class TestDefaultSelection:
    """Tests for selecting a binding without QT_API."""

    _available = {
        'both': (BOTH, PYQT6),
        'only PyQt5': ((PYQT5,), PYQT5),
        'only PyQt6': ((PYQT6,), PYQT6),
        }

    @parametrize('available,expect', _available.values(), ids=_available)
    def test_available(self, available, expect, fake_bindings):
        """Check that the preferred available binding is selected."""
        fake_bindings(available)
        assert get_binding() == expect

    def test_fall_back_to_pyqt5(self, fake_bindings, mocker):
        """Check that PyQt5 is tried if PyQt6 cannot be imported."""
        mock_import = fake_bindings((PYQT5,))
        get_binding()
        assert mock_import.mock_calls == [mocker.call('PyQt6.QtCore'),
                                          mocker.call('PyQt5.QtCore')]

    @parametrize(imported=BOTH)
    def test_already_imported(self, imported, fake_bindings):
        """Check that an already-imported binding is preferred."""
        fake_bindings(BOTH, imported=(imported,))
        assert get_binding() == imported

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


class TestQtApi:
    """Tests for selecting a binding via the QT_API variable."""

    _requested = {
        'pyqt5': PYQT5,
        'PyQt5': PYQT5,
        ' PYQT5 ': PYQT5,
        'pyqt6': PYQT6,
        'PyQt6': PYQT6,
        }

    @parametrize('qt_api,expect', _requested.items(), ids=_requested)
    def test_requested(self, qt_api, expect, fake_bindings, monkeypatch):
        """Check that the binding in QT_API is used."""
        monkeypatch.setenv(QT_API_ENV, qt_api)
        fake_bindings(BOTH)
        assert get_binding() == expect

    def test_requested_overrides_preference(self, fake_bindings,
                                            monkeypatch):
        """Check that QT_API wins over the preferred binding."""
        monkeypatch.setenv(QT_API_ENV, 'pyqt5')
        mock_import = fake_bindings(BOTH)
        assert get_binding() == PYQT5
        mock_import.assert_called_once_with('PyQt5.QtCore')

    def test_requested_missing(self, fake_bindings, monkeypatch):
        """Check complaints when the requested binding is missing."""
        monkeypatch.setenv(QT_API_ENV, 'pyqt5')
        mock_import = fake_bindings((PYQT6,))
        with pytest.raises(QtBindingNotFoundError, match=QT_API_ENV):
            get_binding()
        mock_import.assert_called_once_with('PyQt5.QtCore')

    @parametrize(qt_api=('pyside6', 'pyqt4', 'qt5'))
    def test_unsupported(self, qt_api, fake_bindings, monkeypatch):
        """Check complaints for an unsupported QT_API value."""
        monkeypatch.setenv(QT_API_ENV, qt_api)
        mock_import = fake_bindings(BOTH)
        with pytest.raises(QtBindingNotFoundError, match=qt_api):
            get_binding()
        mock_import.assert_not_called()

    @parametrize(qt_api=('', '   '))
    def test_empty_ignored(self, qt_api, fake_bindings, monkeypatch):
        """Check that an empty QT_API is treated as if it were unset."""
        monkeypatch.setenv(QT_API_ENV, qt_api)
        fake_bindings(BOTH)
        assert get_binding() == PYQT6
