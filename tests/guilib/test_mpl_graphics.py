"""Tests for module mpl_graphics of viperleed.guilib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-26'  # Parts of this were in basewidgets
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.guilib.mpl_graphics import MatplotLibBackend
from viperleed.guilib.mpl_graphics import NO_CAIRO
from viperleed.guilib.mpl_graphics import find_and_import_backend
from viperleed.guilib.mpl_graphics import get_cairo_version
from viperleed.guilib.mpl_graphics import import_figure_canvas
from viperleed.guilib.mpl_graphics import import_matplotlib


_MODULE = 'viperleed.guilib.mpl_graphics'

PY37 = (3, 7)
PY38 = (3, 8)
PY310 = (3, 10)


@fixture(name='mock_pyqt')
def fixture_mock_pyqt(mocker):
    """Fake the presence/absence of PyQt."""
    def _mock(present):
        mocker.patch(f'{_MODULE}.has_pyqt', return_value=present)
    return _mock


class TestMatplotLibBackend:
    """Tests for the MatplotLibBackend enumeration."""

    def test_enum_values(self):
        """Check the enumeration values."""
        # pylint: disable=magic-value-comparison
        assert MatplotLibBackend.AGG.module is None
        assert MatplotLibBackend.AGG.mpl_use == 'qt5agg'
        assert MatplotLibBackend.CAIRO.module == 'mplcairo'
        assert MatplotLibBackend.CAIRO.mpl_use == 'module://mplcairo.qt'
        assert MatplotLibBackend.DEFAULT is MatplotLibBackend.CAIRO


# TODO: these tests will need adaptation when we'll fix mplcairo issues
class TestFindAndImportBackend:
    """Tests for the find_and_import_backend function."""

    _cairo_version = {
        (1, 17, 2): MatplotLibBackend.AGG,
        (1, 17, 9999): MatplotLibBackend.AGG,
        (1, 18, 0): MatplotLibBackend.CAIRO,
        (1, 18, 2): MatplotLibBackend.CAIRO,
        (1, 19, 0): MatplotLibBackend.CAIRO,
        }

    @parametrize('cairo_version,expect', _cairo_version.items())
    def test_import_success(self, cairo_version, expect, mocker):
        """Check the outcome of a successful backend import."""
        mocker.patch('importlib.import_module')
        mocker.patch(f'{_MODULE}.get_cairo_version',
                     return_value=cairo_version)
        backend = find_and_import_backend()
        assert backend is expect

    def test_no_module_to_import(self, mocker):
        """Check result when the backend needs no import."""
        mock_backend = mocker.patch(f'{_MODULE}.MatplotLibBackend')
        mock_backend.DEFAULT.module = None
        backend = find_and_import_backend()
        assert backend is mock_backend.DEFAULT

    def test_import_failure(self, mocker):
        """Check fallback when importing a backend module fails."""
        mock_backend = mocker.patch(f'{_MODULE}.MatplotLibBackend')
        mock_backend.DEFAULT.module = 'fake_module'
        mock_backend.DEFAULT.mpl_use = 'fake_use'
        mocker.patch('importlib.import_module', side_effect=ImportError)
        backend = find_and_import_backend()
        assert backend is mock_backend.AGG


class TestGetCairoVersion:
    """Tests for the get_cairo_version function."""

    def test_import_fails(self, mocker):
        """Check expected version when importing mplcairo fails."""
        mocker.patch('importlib.import_module', side_effect=ImportError)
        result = get_cairo_version()
        assert result == NO_CAIRO

    def test_invalid_version(self, mocker):
        """Check result when parsing a version string fails."""
        mock_cairo = mocker.MagicMock()
        mock_cairo.get_versions.return_value = {'cairo': 'a.b.c'}
        mocker.patch('importlib.import_module', return_value=mock_cairo)
        result = get_cairo_version()
        assert result == NO_CAIRO

    def test_no_version_info(self, mocker):
        """Check expected version when the version information is missing."""
        mock_cairo = mocker.MagicMock()
        mock_cairo.get_versions.return_value = {'other_key': 'other_value'}
        mocker.patch('importlib.import_module', return_value=mock_cairo)
        result = get_cairo_version()
        assert result == NO_CAIRO

    _valid = {
        'no at': ('1.6.35', (1, 6, 35)),
        'at': ('3.2.1 @ location', (3, 2, 1)),
        }

    @parametrize('version_str,expect', _valid.values(), ids=_valid)
    def test_valid_version(self, version_str, expect, mocker):
        """Check expected outcome with valid version information."""
        mock_cairo = mocker.MagicMock()
        mock_cairo.get_versions.return_value = {'cairo': version_str}
        mocker.patch('importlib.import_module', return_value=mock_cairo)
        result = get_cairo_version()
        assert result == expect


class TestImportMatplotlib:
    """Tests for the import_matplotlib function."""

    def test_no_pyqt(self, mock_pyqt):
        """Check that no import is done without PyQt."""
        mock_pyqt(False)
        assert import_matplotlib() is None

    def test_with_pyqt(self, mock_pyqt, mocker):
        """Check implementation when PyQt is present."""
        mock_pyqt(True)
        mock_backend = mocker.MagicMock()
        mock_matplotlib = mocker.MagicMock()
        mocker.patch(f'{_MODULE}.find_and_import_backend',
                     return_value=mock_backend)
        mock_import = mocker.patch('importlib.import_module',
                                   return_value=mock_matplotlib)
        result = import_matplotlib()
        assert result is mock_backend
        mock_import.assert_called_once_with('matplotlib')
        mock_matplotlib.use.assert_called_once_with(mock_backend.mpl_use)


class TestImportFigureCanvas:
    """Tests for the import_figure_canvas function."""

    def test_no_pyqt_raises(self, mock_pyqt):
        """Check complaints when called without PyQt."""
        mock_pyqt(False)
        with pytest.raises(ImportError, match='PyQt5 not found'):
            import_figure_canvas()

    _canvas_for_backend = {
        MatplotLibBackend.AGG: (
            'matplotlib.backends.backend_qt5agg.FigureCanvas'
            ),
        MatplotLibBackend.CAIRO: 'mplcairo.qt.FigureCanvasQTCairo',
        }

    @parametrize('backend,imported', _canvas_for_backend.items())
    def test_valid_backend(self, backend, imported, mock_pyqt, mocker):
        """Check the expected result when importing a valid backend."""
        mock_pyqt(True)
        mock_module = mocker.MagicMock()
        mocker.patch(f'{_MODULE}.import_matplotlib', return_value=backend)
        mock_import = mocker.patch('importlib.import_module',
                                   return_value=mock_module)
        result = import_figure_canvas()
        module_name, canvas_name = imported.rsplit('.', maxsplit=1)
        assert result is getattr(mock_module, canvas_name)
        mock_import.assert_called_once_with(module_name)

    def test_invalid_backend(self, mock_pyqt, mocker):
        """Check complaints when an invalid backend is selected."""
        mock_pyqt(True)
        mocker.patch(f'{_MODULE}.import_matplotlib',
                     return_value=mocker.MagicMock())
        with pytest.raises(ImportError, match='Unsupported backend'):
            import_figure_canvas()
