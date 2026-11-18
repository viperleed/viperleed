"""Tests for the public functions of module detect_graphics."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-17'
__license__ = 'GPLv3+'

import os
from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.gui.detect_graphics import PyQtSanity
from viperleed.gui.detect_graphics import find_missing_qt_dependencies
from viperleed.gui.detect_graphics import has_graphics
from viperleed.gui.detect_graphics import has_pyqt
from viperleed.gui.detect_graphics import suppress_file_permission_warnings
from viperleed.gui.qt_compat.binding import PYQT5
from viperleed.gui.qt_compat.binding import PYQT6
from viperleed.gui.qt_compat.binding import QtBindingNotFoundError

from .conftest import MODULE


@parametrize(pyqt_present=(True, False))
@parametrize(has_display=(True, False))
def test_has_graphics(pyqt_present, has_display, mocker):
    """Ensure correct identification of graphics capabilities."""
    mocker.patch(f'{MODULE}.has_pyqt', return_value=pyqt_present)
    check_sanity = mocker.patch(
        f'{MODULE}.check_pyqt_sanity',
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
    mocker.patch(f'{MODULE}.check_pyqt_sanity', return_value=pyqt_sanity)
    pyqt_present = pyqt_sanity in {PyQtSanity.OK, PyQtSanity.NO_DISPLAY}
    assert has_pyqt() == pyqt_present


class TestFindMissingDependencies:
    """Tests for the find_missing_qt_dependencies function."""

    @fixture(autouse=True)
    def mock_sane_pyqt(self, mocker):
        """Do not start a real QApplication to check Qt."""
        # Unmocked, the first call would run a subprocess that starts
        # a QApplication, which aborts (and makes pytest's faulthandler
        # print "Fatal Python error") on systems without a display.
        mocker.patch(f'{MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.OK,
                     binding=PYQT5)

    @parametrize(sanity=(PyQtSanity.NOT_FOUND,
                         PyQtSanity.INVALID_QT_API,
                         PyQtSanity.TOO_OLD))
    def test_raises(self, sanity, mocker):
        """Check complaints when calling the function without Qt."""
        mocker.patch(f'{MODULE}.check_pyqt_sanity', return_value=sanity)
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
        mocker.patch(f'{MODULE}.QtDependencyFinder',
                     return_value=mock_finder)
        mock_finder.find_missing_dependencies.side_effect = NotImplementedError
        missing = find_missing_qt_dependencies()
        assert not missing

    def test_missing_deps(self, mocker):
        """Check the expected result when missing dependencies are found."""
        mocker.patch('sys.platform', 'some-supported-platform')
        mock_finder = mocker.MagicMock()
        mock_finder_cls = mocker.patch(f'{MODULE}.QtDependencyFinder',
                                       return_value=mock_finder)
        missing = find_missing_qt_dependencies()
        assert missing is mock_finder.find_missing_dependencies.return_value
        mock_finder_cls.assert_called_once_with(PYQT5)

    def test_binding_not_imported(self, mocker):
        """Check that the finder looks for a binding if none worked."""
        mocker.patch('sys.platform', 'some-supported-platform')
        mocker.patch(f'{MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.IMPORT_ERROR,
                     binding=None)
        mock_finder_cls = mocker.patch(f'{MODULE}.QtDependencyFinder')
        find_missing_qt_dependencies()
        mock_finder_cls.assert_called_once_with(None)

    def test_no_binding_installed(self, mocker):
        """Check the result when no binding is installed at all."""
        mocker.patch('sys.platform', 'some-supported-platform')
        mocker.patch(f'{MODULE}.check_pyqt_sanity',
                     return_value=PyQtSanity.IMPORT_ERROR,
                     binding=None)
        mocker.patch(f'{MODULE}.QtDependencyFinder',
                     side_effect=QtBindingNotFoundError)
        assert not find_missing_qt_dependencies()


class TestSuppressFilePermissionWarnings:
    """Tests for the suppress_file_permission_warnings function."""

    def test_no_xdg_runtime_dir(self, mocker):
        """Check no changes if XDG_RUNTIME_DIR is not set (e.g., not UNIX)."""
        mocker.patch(f'{MODULE}._UNIX_RUNTIME_DIR', '')
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

        mocker.patch(f'{MODULE}._UNIX_RUNTIME_DIR', str(tmp_path))
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

        mocker.patch(f'{MODULE}._UNIX_RUNTIME_DIR', str(fails))

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
