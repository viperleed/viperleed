"""Tests for module leedbase of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-31'
__license__ = 'GPLv3+'

from fnmatch import fnmatch  # To simulate globbing
from pathlib import Path

import pytest
from pytest_cases import fixture

from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.lib.leedbase import get_tensor_indices
from viperleed.calc.lib.leedbase import getMaxTensorIndex


_MODULE = 'viperleed.calc.lib.leedbase'


@fixture(name='mock_tensor')
def factory_mock_tensor(mocker):
    """Return a fake Tensor file with a given stem."""
    def _make(stem, **kwargs):
        suffix = kwargs.get('suffix', '')
        # Only one between "is_file" and "is_dir" can be True:
        is_file, is_dir = kwargs.get('is_file'), kwargs.get('is_dir')
        if is_file:
            kwargs['is_dir'] = False
        elif is_dir:
            kwargs['is_file'] = False
        else:  # Neither file nor directory (e.g., symlink)
            kwargs['is_file'] = kwargs['is_dir'] = False
        # Replace the relevant is_file/is_dir with a callable
        for key in ('is_file', 'is_dir'):
            kwargs[key] = mocker.MagicMock(return_value=kwargs[key])
        mocked = mocker.MagicMock(spec=Path, stem=stem, **kwargs)
        mocked.name = stem+suffix
        return mocked
    return _make


def mock_glob(*fake_files):
    """Return a callable that yields items in fake_files matching a pattern."""
    def _iter(pattern):
        yield from (f for f in fake_files if fnmatch(f.name, pattern))
    return _iter


class AutousePatchPath:
    """A base class for tests using patch_path as a fixture for all tests."""

    @fixture(name='patch_path', autouse=True)
    def fixture_patch_path(self, mock_path, mocker):
        """Patch leedbase.Path to return `mock_path` instead."""
        mock_path.resolve = mocker.MagicMock(return_value=mock_path)
        return mocker.patch(f'{_MODULE}.Path', return_value=mock_path)


class TestGetMaxTensorIndex(AutousePatchPath):
    """Tests for the getMaxTensorIndex function."""

    def test_no_files(self, mock_path):
        """Test getMaxTensorIndex without tensor files/directories."""
        mock_path.is_dir.return_value = True
        mock_path.glob.return_value = []
        index = getMaxTensorIndex()
        assert index == 0

    def test_files_and_folders(self, mock_path, mock_tensor):
        """Test getMaxTensorIndex with files/folders with different indices."""
        mock_path.is_dir.return_value = True
        mock_path.glob.side_effect = mock_glob(
            mock_tensor(f'{DEFAULT_TENSORS}_001', is_dir=True),
            mock_tensor(f'{DEFAULT_TENSORS}_002', suffix='.zip', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_003', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_004', suffix='.zip', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_005', is_dir=True),
            mock_tensor(f'{DEFAULT_TENSORS}_bad_file', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_bad_folder', is_dir=True),
            mock_tensor(f'{DEFAULT_TENSORS}_009'),  # Not a file/folder
            mock_tensor(f'Not_a_tensor_012', is_dir=True),
            )
        assert getMaxTensorIndex() == 5

    def test_zip_only(self, mock_path, mock_tensor):
        """Test getMaxTensorIndex with zip_only=True."""
        mock_path.is_dir.return_value = True
        mock_path.glob.side_effect = mock_glob(
            mock_tensor(f'{DEFAULT_TENSORS}_001', suffix='.zip', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_004', suffix='.zip', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_008', is_dir=True),   # A folder
            mock_tensor(f'{DEFAULT_TENSORS}_009', is_file=True),  # A non-zip file
            )
        result = getMaxTensorIndex(zip_only=True)
        assert result == 4

    def test_mixed_files(self, mock_path, mock_tensor):
        """Test getMaxTensorIndex with a mix of .zip, folders, and files."""
        mock_path.is_dir.return_value = True
        mock_path.glob.side_effect = mock_glob(
            mock_tensor(f'{DEFAULT_TENSORS}_002', is_dir=True),
            mock_tensor(f'{DEFAULT_TENSORS}_005', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_010', suffix='.zip', is_file=True),
            )
        index = getMaxTensorIndex()
        assert index == 10


class TestGetTensorIndices(AutousePatchPath):
    """Tests for the get_tensor_indices function."""

    def test_files_and_folders(self, mock_path, mock_tensor):
        """Test get_tensor_indices with tensor files/directories present."""
        mock_path.is_dir.return_value = True
        mock_path.glob.side_effect = mock_glob(
            mock_tensor(f'{DEFAULT_TENSORS}_001', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_002', suffix='.zip', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_003', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_004', suffix='.zip', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_005', suffix='.zip', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_005', is_dir=True),  # Unzipped
            mock_tensor(f'{DEFAULT_TENSORS}_123456',
                        suffix='.zip',
                        is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_bad', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_099_not_int', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_bad_folder', is_dir=True),
            mock_tensor(f'{DEFAULT_TENSORS}_123_not_int_folder', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_009'),  # Not a file/folder
            mock_tensor(f'Not_a_tensor_012', is_dir=True),
            )
        result = set(get_tensor_indices())
        assert result == {1, 2, 3, 4, 5, 123456}

    def test_invalid_files(self, mock_path, mock_tensor):
        """Test get_tensor_indices with files that do not match the pattern."""
        mock_path.is_dir.return_value = True
        mock_path.glob.side_effect = mock_glob(
            mock_tensor(f'{DEFAULT_TENSORS}_abc', is_dir=True),
            mock_tensor('SomeOtherFile_123', is_file=True),
            )
        with pytest.raises(StopIteration):
            next(get_tensor_indices())

    def test_no_files(self, mock_path):
        """Test get_tensor_indices when Tensors folder is empty."""
        mock_path.is_dir.return_value = True
        mock_path.glob.return_value = []
        with pytest.raises(StopIteration):
            next(get_tensor_indices())

    def test_no_tensor_folder(self, mock_path):
        """Test get_tensor_indices when the Tensors folder does not exist."""
        mock_path.is_dir.return_value = False
        with pytest.raises(StopIteration):
            next(get_tensor_indices())

    def test_zip_files_only(self, mock_path, mock_tensor):
        """Test get_tensor_indices with zip_only=True."""
        mock_path.is_dir.return_value = True
        mock_path.glob.side_effect = mock_glob(
            mock_tensor(f'{DEFAULT_TENSORS}_001', suffix='.zip', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_002', suffix='.zip', is_file=True),
            mock_tensor(f'{DEFAULT_TENSORS}_003', is_dir=True),  # A folder
            )
        result = set(get_tensor_indices(zip_only=True))
        assert result == {1, 2}
