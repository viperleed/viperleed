"""Tests for module folder of viperleed.calc.bookkeeper.history."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-19'
__license__ = 'GPLv3+'

from contextlib import contextmanager
from unittest.mock import patch

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.errors import MetadataMismatchError
from viperleed.calc.bookkeeper.history.folder import HistoryFolder
from viperleed.calc.bookkeeper.history.folder import IncompleteHistoryFolder
from viperleed.calc.lib.dataclass_utils import set_frozen_attr


_MODULE = 'viperleed.calc.bookkeeper.history.folder'
_FAKE_HASH = 'this is a fake hash value'
_VALID_FOLDER_NAME = 't001.r002_other_stuff'


class MockMetaFile:
    """A fake BookkeeperMetaFile."""

    def __init__(self, path):
        """Initialize from a path."""
        self.path = path
        self.hash_ = _FAKE_HASH
        self.parent = None

    @property
    def file(self):
        """Return the path to this fake BookkeeperMetaFile."""
        return self.path / 'fake_file'

    def read(self):
        """Don't read anything."""

    def compute_hash(self):
        """Don't compute anything."""


@fixture(name='incomplete_folder')
def mock_incomplete_folder(mock_path):
    """Return a mocked IncompleteHistoryFolder."""
    mock_path.name = _VALID_FOLDER_NAME
    mock_path.is_dir.return_value = True
    return IncompleteHistoryFolder(mock_path)


@fixture(name='history_folder')
def mock_history_folder(incomplete_folder, metafile):
    """Return a mocked HistoryFolder."""
    with metafile() as mock_meta:
        folder = HistoryFolder(incomplete_folder.path)
        set_frozen_attr(folder, 'metadata', mock_meta(folder.path))
        return folder


@fixture(name='metafile')
def mock_metafile(monkeypatch):
    """Return a context with a fake BookkeeperMetaFile."""
    @contextmanager
    def _context(**kwargs):
        fake_meta = kwargs.get('fake_meta', MockMetaFile)
        with monkeypatch.context() as patch_:
            patch_.setattr(
                'viperleed.calc.bookkeeper.history.meta.BookkeeperMetaFile',
                fake_meta,
                )
            patch_.setattr(f'{_MODULE}.BookkeeperMetaFile', fake_meta)
            yield fake_meta
    return _context


class TestIncompleteHistoryFolder:
    """Tests for the IncompleteHistoryFolder dataclass."""

    _cls = IncompleteHistoryFolder

    def test_init(self, mock_path, incomplete_folder):
        """Test __post_init__ of IncompleteHistoryFolder."""
        # pylint: disable=magic-value-comparison
        assert incomplete_folder.tensor_num == 1
        assert incomplete_folder.job_num == 2
        assert incomplete_folder.path == mock_path

    def test_invalid_name(self, mock_path):
        """Test invalid folder name raises ValueError in _analyze_path."""
        mock_path.name = 'invalid_name'
        with pytest.raises(ValueError, match='Invalid'):
            self._cls(mock_path)

    def test_exists(self, incomplete_folder):
        """Test exists property of IncompleteHistoryFolder."""
        assert incomplete_folder.exists is True

    @parametrize(file=(True, False))
    def test_copy_file_or_directory(self, file, incomplete_folder, mock_path):
        """Test copy_file_or_directory method for a file."""
        mock_path.is_file.return_value = file
        mock_shutil = 'shutil.copy2' if file else 'shutil.copytree'

        with patch(mock_shutil) as mock_copy:
            incomplete_folder.copy_file_or_directory(mock_path)
            target = incomplete_folder.path / mock_path.name
            mock_copy.assert_called_once_with(mock_path, target)

    def test_copy_file_fails(self, incomplete_folder, mock_path):
        """Check complaints when failing to copy a file."""
        patch_shutil = patch('shutil.copy2', side_effect=OSError)
        patch_log = patch(f'{_MODULE}.LOGGER.error')
        with patch_shutil, patch_log as mock_error:
            incomplete_folder.copy_file_or_directory(mock_path, 'nowhere')
            mock_error.assert_called_once()


class TestHistoryFolder(TestIncompleteHistoryFolder):
    """Tests for the HistoryFolder class."""

    _cls = HistoryFolder

    # pylint: disable=arguments-renamed  # It's a fixture
    def test_init(self, mock_path, history_folder):
        """Test HistoryFolder __post_init__."""
        assert history_folder.hash_ == _FAKE_HASH
        assert not history_folder.parent
        super().test_init(mock_path, history_folder)

    @parametrize(is_file=(True, False))
    def test_has_metadata(self, mock_path, is_file, history_folder):
        """Test has_metadata property of HistoryFolder."""
        with patch.object(mock_path, 'is_file', return_value=is_file):
            assert history_folder.has_metadata == is_file

    def test_check_metadata_ok(self, history_folder, metafile):
        """Test check_metadata method without mismatch."""
        with metafile():
            history_folder.check_metadata()

    def test_check_metadata_raises(self, history_folder, metafile):
        """Test check_metadata raises MetadataMismatchError on mismatch."""
        class _OtherHashMeta(MockMetaFile):
            def __init__(self, *args):
                super().__init__(*args)
                self.hash_ = 'another different hash'

        with metafile(fake_meta=_OtherHashMeta):
            with pytest.raises(MetadataMismatchError):
                history_folder.check_metadata()

    def test_analyze_path_not_directory(self, mock_path):
        """Test _analyze_path raises ValueError if path is not a directory."""
        mock_path.is_dir.return_value = False
        with pytest.raises(ValueError):
            HistoryFolder(mock_path)

    def test_analyze_path_metadata_not_found(self, mock_path, metafile):
        """Test _analyze_path logs warning if metadata file is missing."""
        class _RaiseOnRead(MockMetaFile):
            def read(self):
                raise FileNotFoundError
        patch_ = patch(f'{_MODULE}.LOGGER.warning')
        with metafile(fake_meta=_RaiseOnRead), patch_ as mock_warning:
            mock_path.name = _VALID_FOLDER_NAME
            HistoryFolder(mock_path)
            mock_warning.assert_called_once()
