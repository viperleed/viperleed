"""Tests for module folder of viperleed.calc.bookkeeper.history."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-19'
__license__ = 'GPLv3+'

from operator import attrgetter

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.errors import CantRemoveEntryError
from viperleed.calc.bookkeeper.history.errors import MetadataMismatchError
from viperleed.calc.bookkeeper.history.folder import HistoryFolder
from viperleed.calc.bookkeeper.history.folder import IncompleteHistoryFolder

from ....helpers import not_raises


_MODULE = 'viperleed.calc.bookkeeper.history.folder'
_FAKE_HASH = 'this is a fake hash value'
_VALID_FOLDER_NAME = 't001.r002_other_stuff'


class MockMetaFile:
    """A fake BookkeeperMetaFile."""

    def __init__(self, path):
        """Initialize from a path."""
        self.folder = path
        self.file = path / 'fake_file'
        self.hash_ = _FAKE_HASH
        self.parent = None

    def read(self):
        """Don't read anything."""

    def compute_hash(self):
        """Don't compute anything."""

    def write(self):
        """Don't write anything."""


@fixture(name='mock_entry')
def mock_history_info_entry(mocker):
    """Return a fake HistoryInfoEntry."""
    entry = mocker.MagicMock()
    entry.folder_name.value = _VALID_FOLDER_NAME
    entry.tensor_nums.no_tensors = False
    entry.tensor_nums.value = (1, 2, 5)
    entry.job_nums.value = (2, 3, 10)
    return entry


@fixture(name='incomplete_folder')
def mock_incomplete_folder(mock_path):
    """Return a mocked IncompleteHistoryFolder."""
    mock_path.name = _VALID_FOLDER_NAME
    mock_path.is_dir.return_value = True
    return IncompleteHistoryFolder(mock_path)


@fixture(name='logs')
def mock_logs(mocker):
    """Replace LogFiles with a mock."""
    return mocker.patch(f'{_MODULE}.LogFiles')


@fixture(name='history_folder')
# pylint: disable-next=unused-argument  # logs. Can't mark.usefixtures
def mock_history_folder(incomplete_folder, patch_metafile, logs):
    """Return a mocked HistoryFolder."""
    patch_metafile()
    folder = HistoryFolder(incomplete_folder.path)
    return folder


@fixture(name='patch_metafile')
def mock_metafile(mocker):
    """Replace BookkeeperMetaFile with another class."""
    def _patch(**kwargs):
        fake_meta = kwargs.get('fake_meta', MockMetaFile)
        mocker.patch(f'{_MODULE}.BookkeeperMetaFile', fake_meta)
    return _patch


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
    def test_copy_file_or_directory(self, file, incomplete_folder,
                                    mock_path, mocker):
        """Test copy_file_or_directory method for a file."""
        mock_path.is_file.return_value = file
        mock_copy = mocker.patch('shutil.copy2' if file else 'shutil.copytree')
        incomplete_folder.copy_file_or_directory(mock_path)
        target = incomplete_folder.path / mock_path.name
        mock_copy.assert_called_once_with(mock_path, target)

    @parametrize(with_name=(None, 'nowhere'))
    def test_copy_file_fails(self, with_name,
                             incomplete_folder,
                             mock_path,
                             mocker):
        """Check complaints when failing to copy a file."""
        mocker.patch('shutil.copy2', side_effect=OSError)
        mock_error = mocker.patch(f'{_MODULE}.LOGGER.error')
        incomplete_folder.copy_file_or_directory(mock_path, with_name)
        mock_error.assert_called_once()
        error_msg, *_ = mock_error.call_args[0]
        assert mock_path.name in error_msg
        if with_name:
            assert with_name in error_msg


class TestHistoryFolder(TestIncompleteHistoryFolder):
    """Tests for the HistoryFolder class."""

    _cls = HistoryFolder

    # pylint: disable=arguments-renamed  # It's a fixture
    def test_init(self, mock_path, history_folder):
        """Test HistoryFolder __post_init__."""
        assert history_folder.hash_ == _FAKE_HASH
        assert not history_folder.parent
        history_folder.logs.collect.assert_called_once()
        super().test_init(mock_path, history_folder)

    @parametrize(is_file=(True, False))
    def test_has_metadata(self, is_file, history_folder, mocker):
        """Test has_metadata property of HistoryFolder."""
        mocker.patch.object(history_folder.metadata.file,
                            'is_file',
                            return_value=is_file)
        assert history_folder.has_metadata == is_file

    def test_analyze_path_not_directory(self, mock_path):
        """Test _analyze_path raises ValueError if path is not a directory."""
        mock_path.is_dir.return_value = False
        with pytest.raises(ValueError):
            HistoryFolder(mock_path)

    def test_analyze_path_metadata_not_found(self, mock_path,
                                             patch_metafile, mocker):
        """Test _analyze_path logs warning if metadata file is missing."""
        class _RaiseOnRead(MockMetaFile):
            def read(self):
                raise FileNotFoundError
        mock_warning = mocker.patch(f'{_MODULE}.LOGGER.warning')
        patch_metafile(fake_meta=_RaiseOnRead)
        mock_path.name = _VALID_FOLDER_NAME
        HistoryFolder(mock_path)
        mock_warning.assert_called_once()

    @parametrize(meta_missing=(True, False))
    def test_fix(self, meta_missing, history_folder, mocker):
        """Check that .fix()ing a folder writes a metafile."""
        # Fake the absence of the metadata file
        meta = history_folder.metadata
        mocker.patch.object(meta.file,
                            'is_file',
                            return_value=not meta_missing)
        mocker.patch.object(meta, 'write')
        history_folder.fix()

        assert_called = (meta.write.assert_called_once if meta_missing
                         else meta.write.assert_not_called)
        assert_called()


class TestHistoryFolderConsistency:
    """Collection of tests for consistency checks of a HistoryFolder."""

    def test_check_metadata_ok(self, history_folder):
        """Test check_metadata method without mismatch."""
        history_folder.check_metadata()

    def test_check_metadata_raises(self, history_folder, patch_metafile):
        """Test check_metadata raises MetadataMismatchError on mismatch."""
        class _OtherHashMeta(MockMetaFile):
            def __init__(self, *args):
                super().__init__(*args)
                self.hash_ = 'another different hash'

        patch_metafile(fake_meta=_OtherHashMeta)
        with pytest.raises(MetadataMismatchError):
            history_folder.check_metadata()

    def test_check_entry_valid(self, history_folder, mock_entry):
        """Test that no error is raised when folder matches entry."""
        with not_raises(CantRemoveEntryError):
            history_folder.check_consistent_with_entry(mock_entry)

    _mismatched = {
        'folder name': {'folder_name.value': 'another name'},
        'tensor nums': {'tensor_nums.value': (3,)},
        'tensor nums, init': {'tensor_nums.no_tensors': True},
        'job nums': {'job_nums.value': (1,)},
        }

    @parametrize('attrs', _mismatched.values(), ids=_mismatched)
    def test_check_entry_mismatch(self, attrs, history_folder, mock_entry):
        """Test error raised when folder name does not match entry name."""
        for attr_name, attr_value in attrs.items():
            # pylint: disable-next=magic-value-comparison
            if '.' in attr_name:
                nested_attr_name, attr_name = attr_name.rsplit('.', 1)
                to_patch = attrgetter(nested_attr_name)(mock_entry)
            else:
                to_patch = mock_entry
            setattr(to_patch, attr_name, attr_value)
        with pytest.raises(CantRemoveEntryError):
            history_folder.check_consistent_with_entry(mock_entry)
