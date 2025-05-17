"""Tests for module folder of viperleed.calc.bookkeeper.history."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-19'
__license__ = 'GPLv3+'

from operator import attrgetter
from pathlib import Path

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
_VALID_TIMESTAMP = '010203-040506'
_VALID_FOLDER_NAME = f't001.r002_other_stuff_{_VALID_TIMESTAMP}'


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

    _valid = {
        'not moved': _VALID_FOLDER_NAME,
        'moved': f't001.r002_stuff_moved-{_VALID_TIMESTAMP}'
        }

    @parametrize(name=_valid.values(), ids=_valid)
    def test_init(self, name, mock_path):
        """Test __post_init__ of IncompleteHistoryFolder."""
        mock_path.name = name
        folder = self._cls(mock_path)
        # pylint: disable=magic-value-comparison
        assert folder.tensor_num == 1
        assert folder.job_num == 2
        assert folder.path == mock_path
        # pylint: disable-next=no-member                # Inference
        assert folder.timestamp.endswith(_VALID_TIMESTAMP)

    _invalid = {
        'no tensors/job': 'invalid_name',
        'no timestamp': 't001.r002_stuff',
        }

    @parametrize(name=_invalid.values(), ids=_invalid)
    def test_invalid_name(self, name, mock_path):
        """Test invalid folder name raises ValueError in _analyze_path."""
        mock_path.name = name
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

    _valid = TestIncompleteHistoryFolder._valid

    @parametrize(name=_valid.values(), ids=_valid)
    @pytest.mark.usefixtures('logs')
    # pylint: disable-next=arguments-differ
    def test_init(self, name, mock_path, patch_metafile):
        """Test HistoryFolder __post_init__."""
        patch_metafile()
        mock_path.name = name
        folder = self._cls(mock_path)
        assert folder.hash_ == _FAKE_HASH
        assert not folder.parent
        # pylint: disable-next=no-member                # Inference
        folder.logs.collect.assert_called_once()
        super().test_init(name, mock_path)

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


class TestHistoryFolderDomains:
    """Tests for marking history folders in a DOMAINS calculation."""

    @fixture(name='make_fake_path')
    def factory_make_fake_path(self, mocker):
        """Return a fake Path to a history subfolder."""
        folder_name = 't000.r001_010203-040506'
        def _make(str_value=None):
            fake_path = mocker.MagicMock(spec=Path)
            fake_path.name = folder_name
            fake_path.is_dir.return_value = True
            if str_value is not None:
                fake_path.__str__.return_value = str_value
            return fake_path
        return _make

    def test_mark_domain(self, make_fake_path, mocker):
        """Check successful marking of a history folder as a subdomain."""
        mocker.patch(f'{_MODULE}.BookkeeperMetaFile.compute_hash')
        main_folder = HistoryFolder(make_fake_path())
        # pylint: disable-next=protected-access,assigning-non-slot
        main_folder.metadata._hash = 'main_hash'
        main_path = make_fake_path('main_path')
        domain = HistoryFolder(make_fake_path())
        domain.mark_as_domain(main_path, main_folder)
        # pylint: disable-next=no-member  # Inference
        assert domain.metadata.domains == {'main': ('main_path', 'main_hash')}

    _domains = {
        'no subdomains': ((), {}),
        'with subdomains': (
            (('path_1', 'hash_1'), ('path_2', 'hash_2')),
            {'domains': (('path_1', 'hash_1'), ('path_2', 'hash_2'))},
            ),
        'with subdomains none': (
            (('path_1', 'hash_1'), ('path_2', 'hash_2'), ('skipped', None)),
            {'domains': (('path_1', 'hash_1'), ('path_2', 'hash_2'))},
            ),
        }

    @parametrize('domains,expect', _domains.values(), ids=_domains)
    def test_mark_main(self, domains, expect, make_fake_path, mocker):
        """Check successful marking of a history folder as the main one."""
        mocker.patch(f'{_MODULE}.BookkeeperMetaFile.compute_hash')
        domain_paths = []
        domain_folders = []
        for domain_path, domain_hash in domains:
            if domain_hash:
                folder = HistoryFolder(make_fake_path())
                # pylint: disable-next=protected-access,assigning-non-slot
                folder.metadata._hash = domain_hash
            else:
                folder = None
            domain_paths.append(make_fake_path(domain_path))
            domain_folders.append(folder)
        main_path = make_fake_path()
        history_folder = HistoryFolder(main_path)
        history_folder.mark_as_domains_main(domain_paths, domain_folders)
        # pylint: disable-next=no-member  # Inference
        assert history_folder.metadata.domains == expect

    def test_mark_main_raises(self, history_folder):
        """Check complaints when inconsistent paths and folders are given."""
        with pytest.raises(ValueError, match='Inconsistent number'):
            history_folder.mark_as_domains_main((1, 2, 3), (1, 2))


class TestHistoryFolderConsistency:
    """Collection of tests for consistency checks of a HistoryFolder."""

    def test_check_metadata_ok(self, history_folder):
        """Test check_metadata method without mismatch."""
        with not_raises(MetadataMismatchError):
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
