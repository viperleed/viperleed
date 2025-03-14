"""Tests for module explorer of viperleed.calc.bookkeeper.history."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-19'
__license__ = 'GPLv3+'

import ast
from operator import attrgetter
from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.explorer import HistoryExplorer
from viperleed.calc.bookkeeper.history.errors import CantRemoveEntryError
from viperleed.calc.bookkeeper.history.errors import MetadataMismatchError
from viperleed.calc.bookkeeper.history.folder import HistoryFolder
from viperleed.calc.constants import DEFAULT_HISTORY

from ....helpers import not_raises
from ....helpers import raises_exception

_MODULE = 'viperleed.calc.bookkeeper.history.explorer'


@fixture(name='history')
def mock_explorer(mock_path):
    """Return a HistoryExplorer instance with a mocked root path."""
    return HistoryExplorer(mock_path)


@fixture(name='patched_folder')
def fixture_patched_folder(mocker):
    """Return a fake HistoryFolder class."""
    return mocker.patch(f'{_MODULE}.HistoryFolder')


class TestHistoryExplorer:
    """Tests for the HistoryExplorer class (except collection)."""

    def test_init(self, history, mock_path):
        """Test initialization of HistoryExplorer."""
        # pylint: disable=protected-access                # OK in tests
        assert history.path == (mock_path / DEFAULT_HISTORY)
        assert not any(history._subfolders)
        assert all(m is None for m in history._maps.values())
        assert history._new_calc_run_folder is None
        assert history._info is None

    @parametrize(backed_up=(Path('.bak'), None))
    def test_fix(self, backed_up, history, patched_folder, mocker):
        """Check correct behavior of the .fix() method."""
        if backed_up:
            mocker.patch('pathlib.Path.read_text', return_value='some text')
        mocker.patch.object(history, '_info')
        mock_log = mocker.patch(f'{_MODULE}.LOGGER.info')
        fake_folder = patched_folder()
        # pylint: disable=protected-access                # OK in tests
        history._subfolders = [fake_folder]
        called = (
            mocker.patch.object(history,
                                '_backup_info_file',
                                return_value=backed_up),
            mocker.patch.object(history.info, 'fix'),
            # Only works for a single folder, as calling patched_folder
            # repeatedly always returns the same MagicMock object.
            mocker.patch.object(fake_folder, 'fix'),
            )
        history.fix()
        for fixer in called:
            fixer.assert_called_once()
        assert mock_log.call_count == (2 if backed_up else 0)

    def test_discard_last(self, history, mocker):
        """Test the discard_most_recent_run method."""
        rmtree = mocker.patch('shutil.rmtree')
        fake_paths = tuple(range(5))
        mocker.patch.object(history,
                            'list_paths_to_discard',
                            return_value=fake_paths)
        history.discard_most_recent_run()
        rmtree.assert_has_calls(mocker.call(p) for p in fake_paths)

    def test_discard_last_raises(self, history, mocker):
        """Check that exceptions while removing folders are propagated."""
        mocker.patch.object(history,
                            'list_paths_to_discard',
                            return_value=(history.path/'some_subfolder',))
        with raises_exception('shutil.rmtree', OSError):
            history.discard_most_recent_run()

    def test_list_to_discard(self, history, mocker):
        """Test the list_paths_to_discard method."""
        folders = *_, last = [
            mocker.MagicMock(spec=HistoryFolder, path=mocker.MagicMock())
            for _ in range(3)
            ]
        for i, folder in enumerate(folders):
            folder.path.name = f'folder_{i}'
        # pylint: disable=protected-access                # OK in tests
        history._maps['hash_to_parent'] = {last.hash_: last}
        history._maps['main_hash_to_folders'] = {last.hash_: set(folders)}
        history._subfolders = folders
        # pylint: enable=protected-access
        expect = tuple(f.path for f in folders)
        assert history.list_paths_to_discard() == expect

    def test_prepare_info_file(self, history, mocker):
        """Test prepare_info_file method."""
        mock_info_file = mocker.patch(f'{_MODULE}.HistoryInfoFile')
        history.prepare_info_file()
        mock_info_file.assert_called_once_with(history.root,
                                               create_new=True)
        assert history.info == mock_info_file.return_value

    def test_register_folder(self, history, mock_path, mocker):
        """Test register_folder method."""
        history.collect_subfolders()
        mock_append = mocker.patch.object(history, '_append_existing_folder')
        mock_update = mocker.patch.object(history, '_update_maps')
        history.register_folder(mock_path)
        mock_append.assert_called_once_with(mock_path, insert_sorted=True)
        mock_update.assert_called_once()

    def test_append_existing_folder_invalid_path(self, history,
                                                 mock_path, mocker):
        """Test _append_existing_folder with invalid path raises ValueError."""
        mocker.patch.object(mock_path, 'parent')
        with pytest.raises(ValueError, match='Not a subfolder'):
            # pylint: disable-next=protected-access       # OK in tests
            history._append_existing_folder(mock_path)

    def test_append_existing_folder_valid_path(self, history, mock_path,
                                               patched_folder, mocker):
        """Test _append_existing_folder with valid path."""
        history.collect_subfolders()
        mock_path.parent = history.path
        patched_folder.return_value = folder = mocker.MagicMock()
        # pylint: disable-next=protected-access           # OK in tests
        history._append_existing_folder(mock_path)
        # pylint: disable-next=protected-access           # OK in tests
        assert folder in history._subfolders

    @parametrize(info_exists=(True, False))
    def test_backup_info(self, info_exists, history, mocker):
        """Check the backing up of the history.info file."""
        mocker.patch(f'{_MODULE}.HistoryInfoFile')
        history.prepare_info_file()
        info = history.info.path
        mocker.patch.object(info, 'is_file', return_value=info_exists)
        # Avoid infinite loop, and patch .replace
        mocker.patch('pathlib.Path.is_file',
                     new=lambda p: not p.name.endswith('4'))
        mock_replace = mocker.patch.object(info, 'replace')
        # pylint: disable-next=protected-access           # OK in tests
        backed_up = history._backup_info_file()
        if info_exists:
            assert backed_up == Path(f'{info}.bak4')
            mock_replace.assert_called_once_with(backed_up)
        else:
            assert backed_up is None
            mock_replace.assert_not_called()

    _too_early_attr = {
        'info': 'prepare_info_file',
        'new_folder': 'find_new_history_directory',
        'last_folder': None,
        'max_run_per_tensor': None,
        }
    _too_early_call = {
        'fix': 'prepare_info_file',
        'find_new_history_directory(None, None)': None,
        'register_folder': None,
        '_backup_info_file': 'prepare_info_file',
        '_find_name_for_new_history_subfolder(None, None)': None,
        }

    @parametrize('attr,updater', _too_early_attr.items(), ids=_too_early_attr)
    def test_too_early_attr(self, attr, updater, history):
        """Check complaints when accessing attributes too early."""
        # pylint: disable-next=protected-access           # OK in tests
        updater = updater or HistoryExplorer._updater
        with pytest.raises(AttributeError, match=rf'.*{updater}.*'):
            attrgetter(attr)(history)

    @parametrize('method_name,updater',
                 _too_early_call.items(),
                 ids=_too_early_call)
    def test_too_early_method_call(self, method_name, updater, history):
        """Check that accessing attributes before update_from_cwd fails."""
        # pylint: disable-next=magic-value-comparison
        if '(' not in method_name:
            args = tuple()
        else:
            method_name, args_str = method_name.split('(')
            args_str = args_str.replace(')', '') + ','
            args = ast.literal_eval(args_str)
        method = attrgetter(method_name)(history)
        # pylint: disable-next=protected-access           # OK in tests
        updater = updater or HistoryExplorer._updater
        with pytest.raises(AttributeError, match=rf'.*{updater}.*'):
            method(*args)


class TestHistoryExplorerCollection:
    """Tests for collection of history subfolders and related methods."""

    @staticmethod
    def _check_collected_nothing(history):
        """Check that no folders were collected."""
        # pylint: disable-next=protected-access           # OK in tests
        assert not any(history._subfolders)
        # pylint: disable-next=protected-access           # OK in tests
        maps = history._maps.values()
        assert all(not m for m in maps)
        assert all(m is not None for m in maps)

    def test_collect_subfolders_empty(self, history, patched_folder):
        """Test collect_subfolders when there are no subfolders."""
        history.path.iterdir.return_value = []
        history.collect_subfolders()
        patched_folder.assert_not_called()
        self._check_collected_nothing(history)

    def test_collect_subfolders_no_history(self, history,
                                           patched_folder,
                                           mocker):
        """Test collect_subfolders when there are no subfolders."""
        mocker.patch.object(history.path, 'is_dir', return_value=False)
        mock_warn = mocker.patch(f'{_MODULE}.LOGGER.warning')
        history.collect_subfolders()
        patched_folder.assert_not_called()
        self._check_collected_nothing(history)
        mock_warn.asset_called_once()

    def test_collect_subfolders_with_folders(self, history,                     # TODO: this needs some realistic cases
                                             patched_folder,
                                             mocker):
        """Test collect_subfolders with valid folders."""
        directories = [mocker.MagicMock(spec=Path) for _ in range(2)]
        for i, directory in enumerate(directories):
            directory.name = f't00{i}.r001_rest'
            directory.parent = history.path
        mocker.patch.object(history.path, 'iterdir', return_value=directories)
        patched_folder.return_value.parent = None
        history.collect_subfolders()
        # pylint: disable-next=protected-access           # OK in tests
        assert len(history._subfolders) == len(directories)

    def test_collect_subfolders_funny(self, history, mocker):
        """Test collect_subfolders with valid folders."""
        directories = [mocker.MagicMock(spec=Path) for _ in range(2)]
        for i, directory in enumerate(directories):
            directory.name = f't00{i}.r001_rest'
        mocker.patch.object(history.path, 'iterdir', return_value=directories)
        mocker.patch.object(history,
                            '_append_existing_folder',
                            side_effect=ValueError)
        history.collect_subfolders()
        self._check_collected_nothing(history)

    def test_find_name_for_new_history_subfolder(self, history):
        """Test _find_name_for_new_history_subfolder."""
        # pylint: disable-next=protected-access           # OK in tests
        history._maps['jobs_for_tensor'] = {1: {1}}
        args = 1, 'test'
        # pylint: disable-next=protected-access           # OK in tests
        result = history._find_name_for_new_history_subfolder(*args)
        assert result == f't00{args[0]}.r002_{args[1]}'

    def test_find_new_history_directory(self, history, mocker):
        """Test find_new_history_directory creates a new folder."""
        history.collect_subfolders()
        mocker.patch.object(history,
                            '_find_name_for_new_history_subfolder',
                            return_value='new_folder')
        mock_folder = mocker.patch(f'{_MODULE}.IncompleteHistoryFolder')
        history.find_new_history_directory(1, 'test')
        mock_folder.assert_called_once()
        assert history.new_folder

    _has_subfolders = {
        r't001\.r002.*_timestamp': True,
        r't001\.r002.*_other_timestamp': False,
        }

    @parametrize('name,expect', _has_subfolders.items())
    def test_has_subfolder(self, name, expect, history, mocker):
        """Check the expected result of the has_subfolder method."""
        subfolder = mocker.MagicMock()
        # We cant pass 'name' as a kwarg to MagicMock, as it
        # is a special initialization argument for MagicMock.
        subfolder.name = 't001.r002.more-text_timestamp'
        # pylint: disable-next=protected-access           # OK in tests
        history._subfolders = [subfolder]
        result = history.has_subfolder(name)
        assert result == expect

    def test_last_folder_no_subfolders(self, history):
        """Test last_folder when there are no subfolders."""
        history.collect_subfolders()
        assert history.last_folder is None
        assert not any(history.last_folder_and_siblings)

    def test_last_folder_with_subfolders(self, history, mocker):
        """Test last_folder returns the most recent folder."""
        mock_folder_last = mocker.MagicMock()
        sibling_folders = [mocker.MagicMock() for _ in range(3)]
        mock_folders = (
            mocker.MagicMock(),  # Some other folder
            *sibling_folders,
            mock_folder_last,
            )
        # pylint: disable=protected-access                # OK in tests
        history._subfolders = mock_folders
        history._maps['hash_to_parent'] = {mock_folder_last.hash_:
                                           mock_folder_last}
        history._maps['main_hash_to_folders'] = {
            mock_folder_last.hash_: {mock_folder_last, *sibling_folders},
            }
        # pylint: enable=protected-access
        assert history.last_folder is mock_folder_last
        # last_folder_and_siblings is unsorted!
        last_and_siblings = set(history.last_folder_and_siblings)
        assert last_and_siblings == {mock_folder_last, *sibling_folders}

    def test_max_run_per_tensor_empty(self, history):
        """Test max_run_per_tensor when no jobs exist."""
        history.collect_subfolders()
        max_runs = history.max_run_per_tensor
        assert not max_runs
        assert not max_runs['non existing key']  # == 0

    def test_max_run_per_tensor_with_jobs(self, history):
        """Test max_run_per_tensor when jobs exist for tensors."""
        # pylint: disable-next=protected-access           # OK in tests
        history._maps['jobs_for_tensor'] = {1: {1, 2}, 2: {3}}
        expected_result = {1: 2, 2: 3}
        assert history.max_run_per_tensor == expected_result

    def test_update_maps(self, history, mocker):
        """Test _update_maps method."""
        folders = mocker.MagicMock(), mocker.MagicMock()
        main = folders[0]
        main.parent = None
        # pylint: disable-next=protected-access           # OK in tests
        history._maps['main_hash_to_folders'] = {None: folders}
        # pylint: disable-next=protected-access           # OK in tests
        history._update_maps()
        expect = {f.hash_: main for f in folders}
        # pylint: disable-next=protected-access           # OK in tests
        assert history._maps['hash_to_parent'] == expect


class TestHistoryExplorerConsistencyCheck:
    """Tests for consistency checks between history folder and info entry."""

    @fixture(name='add_subfolder')
    def fixture_add_fake_subfolder(self, mocker):
        """Add a fake subfolder to a HistoryExplorer."""
        def _add(explorer):
            last_folder = mocker.MagicMock(spec=HistoryFolder)
            # pylint: disable-next=protected-access       # OK in tests
            explorer._subfolders = [last_folder]
            # pylint: disable-next=protected-access       # OK in tests
            explorer._maps = {
                'hash_to_parent': {last_folder.hash_: last_folder}
                }
            return last_folder
        return _add

    def test_missing_last_folder(self, history, add_subfolder):
        """Check complaints when last_folder does not exist."""
        last_folder = add_subfolder(history)
        last_folder.exists = False
        with pytest.raises(FileNotFoundError):
            history.check_last_folder_consistent()

    def test_no_last_folder(self, history):
        """Check complaints when no last_folder is found."""
        history.collect_subfolders()
        with pytest.raises(FileNotFoundError):
            history.check_last_folder_consistent()

    def test_metadata_mismatch(self, history, add_subfolder):
        """Check complaints when metadata are outdated."""
        last_folder = add_subfolder(history)
        last_folder.exists = True
        last_folder.check_metadata.side_effect = MetadataMismatchError
        with pytest.raises(MetadataMismatchError):
            history.check_last_folder_consistent()

    def test_inconsistent_entry(self, history, add_subfolder, mocker):
        """Check complaints if the last folder and entry are inconsistent."""
        mocker.patch.object(history, '_info')
        last_folder = add_subfolder(history)
        last_folder.exists = True
        last_folder.check_consistent_with_entry.side_effect = (
            CantRemoveEntryError
            )
        with pytest.raises(CantRemoveEntryError):
            history.check_last_folder_consistent()

    def test_consistent(self, history, add_subfolder, mocker):
        """Test successful call to check_last_folder_consistent."""
        mocker.patch.object(history, '_info')
        last_folder = add_subfolder(history)
        last_folder.exists = True
        with not_raises(Exception):
            history.check_last_folder_consistent()
        last_folder.check_metadata.assert_called_once()
        last_folder.check_consistent_with_entry.assert_called_once_with(
            history.info.last_entry
            )
