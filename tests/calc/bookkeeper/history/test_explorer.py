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
from viperleed.calc.bookkeeper.history.meta import _METADATA_NAME
from viperleed.calc.constants import DEFAULT_HISTORY

from ....helpers import filesystem_from_dict
from ....helpers import filesystem_to_dict
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

    @parametrize(fix_file=(True, False))
    @parametrize(fix_folders=(True, False))
    def test_fix(self, fix_file, fix_folders, history, mocker):
        """Check the expected behavior of the .fix() method."""
        mocker.patch.object(history, '_info')
        called = (
            mocker.patch.object(history,
                                '_fix_info_file',
                                return_value=fix_file),
            mocker.patch.object(history,
                                '_fix_subfolders',
                                return_value=fix_folders),
            )
        if fix_file or fix_folders:
            called += (mocker.patch(f'{_MODULE}.LOGGER.info'),)
        result = history.fix()
        assert result == fix_file or fix_folders
        for fixer in called:
            fixer.assert_called_once()

    _backup = {  # backed_up, same_text
        'no backup': (None, False),
        'backup different': (Path('.bak'), False),
        'backup same': (Path('.bak'), True),
        }

    @parametrize('backed_up,same_text', _backup.values(), ids=_backup)
    def test_fix_info(self, backed_up, same_text, history, mocker):
        """Check correct behavior of the _fix_info_file method."""
        back_up_text = 'some text'
        if backed_up:
            mocker.patch('pathlib.Path.read_text', return_value=back_up_text)
        info = mocker.patch.object(history, '_info')
        if same_text:
            info.raw_contents = back_up_text
        mock_log = mocker.patch(f'{_MODULE}.LOGGER.info')
        called = (
            mocker.patch.object(history,
                                '_backup_info_file',
                                return_value=backed_up),
            mocker.patch.object(history.info, 'fix'),
            )
        if same_text:
            called += mocker.patch('pathlib.Path.unlink')
        # pylint: disable-next=protected-access           # OK in tests
        result = history._fix_info_file()
        assert result == bool(backed_up and not same_text)
        for fixer in called:
            fixer.assert_called_once()
        assert mock_log.call_count == (1 if result else 0)

    @fixture(name='misses_medatata')
    def _fixture_misses_metadata(self, tmp_path):
        """Create a history tree with folders without metadata."""
        with_metadata = {
            't000.r001_000000-000001': {_METADATA_NAME: (
                '[archived]\nhash=7e6c9f1133c4ceae6e51bee48ae49c7f'
                )},
            # The next one should not have its metadata edited,
            # even if it does not have the correct parent
            't000.r001.001_DS_000000-000001': {_METADATA_NAME: (
                '[archived]\nhash=a105e6c5df761776065a38d3ac73668b'
                )},
            }
        no_metadata_only_main = {
            't000.r002_000000-000002': {},
            }
        no_metadata_with_siblings = {
            't001.r001_000000-000003': {},
            't000.r003.001_DS_000000-000003': {},
            't000.r003.002_DS_000000-000003': {},
            't001.r001.001_RDS_000000-000003': {},
            }
        no_metadata_duplicate_times = {
            't002.r001_000000-000004': {},
            't003.r001_000000-000004': {},
            't001.r002.001_DS_000000-000004': {},
            't002.r002.099_DS_000000-000004': {},
            }
        no_metadata_moved = {
            't004.r001_moved-000000-000005': {},
            't003.r002.r001_DS_moved-000000-000005': {},
            }
        tree = {'history': {
            **with_metadata,
            **no_metadata_only_main,
            **no_metadata_with_siblings,
            **no_metadata_duplicate_times,
            **no_metadata_moved,
            }}
        filesystem_from_dict(tree, tmp_path)

        expect_metadata = {
            # with_metadata: unchanged
            't000.r001_000000-000001': {
                'hash_': '7e6c9f1133c4ceae6e51bee48ae49c7f',
                'parent': None,
                },
            't000.r001.001_DS_000000-000001': {
                'hash_': 'a105e6c5df761776065a38d3ac73668b',
                'parent': None,
                },
            # no_metadata_only_main: added metadata, no parent
            't000.r002_000000-000002': {
                'hash_': 'f2606fc9ad8c6d6ce810069eaadb8f27',
                'parent': None,
                },
            # no_metadata_with_siblings: added metadata, marked parent
            't001.r001_000000-000003': {
                'hash_': '0b698119bf547e7dfcd68944fa6372b9',
                'parent': None,
                },
            't000.r003.001_DS_000000-000003': {
                'hash_': '88ec4994cc47ab9d8cb1b588728735e5',
                'parent': '0b698119bf547e7dfcd68944fa6372b9',
                },
            't000.r003.002_DS_000000-000003': {
                'hash_': '73d1cb6866abeb78b0c6336a605a6a58',
                'parent': '0b698119bf547e7dfcd68944fa6372b9',
                },
            't001.r001.001_RDS_000000-000003': {
                'hash_': '4726d09c1ffcd1229f028596871bf746',
                'parent': '0b698119bf547e7dfcd68944fa6372b9',
                },
            # no_metadata_duplicate_times: added metadata, no parent
            't002.r001_000000-000004': {
                'hash_': '2e6bde0dd1fd08517716fd6e5025c875',
                'parent': None,
                },
            't003.r001_000000-000004': {
                'hash_': '340101cb19b66ab4189f9fb525b21c9e',
                'parent': None,
                },
            't001.r002.001_DS_000000-000004': {
                'hash_': '64aedaa473d22240bad2d37c2cee58e3',
                'parent': None,
                },
            't002.r002.099_DS_000000-000004': {
                'hash_': 'f9ce3bdf28a434d0b0bb88a109536f26',
                'parent': None,
                },
            # no_metadata_moved: added metadata, no parent
            't004.r001_moved-000000-000005': {
                'hash_': 'ed112f0513d88962773726d7d1458252',
                'parent': None,
                },
            't003.r002.r001_DS_moved-000000-000005': {
                'hash_': '292f910c4812d65221db922010a7b9c7',
                'parent': None,
                },
            }
        # Ensure we haven't forgotten any
        assert all(f in expect_metadata for f in tree['history'])
        return expect_metadata

    def test_fix_subfolders(self, misses_medatata, tmp_path):
        """Ensure the correct behavior of the _fix_subfolders method."""
        history = HistoryExplorer(tmp_path)
        history.collect_subfolders()
        # pylint: disable-next=protected-access           # OK in tests
        assert history._fix_subfolders()

        # Check that all have a metadata file
        fixed_tree = filesystem_to_dict(tmp_path/'history')
        assert all(_METADATA_NAME in folder for folder in fixed_tree.values())

        # Now check the metadata contents
        for folder_name, meta_attrs in misses_medatata.items():
            # pylint: disable-next=protected-access       # OK in tests
            folder = next(f for f in history._subfolders
                          if f.name == folder_name)
            for attr_name, attr_value in meta_attrs.items():
                assert getattr(folder, attr_name) == attr_value
            # Make sure the parent was written to file,
            # and check that maps were updated correctly
            parent_hash = meta_attrs['parent']
            folder_hash = meta_attrs['hash_']
            # pylint: disable-next=protected-access       # OK in tests
            maps = history._maps
            # pylint: disable=E1135,E1136   # Pylint inference for maps
            if parent_hash:
                meta_contents = fixed_tree[folder_name][_METADATA_NAME]
                assert f'with = {parent_hash}' in meta_contents
                assert folder in maps['main_hash_to_folders'][parent_hash]
                assert maps['hash_to_parent'][folder_hash] is not folder
            else:
                assert parent_hash not in maps['main_hash_to_folders']
                assert maps['hash_to_parent'][folder_hash] is folder

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
        """Check complaints when `methods_name` is called too early."""
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

    def test_collect_subfolders_append_fails(self, history, mocker):
        """Test collect_subfolders when appending folders fails."""
        directories = [mocker.MagicMock(spec=Path) for _ in range(2)]
        for i, directory in enumerate(directories):
            directory.name = f't00{i}.r001_rest'
        mocker.patch.object(history.path, 'iterdir', return_value=directories)
        mocker.patch.object(history,
                            '_append_existing_folder',
                            side_effect=ValueError)
        history.collect_subfolders()
        self._check_collected_nothing(history)

    def test_collect_subfolders_empty(self, history, patched_folder):
        """Test collect_subfolders when there are no subfolders."""
        history.path.iterdir.return_value = []
        history.collect_subfolders()
        patched_folder.assert_not_called()
        self._check_collected_nothing(history)

    def test_collect_subfolders_no_history(self, history,
                                           patched_folder,
                                           mocker):
        """Test collect_subfolders when history does not exist."""
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
