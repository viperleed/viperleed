"""Tests for module explorer of viperleed.calc.bookkeeper.history."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-19'
__license__ = 'GPLv3+'

import ast
from operator import attrgetter
from pathlib import Path
from unittest.mock import MagicMock
from unittest.mock import patch

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.explorer import HistoryExplorer
from viperleed.calc import DEFAULT_HISTORY


_MODULE = 'viperleed.calc.bookkeeper.history.explorer'
patch_folder = patch(f'{_MODULE}.HistoryFolder')


@fixture(name='history')
def mock_explorer(mock_path):
    """Return a HistoryExplorer instance with a mocked root path."""
    return HistoryExplorer(mock_path)


class TestHistoryExplorer:
    """Tests for the HistoryExplorer class."""

    def test_init(self, history, mock_path):
        """Test initialization of HistoryExplorer."""
        # pylint: disable=protected-access                # OK in tests
        assert history.path == (mock_path / DEFAULT_HISTORY)
        assert not history._subfolders
        assert all(m is None for m in history._maps.values())
        assert history._new_calc_run_folder is None
        assert history._info is None

    @staticmethod
    def _check_collected_nothing(history):
        """Check that no folders were collected."""
        # pylint: disable-next=protected-access           # OK in tests
        assert not history._subfolders
        # pylint: disable-next=protected-access           # OK in tests
        maps = history._maps.values()
        assert all(not m for m in maps)
        assert all(m is not None for m in maps)

    def test_collect_subfolders_empty(self, history):
        """Test collect_subfolders when there are no subfolders."""
        history.path.iterdir.return_value = []
        with patch_folder:
            history.collect_subfolders()
            self._check_collected_nothing(history)

    def test_collect_subfolders_with_folders(self, history):                    # TODO: this needs some realistic cases
        """Test collect_subfolders with valid folders."""
        directories = MagicMock(spec=Path), MagicMock(spec=Path)
        for i, directory in enumerate(directories):
            directory.name = f't00{i}.r001_rest'
            directory.parent = history.path
        history.path.iterdir.return_value = directories
        with patch_folder as mock_folder:
            mock_folder.return_value.parent = None
            history.collect_subfolders()
        # pylint: disable-next=protected-access           # OK in tests
        assert len(history._subfolders) == len(directories)

    def test_collect_subfolders_funny(self, history):
        """Test collect_subfolders with valid folders."""
        directories = MagicMock(spec=Path), MagicMock(spec=Path)
        history.path.iterdir.return_value = directories
        for i, directory in enumerate(directories):
            directory.name = f't00{i}.r001_rest'
        raises_valueerror = patch.object(history,
                                         '_append_existing_folder',
                                         side_effect=ValueError)
        with raises_valueerror:
            history.collect_subfolders()
        self._check_collected_nothing(history)

    def test_find_new_history_directory(self, history):
        """Test find_new_history_directory creates a new folder."""
        history.collect_subfolders()
        fake_new_folder_name = 'new_folder'
        patch_find_name = patch.object(history,
                                       '_find_name_for_new_history_subfolder',
                                       return_value=fake_new_folder_name)
        with patch(f'{_MODULE}.IncompleteHistoryFolder') as mock_folder:
            with patch_find_name:
                history.find_new_history_directory(1, 'test')
                mock_folder.assert_called_once()
                assert history.new_folder

    def test_prepare_info_file(self, history):
        """Test prepare_info_file method."""
        with patch(f'{_MODULE}.HistoryInfoFile') as mock_info_file:
            history.prepare_info_file()
            mock_info_file.assert_called_once_with(history.root,
                                                   create_new=True)
            assert history.info == mock_info_file.return_value

    def test_last_folder_no_subfolders(self, history):
        """Test last_folder when there are no subfolders."""
        history.collect_subfolders()
        assert history.last_folder is None

    def test_last_folder_with_subfolders(self, history):
        """Test last_folder returns the most recent folder."""
        mock_folders = MagicMock(), MagicMock()
        mock_folder_last = mock_folders[-1]
        # pylint: disable-next=protected-access           # OK in tests
        history._subfolders = mock_folders
        # pylint: disable-next=protected-access           # OK in tests
        history._maps['hash_to_parent'] = {mock_folder_last.hash_:
                                           mock_folder_last}
        assert history.last_folder is mock_folder_last

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

    def test_register_folder(self, history, mock_path):
        """Test register_folder method."""
        history.collect_subfolders()
        with patch.object(history, '_append_existing_folder') as mock_append:
            with patch.object(history, '_update_maps') as mock_update:
                history.register_folder(mock_path)
                mock_append.assert_called_once_with(mock_path,
                                                    insert_sorted=True)
                mock_update.assert_called_once()

    def test_append_existing_folder_invalid_path(self, history, mock_path):
        """Test _append_existing_folder with invalid path raises ValueError."""
        mock_path.parent = MagicMock(return_value=MagicMock())
        with pytest.raises(ValueError, match='Not a subfolder'):
            # pylint: disable-next=protected-access       # OK in tests
            history._append_existing_folder(mock_path)

    def test_append_existing_folder_valid_path(self, history, mock_path):
        """Test _append_existing_folder with valid path."""
        history.collect_subfolders()
        mock_path.parent = history.path
        with patch_folder as mock_folder:
            mock_folder.return_value = folder = MagicMock()
            # pylint: disable-next=protected-access       # OK in tests
            history._append_existing_folder(mock_path)
            # pylint: disable-next=protected-access       # OK in tests
            assert folder in history._subfolders

    @parametrize(job_exists=(True, False))
    def test_find_name_for_new_history_subfolder_new(self, job_exists,
                                                     history):
        """Test _find_name_for_new_history_subfolder when no job exists."""
        # pylint: disable-next=protected-access           # OK in tests
        history._maps['jobs_for_tensor'] = {1: {1} if job_exists else {}}
        mock_concat = MagicMock(return_value=MagicMock(spec=Path))
        mock_concat.return_value.is_dir.return_value = job_exists
        history.path.__truediv__ = mock_concat
        args = 1, 'test'
        # pylint: disable-next=protected-access           # OK in tests
        result = history._find_name_for_new_history_subfolder(*args)
        assert result == f't00{args[0]}.r001_{args[1]}'

    def test_update_maps(self, history):
        """Test _update_maps method."""
        folders = MagicMock(), MagicMock()
        main = folders[0]
        main.parent = None
        # pylint: disable-next=protected-access           # OK in tests
        history._maps['main_hash_to_folders'] = {None: folders}
        # pylint: disable-next=protected-access           # OK in tests
        history._update_maps()
        expect = {f.hash_: main for f in folders}
        # pylint: disable-next=protected-access           # OK in tests
        assert history._maps['hash_to_parent'] == expect

    _too_early_attr = {
        'info': 'prepare_info_file',
        'new_folder': 'find_new_history_directory',
        'last_folder': None,
        'max_run_per_tensor': None,
        }
    _too_early_call = (
        'find_new_history_directory(None, None)',
        'register_folder',
        '_find_name_for_new_history_subfolder(None, None)',
        )

    @parametrize('attr,updater', _too_early_attr.items(), ids=_too_early_attr)
    def test_too_early_attr(self, attr, updater, history):
        """Check complaints when accessing attributes too early."""
        # pylint: disable-next=protected-access           # OK in tests
        updater = updater or HistoryExplorer._updater
        with pytest.raises(AttributeError, match=rf'.*{updater}.*'):
            attrgetter(attr)(history)

    @parametrize(method_name=_too_early_call)
    def test_too_early_method_call(self, method_name, history):
        """Check that accessing attributes before update_from_cwd fails."""
        if '(' not in method_name:
            args = tuple()
        else:
            method_name, args_str = method_name.split('(')
            args_str = args_str.replace(')', '') + ','
            args = ast.literal_eval(args_str)
        method = attrgetter(method_name)(history)
        with pytest.raises(AttributeError, match=rf'.*collect_subfolders.*'):
            method(*args)
