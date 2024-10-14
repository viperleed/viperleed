"""Tests for module workhistory of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from collections import defaultdict
import logging
from pathlib import Path
from unittest.mock import MagicMock
from unittest.mock import patch
from unittest.mock import call

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.sections.cleanup import PREVIOUS_LABEL
from viperleed.calc.bookkeeper.history.workhistory import WorkhistoryHandler

from ....helpers import make_obj_raise
from ....helpers import raises_exception


patch_rmtree = patch('shutil.rmtree')


@fixture(name='mock_bookkeeper')
def fixture_mock_bookkeeper():
    """Fixture to mock the bookkeeper."""
    bookkeeper = MagicMock()
    bookkeeper.top_level_history_path = Path('/mock/history')
    bookkeeper.timestamp = '20231008'
    jobs = defaultdict(int)
    jobs[1] = 5
    jobs[2] = 3
    bookkeeper.max_job_for_tensor = jobs
    return bookkeeper


@fixture(name='workhistory')
def fixture_workhistory(mock_path, mock_bookkeeper):
    """Fixture for WorkhistoryHandler instance."""
    return WorkhistoryHandler(mock_path, mock_bookkeeper)


@fixture(name='patched_path')
def factory_patched_path(mock_path):
    """Return a fake Path with methods patched as specified."""
    def _patch(**methods_and_return):
        for method_name, return_value in methods_and_return.items():
            method = getattr(mock_path, method_name)
            setattr(method, 'return_value', return_value)
        return mock_path
    return _patch


class TestWorkhistoryHandler:
    """Tests for the WorkhistoryHandler class."""

    def test_init(self, workhistory, mock_path, mock_bookkeeper):
        """Test initialization."""
        assert workhistory.path is mock_path
        assert workhistory.bookkeeper is mock_bookkeeper

    def test_history_property(self, workhistory, mock_bookkeeper):
        """Test history property."""
        assert workhistory.history is mock_bookkeeper.top_level_history_path

    def test_timestamp_property(self, workhistory, mock_bookkeeper):
        """Test timestamp property."""
        assert workhistory.timestamp is mock_bookkeeper.timestamp

    def test_find_current_directories(self, workhistory):
        """Test find_current_directories."""
        directories = [MagicMock(), MagicMock()]
        mock_find = patch.object(workhistory,
                                 '_find_directories',
                                 return_value=directories)
        with mock_find:
            subfolders = workhistory.find_current_directories()
            assert all(PREVIOUS_LABEL not in d.name for d in subfolders)

    _nothing_to_move = {
        'empty workhistory': tuple(),
        'one folder': (1,),
        }

    @patch_rmtree
    @parametrize(folders=_nothing_to_move.values(), ids=_nothing_to_move)
    def test_discard_workhistory_root(self, mock_rmtree, folders,
                                      workhistory, patched_path):
        """Test discard_workhistory_root."""
        mock_path = patched_path(is_dir=True, iterdir=iter(folders))
        workhistory.discard_workhistory_root()
        mock_rmtree.assert_called_once_with(mock_path)

    @patch_rmtree
    def test_discard_workhistory_root_not_there(self, mock_rmtree,
                                                workhistory, patched_path):
        """Test discard_workhistory_root when no folder is present."""
        patched_path(is_dir=False)
        workhistory.discard_workhistory_root()
        mock_rmtree.assert_not_called()

    @patch_rmtree
    def test_move_and_cleanup_nonexistent_folder(self, mock_rmtree,
                                                 workhistory, patched_path):
        """Test move_and_cleanup when the folder doesn't exist."""
        patched_path(is_dir=False)
        tensor_nums = workhistory.move_and_cleanup(None)
        assert not tensor_nums
        mock_rmtree.assert_not_called()

    @patch_rmtree
    def test_move_and_cleanup_with_folders(self, mock_rmtree,
                                           workhistory, patched_path):
        """Test move_and_cleanup where folders are moved."""
        tensors_expect = {1, 2}
        patched_path(is_dir=True, iterdir=iter([1]))
        mock_move = patch.object(workhistory,
                                 '_move_folders_to_history',
                                 return_value=tensors_expect)
        with mock_move:
            tensor_nums = workhistory.move_and_cleanup(None)
            assert tensor_nums == tensors_expect
            mock_rmtree.assert_not_called()

    @patch_rmtree
    def test_discard_previous(self, mock_rmtree, workhistory):
        """Test _discard_previous method."""
        previous_dirs = [MagicMock(), MagicMock()]
        mock_find = patch.object(workhistory,
                                 '_find_directories',
                                 return_value=previous_dirs)
        with mock_find:
            # pylint: disable-next=protected-access       # OK in tests
            workhistory._discard_previous()
            calls = [call(directory) for directory in previous_dirs]
            mock_rmtree.assert_has_calls(calls, any_order=True)

    def test_move_folders_to_history(self, workhistory, monkeypatch):
        """Test _move_folders_to_history."""
        directory = MagicMock()
        directory.name = f't987.r001.003_{workhistory.timestamp}'
        fake_replace = MagicMock()
        mock_replace = patch.object(directory, 'replace', fake_replace)
        mock_find = patch.object(workhistory,
                                 'find_current_directories',
                                 return_value=(directory,))
        mock_meta_cls = MagicMock()
        monkeypatch.setattr(
            'viperleed.calc.bookkeeper.history.workhistory.BookkeeperMetaFile',
            mock_meta_cls,
            )
        main_meta = MagicMock()
        with mock_find, mock_replace:
            # pylint: disable-next=protected-access       # OK in tests
            tensor_nums = workhistory._move_folders_to_history(main_meta)
            fake_replace.assert_called_once()
            assert tensor_nums == {987}
            mock_meta_cls.assert_called_once()  # Created metadata file
            mock_meta = mock_meta_cls.return_value
            # Metadata file processed and written
            mock_meta.compute_hash.assert_called_once()
            mock_meta.collect_from.assert_called_once_with(main_meta)
            mock_meta.write.assert_called_once()

    def test_find_any_directory(self, workhistory, patched_path):
        """Test _find_directories."""
        directories = [MagicMock(), MagicMock()]
        for directory in directories:
            directory.name = 'does not match'
        mock_path = patched_path(glob=directories, iterdir=directories)
        # pylint: disable-next=protected-access           # OK in tests
        result = list(workhistory._find_directories())
        assert not result  # No matching name
        mock_path.glob.assert_not_called()
        mock_path.iterdir.assert_called_once()

    def test_find_directories_with_filter(self, workhistory, patched_path):
        """Test _find_directories."""
        directories = [MagicMock(), MagicMock()]
        contains = '_some_test_string_'
        for i, directory in enumerate(directories):
            directory.name = f't123.r456_{contains}_{i}'
        mock_path = patched_path(glob=directories)
        # pylint: disable-next=protected-access           # OK in tests
        result = list(workhistory._find_directories(contains=contains))
        assert result == directories
        mock_path.glob.assert_called_once_with(f'*{contains}*')


class TestWorkhistoryHandlerRaises:
    """Tests for conditions that raise or log."""

    @staticmethod
    def check_has_error(logger):
        """Ensure logger has at least one ERROR message."""
        assert any(r for r in logger.records if r.levelno == logging.ERROR)

    def test_move_and_cleanup(self, workhistory, patched_path, caplog):
        """Check that failed removal of the workhistory logs errors."""
        patched_path(iterdir=tuple())
        # Make sure we don't raise OSError at _discard_previous
        mock_previous = patch.object(workhistory, '_discard_previous')
        with mock_previous, make_obj_raise('shutil.rmtree', OSError):
            workhistory.move_and_cleanup(None)
        self.check_has_error(caplog)

    def test_discard_previous(self, workhistory, caplog):
        """Check that failing to remove "previous" does not break."""
        mock_find = patch.object(workhistory,
                                 '_find_directories',
                                 return_value=(1,))
        with mock_find, make_obj_raise('shutil.rmtree', OSError):
            # pylint: disable-next=protected-access       # OK in tests
            workhistory._discard_previous()
        self.check_has_error(caplog)

    def test_move_folders_exists(self, workhistory, caplog):
        """Test _move_folders_to_history handling FileExistsError."""
        directory = MagicMock()
        directory.name = f't001.r001.003_{workhistory.timestamp}'
        mock_find = patch.object(workhistory,
                                 'find_current_directories',
                                 return_value=(directory,))
        raises_ = raises_exception(directory, FileExistsError, 'replace')
        with patch_rmtree, mock_find, raises_:
            # pylint: disable-next=protected-access       # OK in tests
            workhistory._move_folders_to_history(None)
        self.check_has_error(caplog)

    def test_move_folders_os_error(self, workhistory, caplog):
        """Test _move_folders_to_history handling OSError."""
        directory = MagicMock()
        directory.name = f't001.r001.003_{workhistory.timestamp}'
        mock_find = patch.object(workhistory,
                                 'find_current_directories',
                                 return_value=(directory,))
        raises_ = make_obj_raise(directory, OSError, 'replace')
        with mock_find, raises_:
            # pylint: disable-next=protected-access       # OK in tests
            tensor_nums = workhistory._move_folders_to_history(None)
            assert not tensor_nums
        self.check_has_error(caplog)
