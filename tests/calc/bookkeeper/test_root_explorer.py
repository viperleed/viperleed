"""Tests for module root_explorer of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-16'
__license__ = 'GPLv3+'

from collections import namedtuple
from contextlib import contextmanager
from itertools import combinations
from itertools import permutations
from operator import attrgetter
from pathlib import Path
from unittest.mock import patch

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.constants import STATE_FILES
from viperleed.calc.bookkeeper.root_explorer import LogFiles
from viperleed.calc.bookkeeper.root_explorer import RootExplorer
from viperleed.calc.bookkeeper.root_explorer import TensorAndDeltaInfo
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_TENSORS

from ...helpers import make_obj_raise
from ...helpers import not_raises
from ...helpers import raises_exception


_MOCK_ROOT = '/mock/root'
_MODULE = 'viperleed.calc.bookkeeper.root_explorer'
_TEST_STRING = 'this is a test value'
patch_discard = patch(f'{_MODULE}.discard_files')


@fixture(name='explorer')
def fixture_explorer(mocker):
    """Return a RootExplorer with a fake bookkeeper."""
    mock_bookkeeper = mocker.MagicMock()
    root_path = Path(_MOCK_ROOT)
    return RootExplorer(root_path, mock_bookkeeper)


@fixture(name='logs')
def fixture_logfiles():
    """Return a LogFiles instance."""
    return LogFiles(Path(_MOCK_ROOT))


@fixture(name='tensors')
def fixture_tensor_and_delta_info():
    """Fixture to initialize TensorAndDeltaInfo with a mock root path."""
    return TensorAndDeltaInfo(Path('/fake/root'))


@fixture(name='mock_history')
def fixture_mock_history(mocker):
    """Return a fake HistoryExplorer."""
    def _mock(max_run_per_tensor):
        return mocker.MagicMock(
            max_run_per_tensor=max_run_per_tensor,
            last_folder_and_siblings=[mocker.MagicMock(tensor_num=t)
                                      for t in max_run_per_tensor],
            )
    return _mock


@fixture(name='mock_path_readlines')
def fixture_mock_readlines(mocker, monkeypatch):
    """Temporarily replace open and readlines for pathlib.Path."""
    @contextmanager
    def _patch():
        mock_readlines = mocker.MagicMock(return_value=(_TEST_STRING,))
        mock_open_file = mocker.MagicMock(readlines=mock_readlines)
        mock_open = mocker.MagicMock()
        mock_open.return_value.__enter__.return_value = mock_open_file
        with monkeypatch.context() as patch_:
            patch_.setattr('pathlib.Path.open', mock_open)
            yield mock_open_file
    return _patch


@fixture(name='patched_path')
def fixture_patched_path(explorer, mocker):
    """Return a version of pathlib.Path with fake read/write_text."""
    MockedPath = namedtuple('MockedPath', ('read', 'write', 'glob'))
    notes_file = explorer.path / 'notes.txt'
    mock_glob = mocker.patch('pathlib.Path.glob',
                             return_value=iter([notes_file]))
    mock_read = mocker.patch('pathlib.Path.read_text',
                             return_value=_TEST_STRING)
    mock_write = mocker.patch('pathlib.Path.write_text')
    yield MockedPath(mock_read, mock_write, mock_glob)


@fixture(name='mock_attributes')
def fixture_mock_attributes(mocker):
    """Replace attributes of obj with MagicMock(s)."""
    def _patch(obj, *attrs):
        for attr_name in attrs:
            mocker.patch.object(obj, attr_name)
    return _patch


def called_or_not(condition):
    """Return the name of an assertion method depending on condition."""
    return 'assert_called' if condition else 'assert_not_called'


@fixture(name='check_methods_called')
def fixture_check_methods_called(mock_attributes):
    """Check that calling method_name also calls other methods."""
    def check_methods_called(obj, method_name, **called):
        mock_attributes(obj, *called)
        method = attrgetter(method_name)(obj)
        method()

        for mocked_attr, called_attr in called.items():
            method = attrgetter(called_attr or mocked_attr)(obj)
            method.assert_called_once()
    return check_methods_called


def mock_exists(collection):
    """Return a callable that can replace pathlib.Path.exists."""
    def _exists(path):
        return path in collection
    return _exists


class TestRootExplorer:
    """Tests for the RootExplorer class."""

    _recent = {
        'no log': None,
        'log file': _TEST_STRING,
        }

    @parametrize(expect=_recent.values(), ids=_recent)
    def test_calc_timestamp(self, explorer, expect, mocker):
        """Test calc_timestamp property."""
        most_recent = (expect if expect is None
                       else mocker.MagicMock(timestamp=expect))
        # pylint: disable-next=protected-access           # OK in tests
        explorer._logs = mocker.MagicMock(most_recent=most_recent)
        assert explorer.calc_timestamp == expect

    _called = {
        # outer_call: {mocked_attr: inner_call or None}
        # if None, take mocked_attr as the one that is called
        'clear_for_next_calc_run' : {
            '_logs': 'logs.discard',
            '_remove_out_and_supp': None,
            '_remove_ori_files': None,
            'collect_info': None,
            },
        'revert_to_previous_calc_run' : {
            '_logs': 'logs.discard',
            '_remove_out_and_supp': None,
            '_replace_state_files_from_ori': None,
            'collect_info': None,
            },
        }

    @parametrize('method_name,called', _called.items(), ids=_called)
    @patch_discard
    def test_calls_methods(self, _, explorer, method_name, called,
                           check_methods_called):
        """Check calling of expected methods."""
        check_methods_called(explorer, method_name, **called)

    def test_collect_info(self, explorer, mock_attributes, mocker):
        """Check calling of expected method when collecting info."""
        mock_logs = mocker.patch(f'{_MODULE}.LogFiles')
        mock_attributes(mock_logs, 'collect')
        mock_attributes(explorer, '_collect_files_to_archive')

        explorer.collect_info()

        explorer.logs.collect.assert_called_once()
        called = {'_collect_files_to_archive', 'collect_info'}
        for method in dir(explorer):
            try:
                called_ok = getattr(explorer, called_or_not(method in called))
            except AttributeError:
                pass
            else:
                called_ok()

    def test_collect_files_to_archive(self, explorer, mocker):
        """Test the _collect_files_to_archive method."""
        # pylint: disable-next=protected-access           # OK in tests
        explorer._logs = mocker.MagicMock(calc=(explorer.path / 'calc.log',))
        mocker.patch.object(explorer, 'workhistory')
        explorer.workhistory.find_current_directories.return_value = (
            explorer.path / 'workhist_folder',
            )
        expected_files = (
            DEFAULT_OUT,
            DEFAULT_SUPP,
            'calc.log',
            'workhist_folder',
            )
        mocker.patch('pathlib.Path.is_file', return_value=True)
        # pylint: disable-next=protected-access           # OK in tests
        explorer._collect_files_to_archive()
        # pylint: disable-next=protected-access           # OK in tests
        to_archive = explorer._files_to_archive
        assert to_archive == tuple(explorer.path / f for f in expected_files)


    _remove_files = {
        '_remove_ori_files': [f'{file}_ori' for file in STATE_FILES],
        '_remove_out_and_supp': (DEFAULT_OUT, DEFAULT_SUPP),
        }

    @parametrize('method_name,files', _remove_files.items(), ids=_remove_files)
    @patch_discard
    def test_files_discarded(self, discard_files,
                             method_name, files,
                             explorer):
        """Check that files and folders are removed."""
        method = getattr(explorer, method_name)
        method()
        removed = (explorer.path / f for f in files)
        discard_files.assert_called_once_with(*removed)

    def test_infer_run_info(self, explorer, mocker):
        """Check correct result of inferring info from log files."""
        mocker.patch.object(explorer, '_logs')
        mocker.patch.object(explorer.logs,
                            'infer_run_info',
                            return_value=_TEST_STRING)
        assert explorer.infer_run_info() is _TEST_STRING

    _archive_files = {
        False: (),
        True: ('file1', 'file2'),
        }

    @parametrize('expect,files', _archive_files.items(), ids=_archive_files)
    def test_needs_archiving(self, explorer, files, expect):
        """Check the needs_archiving property."""
        # pylint: disable-next=protected-access           # OK in tests
        explorer._files_to_archive = files
        assert explorer.needs_archiving == expect

    def test_read_and_clear_notes_file(self, explorer, patched_path):
        """Test successful reading and clearing of a notes file."""
        notes = explorer.read_and_clear_notes_file()
        assert notes == _TEST_STRING
        patched_path.read.assert_called_once_with(encoding='utf-8')
        patched_path.write.assert_called_once_with('', encoding='utf-8')

    def test_read_and_clear_notes_file_no_notes(self, explorer, patched_path):
        """Test successful reading and clearing of a notes file."""
        patched_path.glob.return_value = iter(())
        notes = explorer.read_and_clear_notes_file()
        assert not notes
        patched_path.read.assert_not_called()
        patched_path.write.assert_not_called()

    def test_replace_state_files_from_ori_no_ori(self, explorer):
        """Check there are no complaints if an _ori file is not there."""
        with make_obj_raise('pathlib.Path.replace', FileNotFoundError):
            with not_raises(FileNotFoundError):
                # pylint: disable-next=protected-access   # OK in tests
                explorer._replace_state_files_from_ori()

    @parametrize(out_exists=(True, False))
    @parametrize(file_exists=(True, False))
    def test_update_state_files_from_out(self, out_exists, file_exists,
                                         explorer, mocker):
        """Test moving of files from OUT."""
        mock_copy2 = mocker.patch('shutil.copy2')
        mock_replace = mocker.patch('pathlib.Path.replace')
        mocker.patch('pathlib.Path.is_dir', return_value=out_exists)
        mocker.patch('pathlib.Path.is_file', return_value=file_exists)
        explorer.update_state_files_from_out()
        moved = out_exists and file_exists
        for mock_ in (mock_copy2, mock_replace):
            assert_ok = getattr(mock_, called_or_not(moved))
            assert_ok()


class TestRootExplorerListFiles:
    """Tests for listing files/directories to be removed."""

    def test_to_discard_no_files(self, explorer, mocker):
        """Test list_paths_to_discard when no discardable files exist."""
        explorer.collect_info()
        mocker.patch.object(Path, 'exists', return_value=False)
        paths_to_discard = explorer.list_paths_to_discard()
        assert not paths_to_discard

    def test_to_discard_files_exist(self, explorer, mocker, mock_history):
        """Test list_paths_to_discard when some discardable files exist."""
        mock_log_path = explorer.path / 'log_file.log'
        mock_files = (
            mock_log_path,
            explorer.path / DEFAULT_TENSORS / f'{DEFAULT_TENSORS}_001.zip',
            )
        mock_paths = (
            explorer.path / DEFAULT_OUT,
            # explorer.path / DEFAULT_SUPP,  # No SUPP
            *mock_files
            )
        explorer.collect_info()
        mocker.patch.object(Path, 'exists', new=mock_exists(mock_paths))
        mocker.patch.object(Path, 'is_file', new=mock_exists(mock_files))
        mocker.patch.object(type(explorer.logs), 'files', [mock_log_path])
        mocker.patch.object(explorer, 'history', mock_history({1: 1}))
        paths_to_discard = explorer.list_paths_to_discard()
        assert paths_to_discard == mock_paths

    def test_to_replace_no_ori(self, explorer, mocker):
        """Test list_files_to_replace when no '_ori' files are present."""
        mocker.patch.object(Path, 'exists', return_value=False)
        files_to_replace = explorer.list_files_to_replace()
        assert not files_to_replace

    _ori_files = {}
    for repeats in range(len(STATE_FILES)):
        _ori_files.update({'+'.join(files): files
                           for files in combinations(STATE_FILES, repeats+1)})

    @parametrize(state_files=_ori_files.values(), ids=_ori_files)
    def test_to_replace_ori_files_exist(self, state_files, explorer, mocker):
        """Test list_files_to_replace when '_ori' files are present."""
        mock_ori = [explorer.path / f'{f}_ori' for f in state_files]
        mocker.patch.object(Path, 'exists', new=mock_exists(mock_ori))
        files_to_replace = explorer.list_files_to_replace()
        expected_files_to_replace = {
            explorer.path / f: ori_f
            for f, ori_f in zip(state_files, mock_ori)
            }
        assert dict(files_to_replace) == expected_files_to_replace


class TestRootExplorerRaises:
    """Tests for conditions that causes exceptions or complaints."""

    def test_cannot_replace_state_files_from_ori(self, explorer, mocker):
        """Check complaints when a root file cannot be replaced with _ori."""
        mock_error = mocker.patch(f'{_MODULE}.LOGGER.error')
        with raises_exception('pathlib.Path.replace', OSError):
            mocker.patch('pathlib.Path.exists', return_value=True)
            # pylint: disable-next=protected-access       # OK in tests
            explorer._replace_state_files_from_ori()
            mock_error.assert_called_once()

    _read_notes_raises = {
        'read': '',
        'write': _TEST_STRING,
        }

    @parametrize('method,expect',
                 _read_notes_raises.items(),
                 ids=_read_notes_raises)
    def test_read_and_clear_notes_file(self, method, expect, explorer,
                                       patched_path, mocker):
        """Test that complaints by read/write_text emit errors."""
        complaining = getattr(patched_path, method)
        complaining.side_effect = OSError
        mock_error = mocker.patch(f'{_MODULE}.LOGGER.error')
        notes = explorer.read_and_clear_notes_file()
        assert notes == expect
        mock_error.assert_called_once()

    _attr_needs_update = (
        'calc_timestamp',
        'needs_archiving',
        )
    _method_needs_update = (  # Only those without args
        'clear_for_next_calc_run',
        'infer_run_info',
        'list_paths_to_discard',
        'remove_tensors_and_deltas',
        'revert_to_previous_calc_run',
        '_collect_files_to_archive',
        )

    @parametrize(attr=_attr_needs_update)
    def test_too_early_attribute_access(self, explorer, attr):
        """Check that accessing attributes before update_from_cwd fails."""
        with pytest.raises(AttributeError,
                           match=r'.*collect_(info|subfolders).*'):
            attrgetter(attr)(explorer)

    @parametrize(method_name=_method_needs_update)
    def test_too_early_method_call(self, explorer, method_name):
        """Check that accessing attributes before update_from_cwd fails."""
        method = attrgetter(method_name)(explorer)
        with pytest.raises(AttributeError,
                           match=r'.*collect_(info|subfolders).*'):
            method()


def _combine_log_info(**log_info):
    """Return all combinations of lines and expected values for log_info."""
    _info_keys = tuple(log_info)
    for _repeat in range(2, len(_info_keys)+1):
        for keys in permutations(_info_keys, _repeat):
            all_lines, all_expect = zip(*(log_info[k] for k in keys))
            lines = sum(all_lines, tuple())
            expect = {k: v for dict_ in all_expect for k, v in dict_.items()}
            log_info['+'.join(keys)] = (lines, expect)
    return log_info


class TestLogFiles:
    """Tests for the LogFiles class."""

    def test_collect_calls(self, logs, check_methods_called):
        """Check that calling .collect calls other methods."""
        called = {
            '_collect_logs': None,
            '_read_most_recent': None,
            }
        check_methods_called(logs, 'collect', **called)

    def test_collect_logs(self, logs, mocker):
        """Test that _collect_logs finds files in the root directory."""
        files = {
            '_calc': ('viperleed-calc_1.log', 'tleedm_1.log', 'tleedm_2.log'),
            '_others': ('other_stuff.log',)
            }
        for attr, file_list in files.items():
            # pylint: disable-next=protected-access       # OK in tests
            files[attr] = tuple(logs._path / f for f in file_list)
        mocker.patch('pathlib.Path.glob',
                     return_value=sum(files.values(), tuple()))
        mocker.patch('pathlib.Path.is_file', return_value=True)
        # pylint: disable-next=protected-access           # OK in tests
        logs._collect_logs()
        for attr, expect in files.items():
            assert getattr(logs, attr) == expect

    def test_collect_logs_no_logs(self, logs, mocker):
        """Test that _collect_logs finds files in the root directory."""
        files = {
            '_calc': ('viperleed-calc_1.log', 'tleedm_1.log', 'tleedm_2.log'),
            '_others': ('other_stuff.log',)
            }
        for attr, file_list in files.items():
            # pylint: disable-next=protected-access       # OK in tests
            files[attr] = tuple(logs._path / f for f in file_list)
        mocker.patch('pathlib.Path.glob',
                     return_value=sum(files.values(), tuple()))
        # pylint: disable-next=protected-access           # OK in tests
        logs._collect_logs()
        for attr in files:
            collected = getattr(logs, attr)
            assert collected is not None
            # None of the files above are actual files.
            assert not collected

    @patch_discard
    def test_discard(self, discard_files, logs, check_methods_called):
        """Check deletion of all log files in root."""
        # pylint: disable-next=protected-access           # OK in tests
        files = [logs._path / f'file{i}.log' for i in range(5)]
        logs._calc, logs._others = files[:2], files[2:]

        check_methods_called(logs, 'discard', collect=None)
        discard_files.assert_called_once_with(*files)

    def test_files_property(self, logs):
        """Test the files property of LogFiles."""
        # pylint: disable-next=protected-access           # OK in tests
        files = tuple(logs._path / f
                      for f in ('calc1.log', 'calc2.log', 'other.log'))
        logs._calc, logs._others = files[:2], files[2:]
        assert logs.files == files

    _log_info = _combine_log_info(
        run=(('Executed segments: 123   ',), {'run_info': '123   '}),
        r_ref=(('Final R (refcalc)   : 0.01 ',), {'r_ref': '0.01 '}),
        r_super=(('Final R (superpos)   : 0.02  ',), {'r_super': '0.02  '}),
        )

    @parametrize('lines,expect', _log_info.values(), ids=_log_info)
    def test_infer_run_info(self, logs, lines, expect, mocker):
        """Check correct inference of information from log files."""
        mocker.patch.object(logs, '_calc')
        logs.most_recent = mocker.MagicMock(lines=lines)
        assert logs.infer_run_info() == expect

    def test_infer_run_info_no_log(self, logs, mocker):
        """Check correct result when no log is found in root."""
        # pylint: disable-next=protected-access           # OK in tests
        logs._calc = mocker.MagicMock()
        assert not logs.infer_run_info()

    _calc_logs = {
        'no logs': ((), None),
        'calc, one log': (('viperleed-calc_123456-123456.log',),
                          ('123456-123456', (_TEST_STRING,))),
        'tleedm, one log': (('tleedm_123456-123456.log',),
                            ('123456-123456', (_TEST_STRING,))),
        'calc, more logs': (
            ('viperleed-calc_123456-123456.log',
             'viperleed-calc_234567-234567.log',
             'viperleed-calc_345678-345678.log',
             'viperleed-calc_456789-456789.log'),
            ('456789-456789', (_TEST_STRING,))
            ),
        'mixed logs': (
            ('viperleed-calc_123456-123456.log',
             'tleedm_234567-234567.log',
             'viperleed-calc_345678-345678.log',
             'tleedm_456789-456789.log'),
            ('345678-345678', (_TEST_STRING,))
            ),
        }

    @parametrize('collected,expect', _calc_logs.values(), ids=_calc_logs)
    def test_most_recent(self, logs, collected, expect, mock_path_readlines):
        """Test the _read_most_recent method of LogFiles."""
        # pylint: disable-next=protected-access           # OK in tests
        logs._calc = tuple(logs._path/f for f in collected)
        with mock_path_readlines():
            # pylint: disable-next=protected-access       # OK in tests
            logs._read_most_recent()
        assert logs.most_recent == expect

    def test_most_recent_oserror(self, logs, mock_path_readlines):
        """Check that failure to read lines produces the expected outcome."""
        collected, (expect_time, _) = next(
            log for log in self._calc_logs.values() if log[1]
            )
        # pylint: disable-next=protected-access           # OK in tests
        logs._calc = tuple(logs._path/f for f in collected)
        with mock_path_readlines() as mock_open:
            with make_obj_raise(mock_open, OSError, 'readlines'):
                # pylint: disable-next=protected-access   # OK in tests
                logs._read_most_recent()
            assert logs.most_recent == (expect_time, ())

    _attr_needs_update = (
        'files',
        )
    _method_needs_update = (  # Only those without args
        'infer_run_info',
        '_read_most_recent',
        )

    @parametrize(attr=_attr_needs_update)
    def test_too_early_attribute_access(self, logs, attr):
        """Check that accessing attributes before update_from_cwd fails."""
        with pytest.raises(AttributeError, match=r'.*collect.*'):
            attrgetter(attr)(logs)

    @parametrize(method_name=_method_needs_update)
    def test_too_early_method_call(self, logs, method_name):
        """Check that accessing attributes before update_from_cwd fails."""
        method = attrgetter(method_name)(logs)
        with pytest.raises(AttributeError, match=r'.*collect.*'):
            method()


class TestTensorsAndDeltasInfo:
    """Tests for the TensorsAndDeltasInfo class."""

    tensor_index = 5  # Example index to use for Tensor/Delta files

    @fixture(name='mock_get_tensors', autouse=True)
    def patch_get_tensors(self, mocker):
        """Replace getMaxTensorIndex with a fake one."""
        return mocker.patch(f'{_MODULE}.getMaxTensorIndex',
                            return_value=self.tensor_index)

    @fixture(name='simple_history')
    def fixture_simple_history(self, mock_history):
        """Return a fake HistoryExplorer with one tensor and one run."""
        return mock_history({self.tensor_index: 1})

    def test_collect(self, tensors, mock_get_tensors):
        """Test that collect sets the correct Tensor index."""
        tensors.collect()
        # pylint: disable-next=protected-access       # OK in tests
        mock_get_tensors.assert_called_once_with(home=tensors._root, zip_only=True)
        assert tensors.most_recent == self.tensor_index

    def test_to_discard(self, tensors, simple_history, mocker):
        """Test list_paths_to_discard when Tensor and Delta files exist."""
        expected = (
            f'{DEFAULT_DELTAS}/{DEFAULT_DELTAS}_{self.tensor_index:03d}.zip',
            f'{DEFAULT_TENSORS}/{DEFAULT_TENSORS}_{self.tensor_index:03d}.zip',
            )
        # pylint: disable-next=protected-access           # OK in tests
        expected = tuple(tensors._root / f for f in expected)
        tensors.collect()
        mocker.patch.object(Path, 'is_file', new=mock_exists(expected))
        discard_paths = tensors.list_paths_to_discard(simple_history)
        assert discard_paths == expected

    def test_to_discard_none(self, tensors, simple_history, mocker):
        """Test list_paths_to_discard when no Tensor/Delta files exist."""
        tensors.collect()
        mocker.patch.object(Path, 'is_file', return_value=False)
        assert not tensors.list_paths_to_discard(simple_history)

    def test_to_discard_history_empty(self, tensors, mock_history, mocker):
        """Check that no Tensors/Deltas are discarded with empty history."""
        tensors.collect()
        mocker.patch.object(Path, 'is_file', return_value=True)
        history = mock_history({})
        assert not tensors.list_paths_to_discard(history)

    def test_to_discard_older_runs(self, tensors, mock_history, mocker):
        """Check that no Tensors/Deltas are discarded when more runs exist."""
        tensors.collect()
        mocker.patch.object(Path, 'is_file', return_value=True)
        history = mock_history({self.tensor_index: 3})
        assert not tensors.list_paths_to_discard(history)

    @patch_discard
    def test_discard(self, mock_discard_files,
                     tensors, simple_history, mocker):
        """Test removal of Tensors and Deltas."""
        removed = (
            f'{DEFAULT_DELTAS}/{DEFAULT_DELTAS}_{self.tensor_index:03d}.zip',
            f'{DEFAULT_TENSORS}/{DEFAULT_TENSORS}_{self.tensor_index:03d}.zip',
            )
        # pylint: disable-next=protected-access           # OK in tests
        removed = tuple(tensors._root/f for f in removed)
        tensors.collect()
        mocker.patch.object(Path, 'is_file', return_value=True)
        tensors.discard(simple_history)
        mock_discard_files.assert_called_once_with(*removed)

    _attr_needs_update = (
        'most_recent',
        )
    _method_needs_update = (
        # None for now
        )

    @parametrize(attr=_attr_needs_update)
    def test_too_early_attribute_access(self, tensors, attr):
        """Check that accessing attributes before update_from_cwd fails."""
        with pytest.raises(AttributeError, match=r'.*collect.*'):
            attrgetter(attr)(tensors)

    @parametrize(method_name=_method_needs_update)
    def test_too_early_method_call(self, tensors, method_name):
        """Check that accessing attributes before update_from_cwd fails."""
        method = attrgetter(method_name)(tensors)
        with pytest.raises(AttributeError, match=r'.*collect.*'):
            method()
