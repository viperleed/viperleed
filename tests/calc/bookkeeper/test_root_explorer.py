"""Tests for module root_explorer of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-16'
__license__ = 'GPLv3+'

from collections import namedtuple
from contextlib import contextmanager
from itertools import permutations
from operator import attrgetter
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.constants import STATE_FILES
from viperleed.calc.bookkeeper.root_explorer import LogFiles
from viperleed.calc.bookkeeper.root_explorer import RootExplorer
from viperleed.calc.sections.cleanup import DEFAULT_OUT
from viperleed.calc.sections.cleanup import DEFAULT_SUPP

from ...helpers import make_obj_raise
from ...helpers import not_raises
from ...helpers import raises_exception


_MOCK_ROOT = '/mock/root'
_MODULE = 'viperleed.calc.bookkeeper.root_explorer'
_TEST_STRING = 'this is a test value'
patch_discard = patch(f'{_MODULE}.discard_files')


@fixture(name='explorer')
def fixture_explorer():
    """Return a RootExplorer with a fake bookkeeper."""
    mock_bookkeeper = MagicMock()
    root_path = Path(_MOCK_ROOT)
    return RootExplorer(root_path, mock_bookkeeper)


@fixture(name='logs')
def fixture_logfiles():
    """Return a LogFiles instance."""
    return LogFiles(Path(_MOCK_ROOT))


@fixture(name='mock_path_readlines')
def fixture_mock_readlines(monkeypatch):
    """Temporarily replace open and readlines for pathlib.Path."""
    @contextmanager
    def _patch():
        mock_readlines = MagicMock(return_value=(_TEST_STRING,))
        mock_open_file = MagicMock(readlines=mock_readlines)
        mock_open = MagicMock()
        mock_open.return_value.__enter__.return_value = mock_open_file
        with monkeypatch.context() as patch_:
            patch_.setattr('pathlib.Path.open', mock_open)
            yield mock_open_file
    return _patch


@fixture(name='patched_path')
def fixture_patched_path(explorer, monkeypatch):
    """Return a version of pathlib.Path with fake read/write_text."""
    MockedPath = namedtuple('MockedPath', ('read', 'write', 'glob'))
    notes_file = explorer.path / 'notes.txt'
    with monkeypatch.context() as patch_:
        mock_read = MagicMock(return_value=_TEST_STRING)
        mock_write = MagicMock()
        mock_glob = MagicMock(return_value=iter([notes_file]))
        patch_.setattr('pathlib.Path.glob', mock_glob)
        patch_.setattr('pathlib.Path.read_text', mock_read)
        patch_.setattr('pathlib.Path.write_text', mock_write)
        yield MockedPath(mock_read, mock_write, mock_glob)


def mock_attributes(obj, *attrs):
    """Replace attributes of obj with MagicMock(s)."""
    for attr_name in attrs:
        setattr(obj, attr_name, MagicMock())


def called_or_not(condition):
    """Return the name of an assertion method depending on condition."""
    return 'assert_called' if condition else 'assert_not_called'


def check_methods_called(obj, method_name, **called):
    """Check that calling method_name also calls other methods."""
    mock_attributes(obj, *called)
    method = attrgetter(method_name)(obj)
    method()

    for mocked_attr, called_attr in called.items():
        method = attrgetter(called_attr or mocked_attr)(obj)
        method.assert_called_once()


class TestRootExplorer:
    """Tests for the RootExplorer class."""

    _recent = {
        'no log': None,
        'log file': _TEST_STRING,
        }

    @parametrize(expect=_recent.values(), ids=_recent)
    def test_calc_timestamp(self, explorer, expect):
        """Test calc_timestamp property."""
        most_recent = expect if expect is None else MagicMock(timestamp=expect)
        explorer._logs = MagicMock(most_recent=most_recent)
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
    def test_calls_methods(self, _, explorer, method_name, called):
        """Check calling of expected methods."""
        check_methods_called(explorer, method_name, **called)

    @patch(f'{_MODULE}.LogFiles')
    def test_collect_info(self, mock_logs, explorer):
        """Check calling of expected method when collecting info."""
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

    def test_collect_files_to_archive(self, explorer, monkeypatch):
        """Test the _collect_files_to_archive method."""
        explorer._logs = MagicMock(calc=(explorer.path / 'calc.log',))
        explorer.workhistory = MagicMock()
        explorer.workhistory.find_current_directories.return_value = (
            explorer.path / 'workhist_folder',
            )
        expected_files = (
            'OUT',
            'SUPP',
            'calc.log',
            'workhist_folder',
            )
        with monkeypatch.context() as patch_:
            patch_.setattr('pathlib.Path.is_file',
                           MagicMock(return_value=True))
            explorer._collect_files_to_archive()
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

    def test_infer_run_info(self, explorer):
        """Check correct result of inferring info from log files."""
        explorer._logs = MagicMock()
        explorer.logs.infer_run_info = MagicMock(return_value=_TEST_STRING)
        assert explorer.infer_run_info() is _TEST_STRING

    _archive_files = {
        False: (),
        True: ('file1', 'file2'),
        }

    @parametrize('expect,files', _archive_files.items(), ids=_archive_files)
    def test_needs_archiving(self, explorer, files, expect):
        """Check the needs_archiving property."""
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
                explorer._replace_state_files_from_ori()

    @parametrize(out_exists=(True, False))
    @parametrize(file_exists=(True, False))
    def test_update_state_files_from_out(self, out_exists, file_exists,
                                         explorer, monkeypatch):
        """Test moving of files from OUT."""
        mock_copy2 = MagicMock()
        mock_replace = MagicMock()
        with monkeypatch.context() as patch_:
            set_ = patch_.setattr
            set_('pathlib.Path.is_dir', MagicMock(return_value=out_exists))
            set_('pathlib.Path.is_file', MagicMock(return_value=file_exists))
            set_('shutil.copy2', mock_copy2)
            set_('pathlib.Path.replace', mock_replace)
            explorer.update_state_files_from_out()
        moved = out_exists and file_exists
        for mock_ in (mock_copy2, mock_replace):
            assert_ok = getattr(mock_, called_or_not(moved))
            assert_ok()


class TestRootExplorerRaises:
    """Tests for conditions that causes exceptions or complaints."""

    @patch(f'{_MODULE}.LOGGER.error')
    def test_cannot_replace_state_files_from_ori(self, mock_error, explorer):
        """Check complaints when a root file cannot be replaced with _ori."""
        with raises_exception('pathlib.Path.replace', OSError):
            explorer._replace_state_files_from_ori()
            mock_error.assert_called_once()

    _read_notes_raises = {
        'read': '',
        'write': _TEST_STRING,
        }

    @parametrize('method,expect',
                 _read_notes_raises.items(),
                 ids=_read_notes_raises)
    @patch(f'{_MODULE}.LOGGER.error')
    def test_read_and_clear_notes_file(self, mock_error,
                                       method, expect,
                                       explorer,
                                       patched_path):
        """Test that complaints by read/write_text emit errors."""
        complaining = getattr(patched_path, method)
        complaining.side_effect = OSError
        notes = explorer.read_and_clear_notes_file()
        assert notes == expect
        mock_error.assert_called_once()

    _attr_needs_update = (
        'calc_timestamp',
        'needs_archiving',
        )
    _method_needs_update = (  # Only those without args
        '_collect_files_to_archive',
        )

    @parametrize(attr=_attr_needs_update)
    def test_too_early_attribute_access(self, explorer, attr):
        """Check that accessing attributes before update_from_cwd fails."""
        with pytest.raises(AttributeError, match=rf'.*collect_info.*'):
            attrgetter(attr)(explorer)

    @parametrize(method_name=_method_needs_update)
    def test_too_early_method_call(self, explorer, method_name):
        """Check that accessing attributes before update_from_cwd fails."""
        method = attrgetter(method_name)(explorer)
        with pytest.raises(AttributeError, match=rf'.*collect_info.*'):
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

    def test_collect_calls(self, logs):
        """Check that calling .collect calls other methods."""
        called = {
            '_collect_logs': None,
            '_read_most_recent': None,
            }
        check_methods_called(logs, 'collect', **called)

    def test_collect_logs(self, logs):
        """Test that _collect_logs finds files in the root directory."""
        files = {
            '_calc': ('viperleed-calc_1.log', 'tleedm_1.log', 'tleedm_2.log'),
            '_others': ('other_stuff.log',)
            }
        for attr, file_list in files.items():
            files[attr] = tuple(logs._path / f for f in file_list)
        all_files = sum(files.values(), tuple())
        with patch('pathlib.Path.glob', return_value=all_files):
            with patch('pathlib.Path.is_file', return_value=True):
                logs._collect_logs()
        for attr, expect in files.items():
            assert getattr(logs, attr) == expect

    def test_collect_logs_no_logs(self, logs):
        """Test that _collect_logs finds files in the root directory."""
        files = {
            '_calc': ('viperleed-calc_1.log', 'tleedm_1.log', 'tleedm_2.log'),
            '_others': ('other_stuff.log',)
            }
        for attr, file_list in files.items():
            files[attr] = tuple(logs._path / f for f in file_list)
        all_files = sum(files.values(), tuple())
        with patch('pathlib.Path.glob', return_value=all_files):
            # None of the files above are actual files.
            logs._collect_logs()
        for attr, expect in files.items():
            collected = getattr(logs, attr)
            assert collected is not None
            assert not collected

    @patch_discard
    def test_discard(self, discard_files, logs):
        """Check deletion of all log files in root."""
        files = [logs._path / f'file{i}.log' for i in range(5)]
        logs._calc, logs._others = files[:2], files[2:]

        check_methods_called(logs, 'discard', collect=None)
        discard_files.assert_called_once_with(*files)

    def test_files_property(self, logs):
        """Test the files property of LogFiles."""
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
    def test_infer_run_info(self, logs, lines, expect):
        """Check correct inference of information from log files."""
        logs._calc = MagicMock()
        logs.most_recent = MagicMock()
        logs.most_recent.lines = lines
        assert logs.infer_run_info() == expect

    def test_infer_run_info_no_log(self, logs):
        """Check correct result when no log is found in root."""
        logs._calc = MagicMock()
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
        logs._calc = tuple(logs._path/f for f in collected)
        with mock_path_readlines():
            logs._read_most_recent()
        assert logs.most_recent == expect

    def test_most_recent_oserror(self, logs, mock_path_readlines):
        """Check that failure to read lines produces the expected outcome."""
        collected, (expect_time, _) = next(
            log for log in self._calc_logs.values() if log[1]
            )
        logs._calc = tuple(logs._path/f for f in collected)
        with mock_path_readlines() as mock_open:
            with make_obj_raise(mock_open, OSError, 'readlines'):
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
        with pytest.raises(AttributeError, match=rf'.*collect.*'):
            attrgetter(attr)(logs)

    @parametrize(method_name=_method_needs_update)
    def test_too_early_method_call(self, logs, method_name):
        """Check that accessing attributes before update_from_cwd fails."""
        method = attrgetter(method_name)(logs)
        with pytest.raises(AttributeError, match=rf'.*collect.*'):
            method()
