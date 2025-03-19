"""Tests for module log of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-09-29'
__license__ = 'GPLv3+'

from contextlib import contextmanager
from itertools import permutations
import logging
from operator import attrgetter
from pathlib import Path

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.log import LogFiles
from viperleed.calc.bookkeeper.log import add_bookkeeper_logfile
from viperleed.calc.bookkeeper.log import ensure_has_stream_handler

from ...helpers import make_obj_raise


_MODULE = 'viperleed.calc.bookkeeper.log'
_TEST_STRING = 'this is a test value'


class MockLogger:
    """A fake logging.Logger."""

    def __init__(self):
        """Initialize instance."""
        # We mimic the attributes that we need access to for the
        # functionality in bookkeeper.log and calc.lib.log_utils
        self.handlers = []
        self.parent = None
        self.propagate = False
        self.level = logging.NOTSET

    def addHandler(self, handler):       # pylint: disable=invalid-name
        """Add one handler."""
        self.handlers.append(handler)

    def setLevel(self, level):           # pylint: disable=invalid-name
        """set a new logging level."""
        self.level = level


@fixture(name='logs')
def fixture_logfiles():
    """Return a LogFiles instance."""
    return LogFiles(Path('/mock/root'))


@fixture(name='logger')
def fixture_logger(monkeypatch):
    """Replace the bookkeeper logger with a fake one."""
    with monkeypatch.context() as patch:
        logger = MockLogger()
        patch.setattr('viperleed.calc.bookkeeper.log.LOGGER', logger)
        yield logger
        for handler in logger.handlers:
            handler.close()


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


class TestAddBookkeeperLogfile:
    """Tests for the add_bookkeeper_logfile function."""

    def test_already_there(self, logger, tmp_path):
        """Check that the same file handler is not added twice."""
        add_bookkeeper_logfile(tmp_path)
        handlers_before = logger.handlers.copy()
        add_bookkeeper_logfile(tmp_path)
        assert logger.handlers == handlers_before

    def test_not_there(self, logger, tmp_path):
        """Check addition of a file handler to BOOKIE_LOGFILE."""
        assert not logger.handlers
        add_bookkeeper_logfile(tmp_path)
        assert logger.handlers


class TestEnsureHasStreamHandler:
    """Tests for the ensure_has_stream_handler function."""

    def test_already_there(self, logger):
        """Check no repeated addition of the stdout stream handler."""
        ensure_has_stream_handler()
        handlers_before = logger.handlers.copy()
        ensure_has_stream_handler()
        assert logger.handlers == handlers_before

    def test_not_there(self, logger):
        """Check addition of a console handler to stdout."""
        assert not logger.handlers
        ensure_has_stream_handler()
        assert logger.handlers
        assert logger.propagate
        assert logger.level == logging.INFO


def _combine_log_info(**log_info):
    """Return all combinations of lines and expected values for `log_info`."""
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
            '_infer_calc_version': None,
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
        """Test that _collect_logs finds no existing files."""
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
            assert not any(collected)

    def test_discard(self, logs, check_methods_called, mocker):
        """Check deletion of all log files in root."""
        discard_files = mocker.patch(f'{_MODULE}.discard_files')
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

    _calc_version = {
        'no lines': ((), None),
        'one line': (('This is ViPErLEED version 1.2.3',), '1.2.3'),
        'more lines': (
            ('First line',
             '',
             'This is ViPErLEED version <but this is not a version>',
             'This is ViPErLEED version 0.1.0',
             'Another line',
             'This is ViPErLEED version 5.12.103'),  # Not used
            '0.1.0'),
        }

    @parametrize('lines,expect', _calc_version.values(), ids=_calc_version)
    def test_infer_calc_version(self, lines, expect, logs, mocker):
        """Test the result of the _infer_calc_version method."""
        mocker.patch.object(logs, '_calc')
        logs.most_recent = mocker.MagicMock(lines=lines)
        # pylint: disable-next=protected-access           # OK in tests
        logs._infer_calc_version()
        assert logs.version == expect

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
        """Check that accessing attributes before collection fails."""
        with pytest.raises(AttributeError, match=r'.*collect.*'):
            attrgetter(attr)(logs)

    @parametrize(method_name=_method_needs_update)
    def test_too_early_method_call(self, logs, method_name):
        """Check that calling methods before collection fails."""
        method = attrgetter(method_name)(logs)
        with pytest.raises(AttributeError, match=r'.*collect.*'):
            method()
