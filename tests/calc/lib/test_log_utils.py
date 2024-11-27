"""Tests for the log_utils module of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-05'
__license__ = 'GPLv3+'

from contextlib import contextmanager
import logging
from logging import LogRecord as Record
from unittest.mock import Mock

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.lib.log_utils import at_level
from viperleed.calc.lib.log_utils import close_all_handlers
from viperleed.calc.lib.log_utils import debug_or_lower
from viperleed.calc.lib.log_utils import logging_silent
from viperleed.calc.lib.log_utils import logger_silent
from viperleed.calc.lib.log_utils import prepare_calc_logger
from viperleed.calc.lib.log_utils import CalcLogFormatter

from ...helpers import make_obj_raise
from ...helpers import not_raises
from ...helpers import raises_test_exception


@fixture(name='make_logger', scope='session')
def factory_make_logger():
    """Return a logger at a given level."""
    @contextmanager
    def _make(level=None, name='test_logger'):
        logger = logging.getLogger(name)
        if level is not None:
            logger.setLevel(level)
        try:
            yield logger
        finally:
            # logging uses the garbage collector to
            # clean up loggers. Let's do the same here.
            logger.manager.loggerDict.pop(name, None)
    return _make


class MockHandler(Mock):
    """A mock for the logging.FileHandler class."""

    acquire = Mock()
    flush = Mock()
    close = Mock()
    release = Mock()
    setFormatter = Mock()  # pylint: disable=invalid-name


class TestAtLevel:
    """Tests for the at_level function."""

    def test_level_changed(self, make_logger):
        """Check that the logging level is temporarily changed."""
        with make_logger(logging.WARNING) as logger:
            assert logger.level == logging.WARNING
            with at_level(logger, logging.DEBUG):
                assert logger.level == logging.DEBUG
            assert logger.level == logging.WARNING

    def test_messages_silenced(self, make_logger, caplog):
        """Check that at_level silences messages."""
        with make_logger(logging.DEBUG) as logger:
            msg = 'Warning that should be silenced'
            with at_level(logger, logging.ERROR):
                logger.warning(msg)
            assert msg not in caplog.text


# About the pylint disable: while it is true that this class only has
# one method, the not-so-simple parametrize seems a good-enough reason
# to keep it localized into a test class.
class TestCalcLogFormatter:  # pylint: disable=too-few-public-methods
    """Tests for the CalcLogFormatter class."""

    kwargs = {  # For LogRecord
        'name': 'test',
        'pathname': None,
        'lineno': 0,
        'msg': 'message',
        'args': None,
        'exc_info': None,
        }
    _records = {
        Record(level=1, **kwargs): 'dbg: message',
        Record(level=logging.DEBUG, **kwargs): 'dbg: message',
        Record(level=logging.INFO, **kwargs): 'message',
        Record(level=logging.WARNING, **kwargs): '# WARNING: message',
        Record(level=logging.ERROR, **kwargs): '''\
### ERROR ### in Unknown module:None:0
# message
#############''',
        Record(level=logging.CRITICAL, **kwargs): '''\
### CRITICAL ### in Unknown module:None:0
# message
################''',
        Record(level=48, **kwargs): 'message',
        }

    @parametrize('record,result', _records.items())
    def test_format(self, record, result):
        """Check expected outcome of formatting a record."""
        fmt = CalcLogFormatter()
        assert fmt.format(record) == result


class TestCloseAllHandlers:
    """Tests for the close_all_handlers function."""

    def test_nr_handlers(self, make_logger):
        """Check that the number of handlers is as expected."""
        with make_logger() as logger:
            handlers = logging.StreamHandler(), logging.StreamHandler()
            for handler in handlers:
                logger.addHandler(handler)
            assert len(logger.handlers) == len(handlers)
            close_all_handlers(logger)
            assert not logger.handlers

    _reason = ('Looks like this one is not working as intended. The '
               'problem is that each logger has at least the root '
               'logger as a parent, and caplog seems to hook into '
               'that one. Removing the parent or setting .propagate '
               'to False also removes all the messages captured by '
               'caplog. Patching logging.lastResort with None also '
               'does not help.')

    @pytest.mark.xfail(reason=_reason)
    def test_no_logging(self, make_logger, caplog):
        """Ensure that no log messages are emitted after handler removal."""
        with make_logger(logging.INFO) as logger:
            before_removal, after_removal = 'test before', 'test after'
            logger.addHandler(logging.StreamHandler())
            logger.warning(before_removal)
            assert before_removal in caplog.text
            close_all_handlers(logger)
            logger.warning(after_removal)
            assert after_removal not in caplog.text

    handler_methods= ('flush', 'close', 'acquire', 'release')

    def test_handler_methods_called(self, make_logger):
        """Check that the expected handler methods have been called."""
        handler = MockHandler()
        with make_logger() as logger:
            logger.addHandler(handler)
            close_all_handlers(logger)
        for method_name in self.handler_methods:
            method = getattr(handler, method_name)
            method.assert_called()

    @parametrize(method=handler_methods)
    @parametrize(exc=(OSError, ValueError))
    def test_handler_exc_caught(self, method, exc, make_logger):
        """Check that handler.method raising exc does not raise."""
        handler = MockHandler()
        with make_logger() as logger:
            logger.addHandler(handler)
            with make_obj_raise(handler, exc, method), not_raises(exc):
                close_all_handlers(logger)

    @parametrize(method=handler_methods)
    def test_handler_raises_propagated(self, method, make_logger):
        """Check that no funny exceptions are caught."""
        handler = MockHandler()
        with make_logger() as logger:
            logger.addHandler(handler)
            with raises_test_exception(handler, method):
                close_all_handlers(logger)


class TestDebugOrLower:
    """Tests for the debug_or_lower function."""

    _levels = {
        logging.DEBUG: True,
        logging.INFO: False,
        logging.CRITICAL: False,
        }

    @parametrize('level,result', _levels.items())
    @parametrize(effective=(True, False))
    def test_level_set(self, level, effective, result, make_logger):
        """Check outcome when the root logger is unset."""
        with make_logger(level) as logger:
            assert debug_or_lower(logger, effective) == result

    @parametrize('root_level,result', _levels.items())
    def test_root_level_set(self, root_level, result, make_logger):
        """Check outcome when the root logger is set."""
        with make_logger() as logger, at_level(logger.root, root_level):
            assert debug_or_lower(logger) == result

    def test_not_effective(self, make_logger):
        """Check outcome when the root logger is set."""
        with make_logger() as logger, at_level(logger.root, logging.ERROR):
            assert debug_or_lower(logger, effective=False)


class TestSilent:
    """Tests for logger_silent and logging_silent context managers."""

    outside_context = 'test outside'
    inside_context = 'test inside'

    def test_logger_logs_above_level(self, make_logger, caplog):
        """Check emission of logger messages above the silence level."""
        with make_logger(logging.INFO) as logger:
            logger.warning(self.outside_context)
            assert self.outside_context in caplog.text
            silenced = logger_silent(logger, level=logging.WARNING)
            with silenced():
                logger.error(self.inside_context)
            assert self.inside_context in caplog.text

    def test_logging_logs_above_level(self, caplog):
        """Check emission of messages above the module-silence level."""
        logging.warning(self.outside_context)
        assert self.outside_context in caplog.text
        with logging_silent(level=logging.WARNING):
            logging.error(self.inside_context)
        assert self.inside_context in caplog.text

    def test_no_logger_messages(self, make_logger, caplog):
        """Check that no messages are logged by a logger when silenced."""
        with make_logger() as logger:
            logger.warning(self.outside_context)
            assert self.outside_context in caplog.text
            silenced = logger_silent(logger)
            with silenced():
                logger.warning(self.inside_context)
            assert self.inside_context not in caplog.text

    def test_no_logging_messages(self, caplog):
        """Check that no messages are logged at module level when silenced."""
        logging.warning(self.outside_context)
        assert self.outside_context in caplog.text
        with logging_silent():
            logging.warning(self.inside_context)
        assert self.inside_context not in caplog.text


_with_console = {
    True: ('FileHandler', 'StreamHandler'),
    False: ('FileHandler',),
    }

@parametrize('console,handler_names', _with_console.items(), ids=_with_console)
def test_prepare_calc_logger(console, handler_names, make_logger, monkeypatch):
    """Check preparation of viperleed.calc logger."""
    mock_file = 'mock_log_file'
    with make_logger(logging.DEBUG) as logger, monkeypatch.context() as patch:
        for handler_name in handler_names:
            patch.setattr(logging, handler_name, MockHandler)
        prepare_calc_logger(logger, mock_file, with_console=console)
        assert len(logger.handlers) == len(handler_names)
        assert all(isinstance(h, MockHandler) for h in logger.handlers)
        assert logger.level == logging.INFO
