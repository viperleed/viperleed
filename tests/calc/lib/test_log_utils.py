"""Tests for the log_utils module of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-08-05'
__license__ = 'GPLv3+'

from contextlib import contextmanager
import logging
import logging.handlers
from logging import LogRecord as Record
from unittest.mock import Mock

import pytest
from pytest_cases import fixture
from pytest_cases import fixture_ref
from pytest_cases import parametrize

from viperleed.calc.lib.log_utils import at_level
from viperleed.calc.lib.log_utils import close_all_handlers
from viperleed.calc.lib.log_utils import debug_or_lower
from viperleed.calc.lib.log_utils import get_handlers
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
    def _make(*handlers, level=None, name='test_logger'):
        logger = logging.getLogger(name)
        if level is not None:
            logger.setLevel(level)
        for handler in handlers:
            logger.addHandler(handler)
        try:
            yield logger
        finally:
            # logging uses the garbage collector to
            # clean up loggers. Let's do the same here.
            logger.manager.loggerDict.pop(name, None)
    return _make


class CustomHandler(logging.Handler):
    """A mock subclass for a logging.Handler."""

    def __init__(self, level=logging.NOTSET, custom_attr=None):
        """Initialize handler at a level and with a custom attribute."""
        super().__init__(level)
        self.custom_attr = custom_attr

    emit = Mock()


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
        with make_logger(level=logging.WARNING) as logger:
            assert logger.level == logging.WARNING
            with at_level(logger, logging.DEBUG):
                assert logger.level == logging.DEBUG
            assert logger.level == logging.WARNING

    def test_messages_silenced(self, make_logger, caplog):
        """Check that at_level silences messages."""
        with make_logger(level=logging.DEBUG) as logger:
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
        handlers = logging.StreamHandler(), logging.StreamHandler()
        with make_logger(*handlers) as logger:
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
        handler = logging.StreamHandler()
        with make_logger(handler, level=logging.INFO) as logger:
            before_removal, after_removal = 'test before', 'test after'
            logger.warning(before_removal)
            assert before_removal in caplog.text
            close_all_handlers(logger)
            logger.warning(after_removal)
            assert after_removal not in caplog.text

    handler_methods= ('flush', 'close', 'acquire', 'release')

    def test_handler_methods_called(self, make_logger):
        """Check that the expected handler methods have been called."""
        handler = MockHandler()
        with make_logger(handler) as logger:
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
        with make_logger(level=level) as logger:
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


class TestGetHandlers:
    """Tests for the get_handlers function."""

    _fmt = '%(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(_fmt)
    handlers_and_levels = (
        (logging.StreamHandler(), logging.WARNING),
        (logging.handlers.MemoryHandler(1), logging.ERROR),
        (CustomHandler(custom_attr='value'), logging.INFO),
        )
    handlers = tuple(h for h, _ in handlers_and_levels)
    log_filter = logging.Filter(name='test_filter')

    @staticmethod
    def _check_handler_attrs(handler, expected):
        """Check that handler has all expected attributes."""
        assert all(getattr(handler, attr) == value
                   for attr, value in expected.items())

    def _prepare_handlers(self):
        """Set up self.handlers."""
        for handler, level in self.handlers_and_levels:
            handler.setFormatter(self.formatter)
            handler.addFilter(self.log_filter)
            handler.setLevel(level)

    @fixture(name='parented_logger')
    def fixture_parented_logger(self, make_logger):
        """Yield a logger with a bunch of handlers and a parent."""
        self._prepare_handlers()
        parent_handler, *child_handlers = self.handlers
        with make_logger(parent_handler, name='parent') as parent_logger, \
                make_logger(*child_handlers, name='parent.child') as logger:
            parent_logger.parent = None  # Avoid pytest handlers
            return logger, (*child_handlers, parent_handler)

    @fixture(name='parentless_logger')
    def fixture_parentless_logger(self, make_logger):
        """Yield a logger with a bunch of handlers."""
        self._prepare_handlers()
        with make_logger(*self.handlers) as logger:
            logger.parent = None
            return logger, self.handlers

    _logger_fixtures = (
        fixture_ref('parented_logger'),
        fixture_ref('parentless_logger'),
        )
    _loggers = (
        *_logger_fixtures,
        (logging.Logger('no_handlers_logger'), ()),
        )

    @parametrize('logger,handlers', _loggers)
    def test_get_all_types(self, logger, handlers):
        """Check get_handlers of any type."""
        result = tuple(get_handlers(logger, logging.Handler))
        assert result == handlers

    _handler_types = (
        logging.StreamHandler,
        logging.handlers.MemoryHandler,
        CustomHandler,
        CustomHandler(),
        )

    @parametrize(handler_type=_handler_types)
    @parametrize(logger_fixture=_logger_fixtures)
    def test_get_one_type(self, handler_type, logger_fixture):
        """Check correct result of fetching handlers of one type."""
        logger, _ = logger_fixture
        handlers = tuple(get_handlers(logger, handler_type))
        if isinstance(handler_type, logging.Handler):
            handler_type = type(handler_type)

        # Only one handler per exact type at setup
        # pylint: disable-next=unidiomatic-typecheck
        result = next(h for h in handlers if type(h) is handler_type)
        assert isinstance(result, handler_type)

        # The rest may be a subclass, but not an exact type
        rest = (h for h in handlers if h is not result)
        # pylint: disable-next=unidiomatic-typecheck
        assert all(isinstance(h, handler_type) and type(h) != handler_type
                   for h in rest)

    _with_attrs = {  # handler_type, attrs, n_expect
        'custom value OK': (CustomHandler, {'custom_attr': 'value'}, 1),
        'by level': (logging.Handler, {'level': logging.ERROR}, 1),
        'by formatter': (logging.Handler, {'formatter': formatter}, 3),
        'by filter': (logging.Handler, {'filters': [log_filter]}, 3),
        'wrong value': (CustomHandler, {'custom_attr': 'wrong'}, 0),
        'attr does not exist': (CustomHandler, {'wrong_attr': 'value'}, 0),
        'attr of other': (CustomHandler, {'level': logging.ERROR}, 0),
        'other handler': (logging.FileHandler, {}, 0),
        'not a handler': (MockHandler, {}, 0),
        }

    @parametrize('handler_type,attrs,n_expect',
                 _with_attrs.values(),
                 ids=_with_attrs)
    @parametrize(logger_fixture=_logger_fixtures)
    def test_get_with_attrs(self, handler_type, attrs, n_expect,
                            logger_fixture):
        """Check expected outcome when filtering attributes."""
        logger, _ = logger_fixture
        result = tuple(get_handlers(logger, handler_type, **attrs))
        assert len(result) == n_expect
        for handler in result:
            self._check_handler_attrs(handler, attrs)

    def test_parented_not_propagating(self, parented_logger):
        """Check correct result for a parented, but not propagating logger."""
        logger, (*handlers, parent_handler) = parented_logger
        logger.propagate = False
        result = tuple(get_handlers(logger, logging.Handler))
        assert result == tuple(handlers)
        assert parent_handler not in result


class TestSilent:
    """Tests for logger_silent and logging_silent context managers."""

    outside_context = 'test outside'
    inside_context = 'test inside'

    def test_logger_logs_above_level(self, make_logger, caplog):
        """Check emission of logger messages above the silence level."""
        with make_logger(level=logging.INFO) as logger:
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
    with make_logger(level=logging.DEBUG) as logger, \
            monkeypatch.context() as patch:
        for handler_name in handler_names:
            patch.setattr(logging, handler_name, MockHandler)
        prepare_calc_logger(logger, mock_file, with_console=console)
        assert len(logger.handlers) == len(handler_names)
        assert all(isinstance(h, MockHandler) for h in logger.handlers)
        assert logger.level == logging.INFO
