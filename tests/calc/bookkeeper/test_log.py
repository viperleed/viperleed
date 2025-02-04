"""Tests for module log of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-09-29'
__license__ = 'GPLv3+'

import logging

from pytest_cases import fixture

from viperleed.calc.bookkeeper.log import add_bookkeeper_logfile
from viperleed.calc.bookkeeper.log import ensure_has_stream_handler


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


@fixture(name='logger')
def fixture_logger(monkeypatch):
    """Replace the bookkeeper logger with a fake one."""
    with monkeypatch.context() as patch:
        logger = MockLogger()
        patch.setattr('viperleed.calc.bookkeeper.log.LOGGER', logger)
        yield logger
        for handler in logger.handlers:
            handler.close()


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
