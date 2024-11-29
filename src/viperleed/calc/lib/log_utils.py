"""Module log_utils of viperleed.calc.lib.

Collects functions and classes useful for handling logging features.
Part of the functionality in this module comes from calc.lib.base,
originally created on 2019-06-13.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-05'
__license__ = 'GPLv3+'


from contextlib import contextmanager
import logging
from logging import DEBUG


@contextmanager
def at_level(logger, level):
    """Temporarily set the level of `logger` to `level`."""
    # Fetch the current level, but don't use getEffectiveLevel
    # that walks up the parent tree to fetch up to the root
    previous_level = logger.level
    logger.setLevel(level)
    try:
        yield
    finally:
        logger.setLevel(previous_level)


def close_all_handlers(logger):
    """Close all handlers of `logger`."""
    while logger.handlers:
        handler = logger.handlers[-1]
        logger.removeHandler(handler)
        # The next bit looks similar to logging.shutdown, but only
        # removes and closes the handlers for a specific logger
        try:
            _flush_and_close_handler(handler)
        except (OSError, ValueError):
            # Handler already closed, but still lingering around.
            pass


def debug_or_lower(logger, effective=True):
    """Return whether `logger` has a level less than DEBUG."""
    level = logger.getEffectiveLevel() if effective else logger.level
    return level <= DEBUG


@contextmanager
def logging_silent(level=logging.CRITICAL):
    """Mute all logging messages below `level` for ALL loggers.

    Parameters
    ----------
    level : int, optional
        The maximum logging level in use. Default is
        logging.CRITICAL, which silences all logging
        messages.

    Returns
    -------
    None.
    """
    # Adapted from https://gist.github.com/simon-weber/7853144
    # Uses an undocumented, but public API: logging.Manager.disable
    previous_level = logging.root.manager.disable
    logging.disable(level)
    try:
        yield
    finally:
        logging.disable(previous_level)


def logger_silent(logger, level=logging.CRITICAL):
    """Return a context that silences `logger` messages below `level`."""
    @contextmanager
    def _context():
        with at_level(logger, level):
            yield
    return _context


def prepare_calc_logger(logger, file_name, with_console):
    """Prepare logger to be used with calc.run."""
    logger.setLevel(logging.INFO)
    formatter = CalcLogFormatter()
    file_handler = logging.FileHandler(file_name, mode='w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if not with_console:
        return
    stderr_handler = logging.StreamHandler()  # Uses sys.stderr
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)


class CalcLogFormatter(logging.Formatter):
    """Logging Formatter for level-dependent message formatting."""

    formats = {
        logging.DEBUG: 'dbg: %(msg)s',
        logging.INFO: '%(msg)s',
        logging.WARNING: '# WARNING: %(msg)s',
        logging.ERROR: (
            '### ERROR ### in %(module)s:%(funcName)s:%(lineno)s\n'
            '# %(msg)s\n#############'
            ),
        logging.CRITICAL: (
            '### CRITICAL ### in %(module)s:%(funcName)s:%(lineno)s\n'
            '# %(msg)s\n################'
            ),
        'DEFAULT': '%(msg)s',
        }

    def format(self, record):
        """Use the DEBUG log format for everything at DEBUG level or lower."""
        level = record.levelno
        log_fmt = (self.formats[DEBUG] if level < DEBUG
                   else self.formats.get(level, self.formats['DEFAULT']))
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


@contextmanager
def _acquire_handler(handler):
    """Safely acquire a handler, then release it."""
    try:  # pylint: disable=too-many-try-statements
        handler.acquire()
        yield
    finally:
        handler.release()


def _flush_and_close_handler(handler):
    """Flush all the non-handled record, then close a handler."""
    with _acquire_handler(handler):
        handler.flush()
        handler.close()
