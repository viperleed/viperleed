"""Module constants of viperleed.calc.bookkeeper.

Collects the definitions of a few constants used
in multiple modules of the bookkeeper package.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

import logging

from viperleed.calc.lib.base import parent_name


BOOKIE_LOGFILE = 'bookkeeper.log'  # Persistent among runs
LOGGER = logging.getLogger(parent_name(__name__))


def _get_handlers(handler_type, **with_attrs):
    """Yield handlers attached to LOGGER that are of a given type."""

    def _get_handlers_for(a_logger):
        for handler in a_logger.handlers:
            if not isinstance(handler, handler_type):
                continue
            if not with_attrs:
                yield handler
                continue
            try:
                attrs = {attr: getattr(handler, attr) for attr in with_attrs}
            except AttributeError:
                continue
            attrs_match = (
                value == expected
                for value, expected in zip(attrs, with_attrs.values())
                )
            if all(attrs_match):
                yield handler

    # This implementation mimics the one of logging.Logger.hasHandlers
    logger = LOGGER  # Start from this one but walk up the hierarchy
    while logger:
        yield from _get_handlers_for(logger)
        if not logger.propagate:
            break
        logger = logger.parent


def add_bookeeper_logfile(at_path):
    """Attach a FileHandler to LOGGER for the bookkeeper log file."""
    bookkeeper_log = at_path / BOOKIE_LOGFILE
    has_log = _get_handlers(logging.FileHandler,
                            baseFilename=str(bookkeeper_log))
    if not any(has_log):
        LOGGER.addHandler(logging.FileHandler(bookkeeper_log, mode='a'))


def ensure_has_stream_handler():
    """Attach a stream handler to LOGGER unless one exists already."""
    if not any(_get_handlers(logging.StreamHandler)):
        LOGGER.addHandler(logging.StreamHandler())
    LOGGER.setLevel(logging.INFO)
    LOGGER.propagate = True
