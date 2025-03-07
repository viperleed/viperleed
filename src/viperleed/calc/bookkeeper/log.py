"""Module constants of viperleed.calc.bookkeeper.

Collects the definitions of a few constants used
in multiple modules of the bookkeeper package.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

import logging

from viperleed.calc.lib.log_utils import get_handlers
from viperleed.calc.lib.string_utils import parent_name


BOOKIE_LOGFILE = 'bookkeeper.log'  # Persistent among runs
LOGGER = logging.getLogger(parent_name(__name__))


def add_bookkeeper_logfile(at_path):
    """Attach a FileHandler to LOGGER for the bookkeeper log file."""
    bookkeeper_log = at_path / BOOKIE_LOGFILE
    has_log = get_handlers(LOGGER,
                           logging.FileHandler,
                           baseFilename=str(bookkeeper_log))
    if not any(has_log):
        LOGGER.addHandler(logging.FileHandler(bookkeeper_log, mode='a'))


def ensure_has_stream_handler():
    """Attach a stream handler to LOGGER unless one exists already."""
    if not any(get_handlers(LOGGER, logging.StreamHandler)):
        LOGGER.addHandler(logging.StreamHandler())
    LOGGER.setLevel(logging.INFO)
    LOGGER.propagate = True
