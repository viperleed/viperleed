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


import logging
from logging import DEBUG


class CustomLogFormatter(logging.Formatter):
    """Logging Formatter for level-dependent message formatting."""

    formats = {
        logging.DEBUG: 'dbg: %(msg)s',
        logging.INFO: '%(msg)s',
        logging.WARNING: '# WARNING: %(msg)s',
        logging.ERROR: (
            '### ERROR ### in %(module)s:%(funcName)s:%(lineno)s\n'
            '# %(msg)s \n#############'
            ),
        logging.CRITICAL: (
            '### CRITICAL ### in %(module)s:%(funcName)s:%(lineno)s\n'
            '# %(msg)s \n################'
            ),
        'DEFAULT': '%(msg)s',
        }

    def format(self, record):
        """Debug log format for everything at DEBUG level or lower."""
        level = record.levelno
        log_fmt = (self.formats[DEBUG] if level < DEBUG
                   else self.formats.get(level, self.formats['DEFAULT']))
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
