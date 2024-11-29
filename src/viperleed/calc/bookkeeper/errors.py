"""Module errors of viperleed.calc.bookkeeper.

Defines the base BookkeperError exception and its subclasses.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-11-29'
__license__ = 'GPLv3+'


class BookkeperError(Exception):
    """Base class of all bookkeeper-related errors."""


class _FileNotOlderError(BookkeperError):
    """Exception used internally for file-age checks."""
