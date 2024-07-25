"""Module errors of viperleed.calc.bookkeeper.history_info.

Defines exceptions related to handling the history.info file.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-25'
__license__ = 'GPLv3+'


class HistoryInfoError(Exception):
    """Base class for all errors related to the history.info file."""


class CantRemoveEntryError(HistoryInfoError):
    """An entry cannot be removed."""


class NoHistoryEntryError(HistoryInfoError):
    """There is no entry to process according to the criteria."""


class EntrySyntaxError(HistoryInfoError):
    """Something is wrong with the syntax of an entry."""


class FixableSyntaxError(HistoryInfoError):
    """An entry has some syntax error, but it can be fixed."""


class _PureCommentEntryError(HistoryInfoError):
    """Exception used internally to decide that an entry is comment only."""
