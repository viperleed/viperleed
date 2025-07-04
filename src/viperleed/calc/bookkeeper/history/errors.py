"""Module errors of viperleed.calc.bookkeeper.history.

Defines exceptions related to handling the history.info file.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-07-25'
__license__ = 'GPLv3+'

from enum import Enum

from viperleed.calc.bookkeeper.errors import BookkeeperError


class HistoryInfoError(BookkeeperError):
    """Base class for all errors related to the history.info file."""


class CantDiscardEntryError(HistoryInfoError):
    """An entry cannot be marked as DISCARDED."""


class CantRemoveEntryError(HistoryInfoError):
    """An entry cannot be removed."""


class EntrySyntaxError(HistoryInfoError):
    """Something is wrong with the syntax of an entry."""

    def __init__(self, reason):
        """Initialize instance."""
        if isinstance(reason, Enum):
            reason = reason.value
        super().__init__(reason)


class FieldsScrambledError(HistoryInfoError):
    """Some lines in a history.info entry are in the wrong order."""


class FixableSyntaxError(HistoryInfoError):
    """An entry has some syntax error, but it can be fixed."""

    def __init__(self, reason='', fixed_value=None, action=None):
        """Initialize instance.

        Parameters
        ----------
        reason : str, optional
            The reason for this syntax fix. This is the exception
            message. Available as the .reason attribute. Default
            is an empty string.
        fixed_value : object, optional
            The value that fixes the problem. Available as the
            .fixed_value attribute. Default is None.
        action : FixAction, optional
            Which action should be performed to fix this error.
            Available as the .action attribute. Default is None.

        Returns
        -------
        None.
        """
        self.action = action
        self.fixed_value = fixed_value
        self.reason = reason
        super().__init__(reason)


class FixFailedError(HistoryInfoError):
    """Failed to fix an entry or its fields."""


class NoHistoryEntryError(HistoryInfoError):
    """There is no entry to process according to the criteria."""


class MetadataError(Exception):
    """Base exception for metadata-related errors."""


class MetadataMismatchError(MetadataError):
    """Something is not consistent between the metadata-file and its folder."""
