"""Module file of viperleed.calc.bookkeeper.history.

Defines the HistoryInfoFile class, a handler for the history.info file.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from enum import Enum
from pathlib import Path

from viperleed.calc.lib.log_utils import logging_silent
from viperleed.calc.bookkeeper.history.entry.entry import HistoryInfoEntry
from viperleed.calc.bookkeeper.history.entry.entry import PureCommentEntry

from ..log import LOGGER
from .constants import HISTORY_INFO_NAME
from .constants import HISTORY_INFO_SEPARATOR
from .errors import CantDiscardEntryError
from .errors import CantRemoveEntryError
from .errors import EntrySyntaxError
from .errors import NoHistoryEntryError


class CantRemoveReason(Enum):
    """Error messages explaining why an entry can't be removed."""

    FIXABLE = ('Contains fields with non-standard format (run '
               'bookkeeper --fixup to automatically fix this).')
    HAS_COMMENTS = 'Contains fields with user comments.'
    IS_COMMENT = 'Is a comment-only entry.'
    MISSES_FIELDS = 'Some expected fields were deleted.'
    NOT_UNDERSTOOD = 'Contains invalid fields that could not be interpreted.'

    def __str__(self):
        """Return the reason for removal as a message."""
        return self.value

    @classmethod
    def for_entry(cls, entry):
        """Return a reason suitable for a faulty `entry`."""
        if isinstance(entry, PureCommentEntry):
            return cls.IS_COMMENT
        if entry.has_comments:
            return cls.HAS_COMMENTS
        # This has to come before .was_understood as it's a subset
        if entry.misses_mandatory_fields:
            return cls.MISSES_FIELDS
        if not entry.was_understood:  # Can't be auto-fixed
            return cls.NOT_UNDERSTOOD
        return cls.FIXABLE


class HistoryInfoFile:
    """Deals with the history.info file in a history directory."""

    def __init__(self, file_path, create_new=False):
        """Initialize instance.

        Parameters
        ----------
        file_path : str or Path
            The path to the history.info file.
        create_new : bool, optional
            Whether a empty file should be created
            in case `file_path` does not exist.

        Raises
        ------
        FileNotFoundError
            If `file_path` does not exist and `create_new` is False.
        """
        self.path = Path(file_path)
        if not self.path.is_file() and not create_new:
            raise FileNotFoundError(f'{HISTORY_INFO_NAME} file not '
                                    f'found at {self.path}.')
        if not self.path.is_file():  # create new file
            self.path.touch()

        self.raw_contents = ''
        self._entries = []  # The HistoryInfoEntry(s) from raw_contents
        self._time_format = None  # To keep entries consistent

    @property
    def last_entry(self):
        """Return the most recent entry that was read from the current file."""
        try:
            return self._entries[-1]
        except IndexError:
            return None

    @property
    def last_entry_was_discarded(self):
        """Return whether the last entry is labeled as DISCARDED."""
        return getattr(self.last_entry, 'is_discarded', False)

    def append_entry(self, entry, fix_time_format=True):
        """Write `entry` to the history.info file."""
        skip_separator = (
            self.raw_contents.strip().endswith(HISTORY_INFO_SEPARATOR.strip())
            or not self.raw_contents.strip()
            )
        if fix_time_format:
            entry = entry.with_time_format(self._time_format or 'default')
            self._time_format = entry.time_format
        self._entries.append(entry)
        with open(self.path, 'a', encoding='utf-8') as history_info:
            if not skip_separator:
                history_info.write(HISTORY_INFO_SEPARATOR)
                self.raw_contents += HISTORY_INFO_SEPARATOR
            entry_text = str(entry)
            history_info.write(entry_text)
            self.raw_contents += entry_text

    def discard_last_entry(self):
        """Mark the last entry in the history.info file as discarded."""
        self.may_discard_last_entry()
        last_entry = self.last_entry
        discarded_entry = last_entry.as_discarded()
        if discarded_entry is last_entry:  # Already discarded
            return
        self._do_remove_last_entry()
        self.append_entry(discarded_entry, fix_time_format=False)

    def may_discard_last_entry(self):
        """Raise if the last entry cannot be discarded."""
        last_entry = self.last_entry
        if last_entry is None:
            raise NoHistoryEntryError('No entries to discard.')
        if isinstance(last_entry, PureCommentEntry):
            raise CantDiscardEntryError('It contains only comments.')
        if last_entry.is_discarded:
            # This is not really a case for an exception. The entry
            # is discarded already so we simply don't have to bother
            LOGGER.warning('Last entry is already discarded.')

    def may_remove_last_entry(self):
        """Raise if the last entry can't be removed."""
        last_entry = self.last_entry
        if not last_entry:
            raise NoHistoryEntryError('No entries to remove.')

        if last_entry.can_be_removed:
            return

        # Comment-only entries technically do not have notes
        is_comment = isinstance(last_entry, PureCommentEntry)

        # Check for notes in history.info
        if not is_comment and last_entry.has_notes:
            raise CantRemoveEntryError(
                f'The last entry in {self.path.name} has user '
                'notes. If you really want to purge the last run, '
                'remove the notes first.'
                )

        # Allow removing entries that have out-of-date format
        # but are consistent with the rest of the file
        is_consistent = (
            not is_comment  # PureCommentEntry has no .is_only_outdated
            and last_entry.is_only_outdated
            and last_entry.time_format is self._time_format
            )
        if is_consistent:
            return

        err_ = CantRemoveReason.for_entry(last_entry)
        raise CantRemoveEntryError('Failed to remove last entry from '
                                   f'{self.path.name}: {err_} '
                                   'Please proceed manually.')

    def read(self):
        """Read the current contents of the history.info file."""
        self.raw_contents = self.path.read_text(encoding='utf-8')
        self._read_all_entries()
        self._infer_time_format()

    def remove_last_entry(self):
        """Remove the last entry from the history.info file."""
        self.may_remove_last_entry()
        self._do_remove_last_entry()

    def _do_remove_last_entry(self):
        """Actually remove the last entry assuming its possible."""
        self._entries.pop()
        if self._entries:
            content_without_last, *_ = self.raw_contents.rsplit(
                HISTORY_INFO_SEPARATOR.strip(),
                maxsplit=1
                )
        else:  # Nothing left
            content_without_last = ''
        if content_without_last.endswith('\n'):
            content_without_last = content_without_last[:-len('\n')]
        # Clear file and write back entries
        self.path.write_text(content_without_last, encoding='utf-8')
        self.raw_contents = content_without_last

    def _infer_time_format(self):
        """Parse all entries to find out the format of TIME fields."""
        if self._time_format:  # Only once
            return
        for entry in self._entries:
            try:
                self._time_format = entry.time_format
            except AttributeError:    # Pure comment
                continue
            return

    def _read_all_entries(self):                                                # TODO: should log only once per entry?
        """Parse the raw_contents for all entries and store them."""
        self._entries.clear()
        if not self.raw_contents:
            return
        sep = HISTORY_INFO_SEPARATOR.strip()
        for entry_str in self.raw_contents.split(sep):
            self._entries.append(HistoryInfoEntry.from_string(entry_str))
