"""Module info of viperleed.calc.bookkeeper.history.

Defines the HistoryInfoFile class, a handler for the history.info file.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from enum import Enum
from functools import wraps
import logging
from pathlib import Path
import re

from viperleed.calc.bookkeeper.history.entry.entry import HistoryInfoEntry
from viperleed.calc.bookkeeper.history.entry.entry import PureCommentEntry

from ..mode import BookkeeperMode as Mode
from .constants import HISTORY_INFO_NAME
from .constants import HISTORY_INFO_SEPARATOR
from .errors import CantDiscardEntryError
from .errors import CantRemoveEntryError
from .errors import FixFailedError
from .errors import NoHistoryEntryError


LOGGER = logging.getLogger(__name__)


class CantRemoveReason(Enum):
    """Error messages explaining why an entry can't be removed."""

    FIXABLE = ('Contains fields with non-standard format (run \'bookkeeper '
               f'{Mode.FIX.long_flag}\' to automatically fix this).')
    HAS_COMMENTS = 'Contains fields with user comments.'
    IS_COMMENT = 'Is a comment-only entry.'
    IS_EMPTY = 'Has no contents.'
    MISSES_FIELDS = 'Some expected fields were deleted.'
    NOT_UNDERSTOOD = 'Contains invalid fields that could not be interpreted.'

    def __str__(self):
        """Return the reason for removal as a message."""
        return self.value

    @classmethod
    def for_entry(cls, entry):
        """Return a reason suitable for a faulty `entry`."""
        is_comment = isinstance(entry, PureCommentEntry)
        if is_comment and not entry.raw_comment.strip():
            return cls.IS_EMPTY
        if is_comment:
            return cls.IS_COMMENT
        if entry.has_comments:
            return cls.HAS_COMMENTS
        # This has to come before .was_understood as it's a subset
        if entry.misses_mandatory_fields:
            return cls.MISSES_FIELDS
        if not entry.was_understood:  # Can't be auto-fixed
            return cls.NOT_UNDERSTOOD
        return cls.FIXABLE


def _create_or_raise(func):
    """Create a history.info file if asked or raise FileNotFoundError."""
    @wraps(func)
    def _wrapper(self, *args, **kwargs):
        # pylint: disable-next=protected-access
        if not self.path.is_file() and not self._create_new:
            raise FileNotFoundError(f'{self.path.name} file not '
                                    f'found at {self.path.parent}.')
        if not self.path.is_file():  # create new file
            self.path.touch()
        return func(self, *args, **kwargs)
    return _wrapper


class HistoryInfoFile:
    """Deals with the history.info file in a history directory."""

    def __init__(self, root_path, create_new=False):
        """Initialize instance.

        Parameters
        ----------
        root_path : str or Path
            The path to the folder containing the history.info file.
        create_new : bool, optional
            Whether a empty file should be created in case no
            history.info file exists at `root_path`. Default
            is False.

        Returns
        -------
        None.
        """
        self.path = Path(root_path).resolve() / HISTORY_INFO_NAME
        self._create_new = create_new

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

    @_create_or_raise
    def append_entry(self, entry, fix_time_format=True):
        """Write `entry` to the history.info file."""
        if fix_time_format:
            entry = entry.with_time_format(self._time_format or 'default')
            # NB: there's no need to differentiate whether we already
            # had a _time_format here: if we didn't, this is the first
            # entry, if we did, entry.time_format is self._time_format.
            self._time_format = entry.time_format
        new_text = self._new_separator
        entry_text = str(entry)
        if not self._entries:
            # Remove the leading newline for the first entry
            entry_text = entry_text[1:]
        new_text += entry_text
        self._entries.append(entry)
        self.raw_contents += new_text
        try:  # pylint: disable=too-many-try-statements  # OK with open
            with self.path.open('a', encoding='utf-8') as history_info:
                history_info.write(new_text)
        except OSError:
            LOGGER.error('Failed to append entry to '
                         f'{HISTORY_INFO_NAME} file:\n{new_text}')
            raise

    def discard_last_entry(self):
        """Mark the last entry in the history.info file as discarded."""
        self.may_discard_last_entry()
        last_entry = self.last_entry
        discarded_entry = last_entry.as_discarded()
        if discarded_entry is last_entry:  # Already discarded
            return
        self._do_remove_last_entry()
        self.append_entry(discarded_entry, fix_time_format=False)

    @_create_or_raise
    def fix(self):
        """Save to file fixed versions of all the fixable entries."""
        fixed_entries = []
        for entry in self._entries:
            try:
                fixed = entry.as_fixed()
            except AttributeError:
                assert isinstance(entry, PureCommentEntry)
                fixed = entry
            except FixFailedError as exc:
                LOGGER.error(str(exc))
                fixed = entry
            fixed_entries.append(fixed)
        self._entries = fixed_entries
        old_raw = self.raw_contents
        new_raw = HISTORY_INFO_SEPARATOR.join(str(e) for e in self._entries)
        if not old_raw.startswith('\n') and new_raw.startswith('\n'):
            new_raw = new_raw[1:]
        self.raw_contents = new_raw
        self._time_format = None  # So that we can _infer_time_format
        self._infer_time_format()
        self.path.write_text(self.raw_contents, encoding='utf-8')

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
        if not last_entry or not self.raw_contents.strip():
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

    @_create_or_raise
    def read(self):
        """Read the current contents of the history.info file."""
        self.raw_contents = self.path.read_text(encoding='utf-8')
        self._read_all_entries()
        self._infer_time_format()

    def remove_last_entry(self):
        """Remove the last entry from the history.info file."""
        self.may_remove_last_entry()
        self._do_remove_last_entry()

    @property
    def _new_separator(self):
        """Return a separator useful for appending a new entry."""
        if not self.raw_contents.strip():
            # First entry
            return ''
        if self.raw_contents.endswith(HISTORY_INFO_SEPARATOR):
            # There's already a separator
            return ''
        missing = HISTORY_INFO_SEPARATOR
        if self.raw_contents.endswith(HISTORY_INFO_SEPARATOR.rstrip()):
            # There's part of a separator. Return what's missing
            return re.sub(rf'^{missing.rstrip()}', '', missing, count=1)
        return missing

    @_create_or_raise
    def _do_remove_last_entry(self):
        """Actually remove the last entry assuming it is possible."""
        self._entries.pop()
        if self._entries:
            content_without_last, *_ = self.raw_contents.rsplit(
                HISTORY_INFO_SEPARATOR,
                maxsplit=1,
                )
        else:  # Nothing left
            content_without_last = ''
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
        entries_str = self.raw_contents.split(HISTORY_INFO_SEPARATOR)
        for i, entry_str in enumerate(entries_str, 1):
            try:
                self._entries.append(HistoryInfoEntry.from_string(entry_str))
            except ValueError:  # Some bits of a separator
                if i < len(entries_str):
                    # Only acceptable at the last one
                    raise
                assert entry_str.endswith(HISTORY_INFO_SEPARATOR.rstrip())
                # The last entry has some bits of a separator at the end
                entry_str = re.sub(rf'{HISTORY_INFO_SEPARATOR.rstrip()}$', '',
                                   entry_str, count=1)
                self._entries.append(HistoryInfoEntry.from_string(entry_str))
