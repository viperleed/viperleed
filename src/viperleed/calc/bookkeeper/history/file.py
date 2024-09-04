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

from pathlib import Path

from viperleed.calc.lib.log_utils import logging_silent
from viperleed.calc.bookkeeper.history.entry.entry import HistoryInfoEntry
from viperleed.calc.bookkeeper.history.entry.entry import PureCommentEntry

from ..log import LOGGER
from .constants import HISTORY_INFO_NAME
from .constants import HISTORY_INFO_SEPARATOR
from .errors import CantRemoveEntryError
from .errors import EntrySyntaxError
from .errors import NoHistoryEntryError


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
        self.last_entry = None
        self._time_format = None  # To keep consistency
        self.read()

    @property
    def last_entry_was_discarded(self):
        """Return whether the last entry is labeled as DISCARDED."""
        return getattr(self.last_entry, 'discarded', False)

    def append_entry(self, entry):
        """Write `entry` to the history.info file."""
        skip_separator = (
            self.raw_contents.strip().endswith(HISTORY_INFO_SEPARATOR.strip())
            or not self.raw_contents.strip()
            )
        with open(self.path, 'a', encoding='utf-8') as history_info:
            if not skip_separator:
                history_info.write(HISTORY_INFO_SEPARATOR)
                self.raw_contents += HISTORY_INFO_SEPARATOR
            entry_text = str(entry)
            history_info.write(entry_text)
            self.raw_contents += entry_text
            self.last_entry = entry

    def discard_last_entry(self):
        """Mark the last entry in the history.info file as discarded."""
        if self.last_entry is None:
            raise NoHistoryEntryError('No entries to discard.')
        if isinstance(self.last_entry, PureCommentEntry):
            LOGGER.warning('Cannot discard last entry. '
                           'It contains only comments.')
            return
        if self.last_entry.discarded:
            LOGGER.warning('Last entry is already discarded.')
            return
        discarded_entry = self.last_entry.as_discarded()
        self._do_remove_last_entry()
        self.append_entry(discarded_entry)

    def may_remove_last_entry(self):
        """Raise if the last entry can't be removed."""
        if not self.last_entry:
            raise NoHistoryEntryError('No entries to remove.')

        if self.last_entry.can_be_removed:
            return

        # Comment-only entries technically do not have notes
        is_comment = isinstance(self.last_entry, PureCommentEntry)

        # Check for notes in history.info
        if not is_comment and self.last_entry.has_notes:
            raise CantRemoveEntryError(
                f'The last entry in {self.path.name} has user '
                'notes. If you really want to purge the last run, '
                'remove the notes first.'
                )

        if is_comment:
            err_ = 'is a comment-only field'
        elif self.last_entry.misses_mandatory_fields:
            # This has to come before .was_understood as it's a subset
            err_ = 'some expected fields were deleted'
        elif not self.last_entry.was_understood:  # Can't be auto-fixed
            err_ = 'contains invalid fields that could not be interpreted'
        else:    # Needs fixing, but can be done in --fixup mode                # TODO: implement
            assert not self.last_entry.has_notes  # Checked above
            err_ = ('contains fields with non-standard format (run '
                    'bookkeeper --fixup to automatically fix it)')
        raise CantRemoveEntryError('Failed to remove last entry from '
                                   f'{self.path.name}: {err_}. Please '
                                   'proceed manually.')

    def read(self):
        """Read the current contents of the history.info file."""
        self.raw_contents = self.path.read_text(encoding='utf-8')
        if not self.raw_contents:
            self.last_entry = None
            return
        self._infer_time_format()
        *_, entry_str = self.raw_contents.split(HISTORY_INFO_SEPARATOR.strip())
        entry = HistoryInfoEntry.from_string(entry_str)                         # TODO: how to test this properly? It does not raise. Perhaps a TolerantHistoryInfoEntry?
        if self._time_format:
            try:
                entry = entry.with_time_format(self._time_format)
            except AttributeError:  # Pure comment
                pass
        self.last_entry = entry

    def remove_last_entry(self):
        """Remove the last entry from the history.info file."""
        self.may_remove_last_entry()
        self._do_remove_last_entry()

    def _do_remove_last_entry(self):
        """Actually remove the last entry assuming its possible."""
        if HISTORY_INFO_SEPARATOR.strip() in self.raw_contents:
            content_without_last, *_ = self.raw_contents.rsplit(
                HISTORY_INFO_SEPARATOR.strip(),
                maxsplit=1
                )
        else:  # Only one entry
            content_without_last = ''
        if content_without_last.endswith('\n'):
            content_without_last = content_without_last[:-len('\n')]
        # Clear file and write back entries
        self.path.write_text(content_without_last, encoding='utf-8')
        self.read()  # Re-read to update last_entry

    def _infer_time_format(self):
        """Parse an entry to find out the format of TIME fields."""
        if self._time_format:  # Only once
            return
        sep = HISTORY_INFO_SEPARATOR.strip()
        with logging_silent():
            for entry_str in self.raw_contents.split(sep):
                try:
                    entry = HistoryInfoEntry.from_string(entry_str)
                except EntrySyntaxError:  # Not a valid entry
                    continue
                try:
                    self._time_format = entry.time_format
                except AttributeError:    # Pure comment
                    continue
                return
