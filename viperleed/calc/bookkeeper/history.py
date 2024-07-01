"""Module history of viperleed.calc.bookkeeper."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from .constants import HISTORY_INFO_NAME
from .constants import LOGGER


_HISTORY_INFO_SPACING = 12  # For the leftmost field in history.info
HISTORY_INFO_SEPARATOR = '\n###########\n'


class HistoryInfoError(Exception):
    """Base class for all errors related to the history.info file."""


class NoHistoryEntryError(HistoryInfoError):
    """There is no entry to process according to the criteria."""


class HistoryInfoFile:
    """Deals with the history.info file in a history directory."""

    def __init__(self, file_path, create_new=False):
        """Initialize instance.

        Parameters
        ----------
        file_path : str or Path
            The path to the history.info file.
        create_new : bool, optional
            Whether a empty file should be created in case file_path
            does not exist.

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
        self.read()

    def read(self):
        """Read the current contents of the history.info file."""
        with open(self.path, 'r', encoding='utf-8') as history_info:
            self.raw_contents = history_info.read()
        if not self.raw_contents:
            self.last_entry = None
        else:
            self.last_entry = self._parse_entry(
                self.raw_contents.split(HISTORY_INFO_SEPARATOR.strip())[-1]
                )

    @property
    def last_entry_has_notes(self):
        """Return whether the last entry has associated notes."""
        if self.last_entry is None:
            return False
        return bool(self.last_entry.notes)

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
        last_entry = self.last_entry
        if last_entry.discarded:
            LOGGER.warning('Last entry is already discarded.')
        last_entry.discarded = True
        self.remove_last_entry()
        self.append_entry(last_entry)

    def remove_last_entry(self):
        """Discard the last entry from the history.info file."""
        if self.last_entry is None:
            raise NoHistoryEntryError('No entries to remove.')
        if HISTORY_INFO_SEPARATOR.strip() in self.raw_contents:
            content_without_last, *_ = self.raw_contents.rsplit(
                HISTORY_INFO_SEPARATOR.strip(),
                maxsplit=1
                )
        else:
            # Only one entry
            content_without_last = ''
        if content_without_last.endswith('\n'):
            content_without_last = content_without_last[:-len('\n')]
        # Clear file and write back entries
        self.path.write_text(content_without_last)
        self.read()  # Re-read to update last_entry

    def _parse_entry(self, entry_str):
        # remove leading and trailing whitespace
        entry_str = entry_str.strip()
        # check for 'DISCARDED' at the end
        if entry_str.endswith('DISCARDED'):
            discarded = True
            entry_str = entry_str[:-len('DISCARDED')]
        else:
            discarded = False
        # Notes may also be optional
        if 'Notes:' not in entry_str:
            notes = ''
            general_info = entry_str
        else:
            # split at 'Notes: ' because that is the only one that can be multiline
            general_info, notes = entry_str.split('Notes:', 1)
            notes = notes.replace('Notes: ', '').strip()
        # parse general_info
        general_info = iter(general_info.split('\n'))
        tensors = [int(num) for num in
                   next(general_info).replace('# TENSORS ', '').split()[1:]]
        jobs = [int(num) for num in
                next(general_info).replace('# JOB ID ', '').split()[1:]]
        the_rest = next(general_info)
        # optionals
        job_name, run_info, r_ref, r_super = None, None, None, None
        if the_rest.startswith('# JOB NAME '):
            job_name = the_rest.replace('# JOB NAME ', '').strip()
            the_rest = next(general_info)
        if the_rest.startswith('# RUN '):
            run_info = the_rest.replace('# RUN ', '').strip()
            the_rest = next(general_info)
        if the_rest.startswith('# TIME'):
            time = the_rest.replace('# TIME ', '').strip()
            the_rest = next(general_info)
        if the_rest.startswith('# R REF'):
            r_ref = float(the_rest.replace('# R REF ', ''))
            the_rest = next(general_info)
        if the_rest.startswith('# R SUPER'):
            r_super = float(the_rest.replace('# R SUPER ', ''))
            the_rest = next(general_info)
        folder = the_rest.replace('# FOLDER ', '').strip()

        # this should be all, if there is more, we have a problem
        try:
            # there may be empty lines at the end
            while not next(general_info):
                pass
        except StopIteration:
            pass
        else:
            raise ValueError(f'Error parsing {HISTORY_INFO_NAME} file.')

        return HistoryInfoEntry(tensors, jobs, time, folder, notes, discarded,
                                job_name, run_info, r_ref, r_super)


@dataclass
class HistoryInfoEntry:
    tensor_nums: List[int]
    job_nums: List[int]
    timestamp: str
    folder_name: str
    notes: str
    discarded: bool
    job_name: Optional[str] = None
    run_info: Optional[str] = None
    r_ref: Optional[float] = None
    r_super: Optional[float] = None


    def __str__(self):
        tensor_str = 'None' if self.tensor_nums == [0] else str(self.tensor_nums)[1:-1]
        job_str = str(self.job_nums)[1:-1]
        # translate timestamp if necessary
        time_str = (_translate_timestamp(self.timestamp)
                    if '-' in self.timestamp else self.timestamp)
        folder_str = self.folder_name
        return (
            '\n'
            + '# TENSORS '.ljust(_HISTORY_INFO_SPACING) + tensor_str + '\n'
            + '# JOB ID '.ljust(_HISTORY_INFO_SPACING) + job_str + '\n'
            + ('# JOB NAME '.ljust(_HISTORY_INFO_SPACING) + self.job_name + '\n'
                if self.job_name else '')
            + ('# RUN '.ljust(_HISTORY_INFO_SPACING) + self.run_info + '\n'
                if self.run_info else '')
            + '# TIME '.ljust(_HISTORY_INFO_SPACING) + time_str + '\n'
            + ('# R REF '.ljust(_HISTORY_INFO_SPACING) +f'{self.r_ref:.4f}\n'
                if self.r_ref is not None else '')
            + ('# R SUPER '.ljust(_HISTORY_INFO_SPACING) + f'{self.r_super:.4f}\n'
                if self.r_super is not None else '')
            + '# FOLDER '.ljust(_HISTORY_INFO_SPACING) + folder_str + '\n'
            + 'Notes: ' + self.notes.strip()
            + ('\nDISCARDED' if self.discarded else '')
            )


# This could probably be done by datetime.strptime
def _translate_timestamp(time_stamp):
    """Return a 'DD.MM.YY hh:mm:ss' timestamp from a YYMMDD-hhmmss one."""
    time_stamp = time_stamp.replace('moved-', '')
    if len(time_stamp) != 13:
        print(time_stamp)
        raise ValueError('Error translating timestamp: Invalid length '
                         f'{len(time_stamp)}. Expected 13 characters.')
    year, month, day = time_stamp[:2], time_stamp[2:4], time_stamp[4:6]
    hour, minutes, secs = time_stamp[7:9], time_stamp[9:11], time_stamp[11:13]
    return f'{day}.{month}.{year} {hour}:{minutes}:{secs}'
