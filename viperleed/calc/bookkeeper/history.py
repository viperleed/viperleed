"""Module history of viperleed.calc.bookkeeper."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

import ast
from dataclasses import dataclass, fields
from pathlib import Path
import re
from typing import List, Optional

from .constants import HISTORY_INFO_NAME
from .constants import LOGGER


_DISCARDED = 'DISCARDED'    # Label for entries marked via --discard
_HISTORY_INFO_SPACING = 12  # For the leftmost field in history.info
HISTORY_INFO_SEPARATOR = '\n###########\n'


class HistoryInfoError(Exception):
    """Base class for all errors related to the history.info file."""


class NoHistoryEntryError(HistoryInfoError):
    """There is no entry to process according to the criteria."""


class EntrySyntaxError(HistoryInfoError):
    """Something is wrong with the syntax of an entry."""


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
        if self.last_entry.discarded:
            LOGGER.warning('Last entry is already discarded.')
            return
        discarded_entry = self.last_entry.as_discarded()
        self.remove_last_entry()
        self.append_entry(discarded_entry)

    def read(self):
        """Read the current contents of the history.info file."""
        with open(self.path, 'r', encoding='utf-8') as history_info:
            self.raw_contents = history_info.read()
        if not self.raw_contents:
            self.last_entry = None
        else:
            self.last_entry = HistoryInfoEntry.from_string(
                self.raw_contents.split(HISTORY_INFO_SEPARATOR.strip())[-1]
                )

    def remove_last_entry(self):
        """Discard the last entry from the history.info file."""
        if self.last_entry is None:
            raise NoHistoryEntryError('No entries to remove.')
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


# Map of HistoryInfoEntry-field names to tags in history.info file
_TAG = {
    'tensor_nums': '# TENSORS',
    'job_nums': '# JOB ID',
    'timestamp': '# TIME',
    'folder_name': '# FOLDER',
    'notes': 'Notes:',
    'job_name': '# JOB NAME',
    'run_info': '# RUN',
    'r_ref': '# R REF',
    'r_super': '# R SUPER',
    }

# Regular expression portions for parsing a history entry
_SPACE = r'[ \t]'  # Exclude \n
_COMMA_SEPARATED = rf'\d+({_SPACE}*,{_SPACE}*\d+)*'
_SPACE_SEPARATED = rf'\d+({_SPACE}+\d+)*'
_FLOAT_RE = r'\d+(.\d+)?'
_FOLDER_RE = (rf'{_TAG["folder_name"]}{_SPACE}+'
              r'(?P<folder_name>t\d+\.r\d+_\d{6,6}-\d{6,6})')
_JOB_ID_RE = rf'{_TAG["job_nums"]}{_SPACE}+(?P<job_nums>({_COMMA_SEPARATED})?)'
_JOB_NAME_RE = rf'{_TAG["job_name"]}{_SPACE}+(?P<job_name>.+)'
_NOTES_RE = rf'{_TAG["notes"]}{_SPACE}*(?P<notes>.*(\n.*)*)'
_RUN_ID_RE = rf'{_TAG["run_info"]}{_SPACE}+(?P<run_info>({_SPACE_SEPARATED})?)'
_R_REF_RE = rf'{_TAG["r_ref"]}{_SPACE}+(?P<r_ref>{_FLOAT_RE})'
_R_SUPER_RE = rf'{_TAG["r_super"]}{_SPACE}+(?P<r_super>{_FLOAT_RE})'
_TENSORS_RE = (rf'{_TAG["tensor_nums"]}{_SPACE}+'
               rf'(?P<tensor_nums>({_COMMA_SEPARATED}|None)?)')
_TIME_RE = rf'{_TAG["timestamp"]}{_SPACE}+(?P<timestamp>[-:. \d]+)'


# Make all entries optional and complain later on if the mandatory
# fields are missing. It probably means users have modified stuff.
_ENTRY_RE = re.compile(
    fr'({_TENSORS_RE})?'
    fr'(\n{_JOB_ID_RE})?'
    fr'(\n{_JOB_NAME_RE})?'
    fr'(\n{_RUN_ID_RE})?'
    fr'(\n{_TIME_RE})?'
    fr'(\n{_R_REF_RE})?'
    fr'(\n{_R_SUPER_RE})?'
    fr'(\n{_FOLDER_RE})?'
    fr'(\n{_NOTES_RE})?'
    )


@dataclass(frozen=True)
class HistoryInfoEntry:  # pylint: disable=R0902  # See pylint #9058
    """A container for information in a single "block" of history.info.

    Attributes
    ----------
    tensor_nums : list of int
        The progressive identifiers of Tensors used during this run.
    job_nums : list of int
        The progressive identifiers of the executions that used
        certain Tensors during this run.
    timestamp : str
        The time when this run was started.
    folder_name : str
        The name of the history folder that corresponding to this run.
    notes : str
        The notes that users have added for
        this run. May span multiple lines.
    discarded : bool
        Whether the user decided to mark this run as DISCARDED.
    job_name : str, optional
        A name assigned to this run by the user.
    run_info : str, optional
        The sequence of segments that were executed.
    r_ref : float, optional
        The R factor resulting from a reference-calculation run.
    r_super : float, optional
        The R factor resulting from a superpos-calculation run.
    """

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
        """Return a string version of this entry, ready for writing to file."""
        no_tensors = not self.tensor_nums or self.tensor_nums == [0]
        tensor_str = 'None' if no_tensors else str(self.tensor_nums)[1:-1]
        job_str = str(self.job_nums)[1:-1]
        # translate timestamp if necessary
        time_str = (_translate_timestamp(self.timestamp)
                    if '-' in self.timestamp else self.timestamp)
        return (
            '\n'
            + self._format_field('tensor_nums', value_str=tensor_str)
            + self._format_field('job_nums', value_str=job_str)
            + self._format_field('job_name')
            + self._format_field('run_info')
            + self._format_field('timestamp', value_str=time_str)
            + self._format_field('r_ref', fmt='.4f')
            + self._format_field('r_super', fmt='.4f')
            + self._format_field('folder_name')
            + f'{_TAG["notes"]} {self.notes.strip()}'
            + (f'\n{_DISCARDED}' if self.discarded else '')
            )

    def as_discarded(self):
        """Return a discarded version of this entry."""
        if self.discarded:
            return self
        cls = type(self)
        kwargs = {f.name: getattr(self, f.name) for f in fields(self)}
        kwargs['discarded'] = True
        return cls(**kwargs)

    @classmethod
    def from_string(cls, entry_str):
        """Return an HistoryInfoEntry instance from its string version."""
        match = _ENTRY_RE.match(entry_str.strip())
        if not match:
            raise EntrySyntaxError(f'Could not parse entry:\n{entry_str}')

        kwargs = {'discarded': False}  # Deal with discarded later
        for field in fields(cls):
            name = field.name
            if name in kwargs:
                continue
            tag = _TAG[name].replace('# ', '')
            value = match[name]
            try:
                kwargs[name] = cls._field_from_str(field, value)
            except (ValueError, TypeError, SyntaxError,
                    MemoryError, RecursionError) as exc:
                raise EntrySyntaxError('Could not parse entry field '
                                       f'{tag!r} with value={value!r}. '
                                       f'Info: {exc}') from exc
        # Now deal with "DISCARDED"
        notes = kwargs['notes']
        if notes is not None:
            kwargs['discarded'] = notes.endswith(_DISCARDED)
            if kwargs['discarded']:
                kwargs['notes'] = notes[:-len(_DISCARDED)].rstrip()
        return cls(**kwargs)

    @staticmethod
    def _field_from_str(field, value_str):
        """Return a field of the right type from a its string value.

        Parameters
        ----------
        field : datclasses.Field
            The field whose value should be returned.
        value_str : str or None
            The string version of the `field` value.
            It may be `None` for optional fields.

        Returns
        -------
        field_value : object
            The value for field, as parsed from `value_str`. It may
            be None for an optional field if `value_str` is None.

        Raises
        ------
        ValueError, TypeError, SyntaxError, MemoryError, RecursionError
            If parsing `value_str` via ast.literal_eval fails. Only
            `value_str` for non-string fields are parsed like this.
        """
        # The following trick works because Optional removes 'repeated'
        # entries, so that Optional[Optional[t]] == Optional[t]
        is_optional = field.type == Optional[field.type]
        if is_optional and value_str is None:
            return value_str
        is_optional_str = field.type == Optional[str]
        if is_optional_str or field.type is str:
            return value_str
        if field.type.__origin__ is list and not value_str:
            return []
        if field.type.__origin__ is list:
            return list(ast.literal_eval(value_str + ','))
        return ast.literal_eval(value_str)

    @staticmethod
    def _format_line(tag, value):
        """Return a formatted line from a tag and values."""
        return f'{tag:<{_HISTORY_INFO_SPACING}}{value}\n'

    def _format_field(self, field_name, value_str=None, fmt=''):
        """Return a formatted version of a field."""
        if value_str is None:
            value = getattr(self, field_name)
            if not value:
                return ''
            value_str = f'{value:{fmt}}'
        return self._format_line(_TAG[field_name], value_str)


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
