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
from collections.abc import Sequence
from copy import deepcopy
from dataclasses import dataclass
from dataclasses import field as data_field
from dataclasses import fields
from dataclasses import replace as replace_value
from datetime import datetime
from enum import Enum
from functools import partial
from functools import partialmethod
from pathlib import Path
import re
from typing import Dict, List, Optional, Union

from viperleed.calc.lib.base import logging_silent

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


class FixableSyntaxError(HistoryInfoError):
    """An entry has some syntax error, but it can be fixed."""


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
    def last_entry_has_notes(self):
        """Return whether the last entry has associated notes."""
        if self.last_entry is None:
            return False
        return bool(self.last_entry.has_notes)

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
        self.raw_contents = self.path.read_text(encoding='utf-8')
        if not self.raw_contents:
            self.last_entry = None
            return
        self._infer_time_format()
        *_, entry_str = self.raw_contents.split(HISTORY_INFO_SEPARATOR.strip())
        entry = HistoryInfoEntry.from_string(entry_str)                         # TODO: how to test this properly? It does not raise. Perhaps a TolerantHistoryInfoEntry?
        if self._time_format:
            entry = entry.with_time_format(self._time_format)
        self.last_entry = entry

    def remove_last_entry(self):
        """Remove the last entry from the history.info file."""
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
                self._time_format = entry.time_format
                return


# Map of HistoryInfoEntry-field names to tags in history.info file,
# SORTED in the expected order. !!! Do not change sorting !!!
_TAG = {
    'tensor_nums': '# TENSORS',
    'job_nums': '# JOB ID',
    'job_name': '# JOB NAME',
    'run_info': '# RUN',
    'timestamp': '# TIME',
    'r_ref': '# R REF',
    'r_super': '# R SUPER',
    'folder_name': '# FOLDER',
    'notes': 'Notes:',
    }


# Regular expression portions for parsing a history entry
_SPACE = r'[ \t]'  # Exclude \n
_COMMA_SEPARATED = rf'\d+({_SPACE}*,{_SPACE}*\d+)*'
_COMMA_OR_SPACE_SEPARATED = rf'\d+({_SPACE}*[,{_SPACE}]?{_SPACE}*\d+)*'
_FLOAT_RE = r'\d+(.\d+)?'  # Positive, without '+'
_MULTILINE = r'.*(\n.*)*'
_SPACE_SEPARATED = rf'\d+({_SPACE}+\d+)*'

# Some special ones for the fields
# Notice the double escape of curly braces that we want to appear as
# literal braces: we .format these patterns later on.
# In principle, the format may be a bit more complex (see  bookkeeper           TODO: test
# _get_new_history_directory_name), but I (@michele-riva) think that
# the v1.0.0 modes make the alternative names obsolete.
_FOLDER_NAME_RE = r't\d+\.r\d+_(moved-)?\d{{6}}-\d{{6}}{job_name}'

# Regular expression for matching whole history entries:
# Catch any character past the expected tags (notes: over more lines).
# HistoryInfoEntry deals with the details of each field.
_FIELDS_RE = {field: rf'{tag}{_SPACE}*(?P<{field}>.*)'
              for field, tag in _TAG.items()}
_FIELDS_RE['notes'] = _FIELDS_RE['notes'].replace('.*', _MULTILINE)             # TODO: marked as untested??
_ENTRY_RE = re.compile(
    # All fields optional, just in the right order.
    # HistoryInfoEntry will deal with missing ones
    ''.join(rf'(\n?{field_re})?' for field_re in _FIELDS_RE.values())
    )


# Sentinel for fields that do not exist at all.


# Some standard error messages
_MSG_NOT_UNDERSTOOD_PREFIX = f'{HISTORY_INFO_NAME}: Could not understand'
_MSG_EMPTY = 'Field is empty'
_MSG_NO_FLOAT = 'Field is not a floating-point number'
_MSG_NO_STRING = 'Field is not a string'
_MSG_MISSING = 'Mandatory field not found'


class TimestampFormat(Enum):
    """A collection of known formats for timestamps."""

    # Underscore names are used only for conversion purposes
    # while parsing strings, and never for formatting.
    GERMAN = '%d.%m.%y %H:%M:%S'  # Format used for < v1.0.0
    ISO = '%Y-%m-%d %H:%M:%S'     # Format used for >= v1.0.0
    _CALC = '%y%m%d-%H%M%S'       # Format used by calc in log files
    DEFAULT = ISO

    @property
    def writable(self):
        """Return whether this format can be used for writing."""
        return not self.name.startswith('_')


def is_optional(field, with_type=None):
    """Return whether a datclasses.Field is optional."""
    if with_type is None:
        with_type = field.type
    # The following trick works because Optional removes 'repeated'
    # entries, so that Optional[Optional[t]] == Optional[t]
    return field.type == Optional[with_type]


# TODO: the _fixup methods may be made even more lenient
# and accept some comment text following each field.
# TODO: this class is a bit messy. Perhaps one could have better code
# with an HistoryInfoEntryField base class (of which all fields are
# subclasses) that collects the base functionality, and subclasses
# for the different fields? Each field could then have a .format and
# a .check method that do what is done here now. It may a bit tricky
# to get it right, though.
# TODO: _MISSING notes and not discarded should not print, but
# missing notes with discarded should.
# TODO: comment-only entry
@dataclass(frozen=True)
class HistoryInfoEntry:  # pylint: disable=R0902  # See pylint #9058
    """A container for information in a single "block" of history.info.

    Attributes
    ----------
    tensor_nums : list or str
        The progressive identifiers of Tensors used during this
        run. Also string is tolerated. A string containing a comma-
        or space-separated list of non-negative integers is also
        accepted at object creation. A space-separated list will
        cause the entry to be marked as 'needs fixing'. A string
        that cannot be converted to a list of integers is accepted
        and retained, but causes the entry to be flagged as 'not
        understood' (and thus not fixable). None, 'None', or [None]
        are accepted at object creation, and interpreted as if no
        Tensors were used or produced.
    job_nums : list or str
        The progressive identifiers of the executions that used
        certain Tensors during this run. Also string is tolerated.
        A string containing a comma- or space-separated list of
        non-negative integers is also accepted at object creation.
        A space-separated list will cause the entry to be marked as
        'needs fixing'. A string that cannot be converted to a list
        of integers is accepted, but causes the entry to be flagged
        as 'not understood' (and thus not fixable).
    timestamp : datetime or str
        The time when this run was started. A string representing
        a date-time is also accepted. 'moved-' is discarded from
        the string before parsing. If a valid string is given, but
        does not conform with TimestampFormat.DEFAULT, the value
        is accepted but marked as 'needs fixing'. A string that
        cannot be turned into a date-time is also accepted, but
        causes the entry to be flagged as 'not understood' (and
        thus not fixable).
    folder_name : str
        The name of the history folder that corresponding to this run.
    notes : str, optional
        The notes that users have added for this run. It excludes the
        DISCARDED tag. May span multiple lines. The default is an empty
        string.
    discarded : bool, optional
        Whether the user decided to mark this run as DISCARDED.
        Default is False.
    job_name : str, optional
        A name assigned to this run by the user.
    run_info : str, optional
        The sequence of segments that were executed. Expected format:
        space-separated list of integers. Any non-conforming string
        marks the field as 'not understood' (and thus not fixable).
    r_ref : float, optional
        The R factor resulting from a reference-calculation run.
        Numbers outside the interval [0, 2] are flagged as 'not
        understood'. A string is also accepted at object creation.
        In this case, it is considered acceptable if it can be
        converted to a valid R factor. An invalid string is also
        retained, but will cause the field to be marked as 'not
        understood'.
    r_super : float, optional
        The R factor resulting from a superpos-calculation run.
        Numbers outside the interval [0, 2] are flagged as 'not
        understood'. A string is also accepted at object creation.
        In this case, it is considered acceptable if it can be
        converted to a valid R factor. An invalid string is also
        retained, but will cause the field to be marked as 'not
        understood'.
    """

    tensor_nums: Union[List[int], str] = _MISSING
    job_nums: Union[List[int], str] = _MISSING
    timestamp: Union[datetime, str] = _MISSING
    folder_name: str = _MISSING
    notes: Optional[str] = ''
    discarded: Optional[bool] = False
    job_name: Optional[str] = None
    run_info: Optional[str] = None
    r_ref: Optional[float] = None
    r_super: Optional[float] = None
    # _needs_fixup is {field_name: (reason, fixed)} for all the
    # fields that do not conform to the standard formatting. The
    # fixed item is the attribute value after fixing.
    _needs_fixup: Dict[str, tuple] = data_field(default_factory=dict,
                                                init=False)
    # _not_understood has a similar format (albeit only with reason)
    # for those fields that do not make sense.
    _not_understood: Dict[str, str] = data_field(default_factory=dict,
                                                 init=False)
    # _formatters are used for __str__; Cached for quicker access
    _formatters: Dict[str, callable] = data_field(default_factory=dict,
                                                  init=False)
    _time_format: TimestampFormat = data_field(default=None, init=False)

    __hash__ = None  # Not hash-able, just to be sure

    def __post_init__(self):
        """Determine whether arguments make sense and if they need fixing."""
        for field in fields(self):
            name = field.name
            checker = getattr(self, f'_check_{name}_field', None)
            if not checker:
                continue
            tag = _TAG[name].replace('# ', '')
            try:
                checker()
            except FixableSyntaxError as reason:
                msg = f'Found entry with {reason}.'
                if tag not in msg:
                    msg += f' Field: {tag}.'
                LOGGER.warning(f'{HISTORY_INFO_NAME}: {msg} Consider '
                               'running bookkeeper in --fixup mode.')
            except EntrySyntaxError as exc:
                self._not_understood[name] = str(exc)
                value = getattr(self, name)
                value_msg = (f' with value {value!r}' if value != _MISSING
                             else '')
                LOGGER.warning(
                    f'{HISTORY_INFO_NAME}: Could not understand {tag} '
                    f'field{value_msg}. Reason: {exc}'
                    )
        if not self.was_understood:
            LOGGER.warning(f'{HISTORY_INFO_NAME}: Faulty entry is\n\n%s',
                           self.format_problematic_fields())

    def __str__(self):
        """Return a string version of this entry, ready for writing to file."""
        self._collect_formatters()
        formatted = (fmt() for field, fmt in self._formatters.items()
                     if getattr(self, field) != _MISSING)
        txt = '\n'.join(f for f in formatted if f.rstrip())
        return '\n' + txt

    @property
    def can_be_removed(self):
        """Return whether this entry can be removed from history.info."""
        return (self.was_understood
                and not self.needs_fixing
                and not self.has_notes)

    @property
    def has_notes(self):
        """Return whether this entry has user notes."""
        return self.notes != _MISSING and self.notes

    @property
    def misses_mandatory_fields(self):
        """Return whether this entry has missing fields."""
        return any(getattr(self, f.name) == _MISSING
                   for f in fields(self)
                   if f.init and not is_optional(f))

    @property
    def needs_fixing(self):
        """Return whether any of the entry fields are wrongly formatted."""
        return any(self._needs_fixup)

    @property
    def time_format(self):
        """Return the format used for self.timestamp."""
        return self._time_format or TimestampFormat.DEFAULT

    @property
    def was_understood(self):
        """Return whether all fields in this entry were recognized."""
        return not any(self._not_understood)

    def as_discarded(self):
        """Return a discarded version of this entry."""
        if self.discarded:
            return self
        # replace_value goes via __init__ + __post_init__, and would
        # log the same messages as for this instance. Mute logger.
        with logging_silent():
            return replace_value(self, discarded=True)

    def as_fixed(self):
        """Return a version of this entry with fixable fields replaced."""
        if not self.needs_fixing:
            return self
        kwargs = {attr: value  # Fix only those that were understood
                  for attr, (_, value) in self._needs_fixup.items()}
        # The timestamp needs special care, as there is no
        # 'fixed' value, but it is a matter of format change
        fix_time_fmt = kwargs.pop('timestamp', False)

        # replace_value goes via __init__ + __post_init__, and would
        # log the same messages as for this instance. Mute logger.
        with logging_silent():
            fixed = replace_value(self, **kwargs)
            if fix_time_fmt:
                fixed = fixed.with_time_format(TimestampFormat.DEFAULT)
        if fixed.needs_fixing:  # Safeguard for the future
            raise HistoryInfoError('Failed to fix fields for entry:\n'
                                   + self.format_problematic_fields())
        return fixed

    def format_problematic_fields(self):
        """Return a string with problematic fields highlighted."""
        if not self.needs_fixing and self.was_understood:
            return ''

        def _get_prefix(field):
            if getattr(self, field) == _MISSING:
                return 'missing -> '
            if field in self._needs_fixup:
                return 'fixable -> '
            if field in self._not_understood:
                return 'edited? -> '
            return ' '*11

        self._collect_formatters()
        return ''.join(
            _get_prefix(field) + line + '\n'
            for field, formatter in self._formatters.items()
            for line in formatter().splitlines()
            )

    @classmethod
    def from_string(cls, entry_str):
        """Return an HistoryInfoEntry instance from its string version."""
        cls._from_string_check_single_entry(entry_str)
        kwargs = {'discarded': False}      # Later, as we need notes
        cls._parse_fields_from_string(entry_str, kwargs)
        cls._process_string_notes(kwargs)  # Notes and DISCARD
        return cls(**kwargs)

    def with_time_format(self, fmt):
        """Return a version of this entry with a specific timestamp format."""
        known_fmts = ', '.join(repr(f.name) for f in TimestampFormat
                               if f.writable)
        if isinstance(fmt, str):
            fmt = fmt.upper()
            try:
                fmt = TimestampFormat[fmt]
            except KeyError:
                raise ValueError(f'Unknown timestamp format {fmt!r}. '
                                 f'Should be one of {known_fmts}') from None
        if not isinstance(fmt, TimestampFormat):
            raise TypeError(f'Expected {known_fmts}, or TimestampFormat. '
                            f'Found {type(fmt).__name__!r}')
        if fmt is self._time_format:
            return self
        new = deepcopy(self)
        object.__setattr__(new, '_time_format', fmt)
        # And clean up possible time format issues
        if fmt is fmt.DEFAULT:
            # pylint: disable-next=protected-access  # It's a deepcopy
            new._needs_fixup.pop('timestamp', None)
        return new

    def _check_folder_name_field(self):
        """Check whether the folder_name attribute is fine."""
        attr = 'folder_name'
        self._check_missing(attr)
        self._check_mandatory_field_is_not_empty(attr)
        value = getattr(self, attr)
        if not isinstance(value, str):
            raise EntrySyntaxError(_MSG_NO_STRING)

        job_name = re.escape(
            '' if not self.job_name or not isinstance(self.job_name, str)
            else self.job_name.strip()
            )
        suffix = f'_{job_name}' if job_name else ''
        pattern = _FOLDER_NAME_RE.format(job_name=suffix)
        folder_re = re.compile(pattern)
        match_leading = folder_re.match(value)
        if not match_leading:
            raise EntrySyntaxError(
                'Does not have format '
                'tTTT.rRRR_[moved-]yymmdd-HHMMSS[_job_name]'
                )
        if not folder_re.fullmatch(value):
            # May be fixed by stripping the part that does not match
            fixup_value = match_leading.group()
            reason = f'unexpected text after {fixup_value!r}'
            self._needs_fixup[attr] = (reason, fixup_value)
            raise FixableSyntaxError(reason)

    def _check_job_nums_field(self):
        """Check whether the job_nums attribute is fine."""
        self._check_list_of_int_field('job_nums', accept_none=False)

    def _check_mandatory_field_is_not_empty(self, attr):
        """Complain if a mandatory field is empty."""
        value = getattr(self, attr)
        if isinstance(value, str) and not value.strip():
            object.__setattr__(self, attr, '')
            raise EntrySyntaxError(_MSG_EMPTY)

    def _check_missing(self, attr):
        """Complain if no value was provided at all for `attr`."""
        if getattr(self, attr) == _MISSING:
            raise EntrySyntaxError(_MSG_MISSING)

    def _check_r_factor_field(self, which_r):
        """Check whether one of the R-factor fields is fine."""
        value = getattr(self, which_r)
        if value is None:  # R factors are optional
            return
        if isinstance(value, str) and not value.strip():
            # Also OK to have an empty string for optional fields
            object.__setattr__(self, which_r, None)
            return
        if isinstance(value, str):  # Non-empty
            value = self._fixup_float_string(value)
            object.__setattr__(self, which_r, value)
        if not isinstance(value, float):
            raise EntrySyntaxError(_MSG_NO_FLOAT)
        if not -1e-5 < value < 2+1e-5:
            raise EntrySyntaxError('Number is not between zero and two')

    def _check_r_ref_field(self):
        """Check whether the r_ref attribute is fine."""
        self._check_r_factor_field('r_ref')

    def _check_r_super_field(self):
        """Check whether the r_super attribute is fine."""
        self._check_r_factor_field('r_super')

    def _check_run_info_field(self):
        """Check whether the run_info attribute is fine."""
        value = self.run_info
        if not value:  # OK, it's optional
            return
        if not isinstance(self.run_info, str):
            raise EntrySyntaxError(_MSG_NO_STRING)
        value = value.strip()
        if not value:  # OK, it's optional
            return
        if not re.fullmatch(_SPACE_SEPARATED, value):
            raise EntrySyntaxError('Not a space-separated list of integers.')

    def _check_tensor_nums_field(self):
        """Check whether the tensor_nums attribute is fine."""
        self._check_list_of_int_field('tensor_nums', accept_none=True)

    def _check_list_of_int_field(self, attr, accept_none=False):
        """Check whether the value of `attr` is a valid list of integers."""
        self._check_missing(attr)

        # Below we use object.__setattr__ to circumvent the frozen
        # dataclass. Notice that this is not a problem, as we call
        # this method only in __post_init__, and __hash__ uses a
        # tuple of the values of the attributes.
        self._check_mandatory_field_is_not_empty(attr)
        value = getattr(self, attr)
        is_akin_to_none = (  # Variations of None
            value is None
            or value == [None]
            or (isinstance(value, str) and value.strip() == str(None))
            )
        if accept_none and is_akin_to_none:
            object.__setattr__(self, attr, [])
            return
        if not accept_none and not value:
            object.__setattr__(self, attr, [])
            raise EntrySyntaxError(_MSG_EMPTY)
        if isinstance(value, str):
            self._fixup_list_of_int_string(attr, value)
            value = getattr(self, attr)

        # Check that it is a sequence
        if not isinstance(value, Sequence):
            raise EntrySyntaxError('Not a sequence')
        # Make sure it is indeed a list of integers
        if not all(isinstance(v, int) for v in value):
            raise EntrySyntaxError('Not a list of integers')
        # An finally, that they are non-negative
        if not all(v >= 0 for v in value):
            raise EntrySyntaxError('Contains negative values')

    def _check_timestamp_field(self):
        """Check whether the timestamp attribute is fine."""
        self._check_missing('timestamp')
        value = self.timestamp
        if isinstance(value, datetime):  # OK
            return
        if not isinstance(value, str):
            raise EntrySyntaxError(_MSG_NO_STRING)

        value = value.replace('moved-', '').strip()
        if not value:
            object.__setattr__(self, 'timestamp', '')
            raise EntrySyntaxError(_MSG_EMPTY)

        for fmt in TimestampFormat:
            try:
                object.__setattr__(self,
                                   'timestamp',
                                   datetime.strptime(value, fmt.value))
            except ValueError:  # Not the right format
                continue
            # Conversion successful. Remember the format,
            # as long as it is acceptable for writing
            if fmt.writable:
                object.__setattr__(self, '_time_format', fmt)
            break
        else:  # No valid format
            raise EntrySyntaxError('Field could not be parsed as a date-time')
        # Finally, check whether the format is the current default
        if self.time_format is not TimestampFormat.DEFAULT:
            reason = 'outdated format'
            # No fixed value as one has to fiddle with the _time_format
            # attribute instead. This is done in with_time_format().
            self._needs_fixup['timestamp'] = (reason, None)
            raise FixableSyntaxError(reason)

    def _collect_formatters(self):
        """Update `self._formatters` if it was not done already."""
        if self._formatters:
            return
        formatters = {
            field: getattr(self,
                           f'_format_{field}',
                           partial(self._format_field, field_name=field))
            for field in _TAG
            }
        object.__setattr__(self, '_formatters', formatters)

    def _fixup_float_string(self, value_str):
        """Try converting `value_str` to a float value."""
        value_str = value_str.strip()
        if not re.fullmatch(_FLOAT_RE, value_str):
            raise EntrySyntaxError(_MSG_NO_FLOAT)
        return ast.literal_eval(value_str)  # Should never fail

    def _fixup_list_of_int_string(self, attr, value_str):
        """Try converting `value_str` to a list of int for `attr`."""
        value_str = value_str.strip()
        if not value_str:
            raise EntrySyntaxError(_MSG_EMPTY)
        if not re.fullmatch(_COMMA_OR_SPACE_SEPARATED, value_str):
            raise EntrySyntaxError('Not a list of comma-separated integers')

        comma_separated = re.fullmatch(_COMMA_SEPARATED, value_str)
        if not comma_separated:
            # Make it a comma-separated list of int.
            # First get rid of spaces before commas
            value_str = re.sub(r'\s+,', ',', value_str)
            # Then replace all multi-spaces not preceded by comma
            value_str = re.sub(r'(?<!,)\s+', ',', value_str)

        value = list(ast.literal_eval(value_str + ','))
        if comma_separated:  # Format OK. Store list instead of str
            object.__setattr__(self, attr, value)
            return
        reason = 'some space-separated items, rather than comma separated'
        self._needs_fixup[attr] = (reason, value)
        raise FixableSyntaxError(reason)

    def _format_field(self, field_name, value_str=None, fmt=''):
        """Return a formatted version of a field."""
        value = getattr(self, field_name)
        if value == _MISSING:
            value_str = ''
        elif value_str is None and not value:
            return ''
        elif value_str is None:
            value_str = value if isinstance(value, str) else f'{value:{fmt}}'
        return self._format_line(_TAG[field_name], value_str)

    def _format_job_nums(self):
        """Return a line with JOB ID."""
        # Important: must check str first, as this is the best way
        # to retain as much as possible the value that users may
        # have modified.
        if isinstance(self.job_nums, str):
            jobs_str = self.job_nums
        elif not self.job_nums:
            jobs_str = ''
        else:
            assert isinstance(self.job_nums, Sequence)
            jobs_str = str(list(self.job_nums))[1:-1]
        return self._format_field('job_nums', value_str=jobs_str)

    @staticmethod
    def _format_line(tag, value):
        """Return a formatted line from a tag and values."""
        return f'{tag:<{_HISTORY_INFO_SPACING}}{value}'

    def _format_notes(self):
        """Return a string with multi-line notes."""
        notes_str = self.notes.strip()
        if notes_str == _MISSING:
            return ''
        notes = f'{_TAG["notes"]} {notes_str}'
        discarded = f'{_DISCARDED}' if self.discarded else ''
        if discarded and notes_str:
            notes += '\n'
        return notes + discarded

    _format_r_super = partialmethod(_format_field,
                                    field_name='r_super',
                                    fmt='.4f')
    _format_r_ref = partialmethod(_format_field, field_name='r_ref', fmt='.4f')

    def _format_tensor_nums(self):
        """Return a line with TENSORS."""
        value = self.tensor_nums
        no_tensors = (not value         # Field empty
                      or value == [0])  # INIT-only run
        # Important: must check str first, as this is the best way
        # to retain as much as possible the value that users may
        # have removed
        if isinstance(value, str):
            tensor_str = value
        elif no_tensors:
            tensor_str = 'None'
        else:  # A sequence. Use list as it does not add trailing comma
            assert isinstance(value, Sequence)
            tensor_str = str(list(value))[1:-1]
        return self._format_field('tensor_nums', value_str=tensor_str)

    def _format_timestamp(self):
        """Return a line with TIME."""
        value = self.timestamp
        time_str = (value if isinstance(value, str)
                    else value.strftime(self.time_format.value))
        return self._format_field('timestamp', value_str=time_str)

    @staticmethod
    def _from_string_check_single_entry(entry_str):
        """Make sure that `entry_str` only contains a single entry."""
        if HISTORY_INFO_SEPARATOR.strip() in entry_str:                         # TODO: This, and the splitting we do above, does not cover cases in which users may have used the separator in other spots.
            raise ValueError('Found multiple entries. Split them at '
                             'HISTORY_INFO_SEPARATOR beforehand')
        err_msg = ('Found multiple entries without a separator (i.e., '
                   f'{HISTORY_INFO_SEPARATOR.strip()!r}). Did you delete '
                   f'some? Problematic entry:\n{entry_str}')
        # There are two options: if there is a 'Notes:' field, its
        # multi-line regex will gobble up everything that follows,
        # including following entries. This is not a problem if
        # the 'Notes:' field has been deleted. Do the last one first.
        n_entries = sum(bool(m.group(0))
                        for m in _ENTRY_RE.finditer(entry_str))
        if n_entries > 1:
            raise EntrySyntaxError(err_msg)

        # Now check the second option. Here we have to take off one
        # line at a time from the notes, and check if the rest matches
        notes = _ENTRY_RE.match(entry_str)['notes']
        if notes is None:  # OK
            return
        # Add an extra line to check all the lines in notes
        notes = notes.strip() + '\n' + _MISSING
        while notes and notes != _MISSING:
            match = _ENTRY_RE.match(notes.strip())
            if match and match.group(0):
                raise EntrySyntaxError(err_msg)
            # The extra _MISSING line and the while
            # condition prevent a ValueError here:
            _, notes = notes.split('\n', maxsplit=1)
            notes = notes.strip()

    @classmethod
    def _parse_fields_from_string(cls, entry_str, kwargs):
        """Update `kwargs` by parsing `entry_str`."""
        match = _ENTRY_RE.match(entry_str.strip())
        if not match:
            raise HistoryInfoError(f'Could not parse entry\n{entry_str}\n'
                                   'This is probably a bug. Please report it '
                                   'to the ViPErLEED team')
        for field in fields(cls):
            if not field.init:
                continue
            name = field.name
            known_value = kwargs.get(name)
            if known_value is not None:
                continue
            value = match[name]
            if value is not None:
                kwargs[name] = value.rstrip()
            elif not is_optional(field):  # No match at all
                kwargs[name] = _MISSING

    @staticmethod
    def _process_string_notes(kwargs):
        """Update `kwargs` according to its notes."""
        notes = kwargs['notes']
        if notes is None:
            kwargs['notes'] = _MISSING
            return
        kwargs['discarded'] = notes.endswith(_DISCARDED)
        if kwargs['discarded']:
            kwargs['notes'] = notes[:-len(_DISCARDED)].rstrip()
