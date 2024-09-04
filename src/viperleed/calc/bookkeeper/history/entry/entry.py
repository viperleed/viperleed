"""Module entry of viperleed.calc.bookkeeper.history.entry.

Defines classes that handle the contents of a single 'block'
of the history.info file.
"""

__authors__ = (
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
from functools import partial
from functools import partialmethod
import re
from typing import Dict, List, Optional, Union

from viperleed.calc.bookkeeper.log import LOGGER
from viperleed.calc.lib.dataclass_utils import is_optional_field
from viperleed.calc.lib.dataclass_utils import set_frozen_attr
from viperleed.calc.lib.log_utils import logging_silent

from ..constants import HISTORY_INFO_NAME
from ..constants import HISTORY_INFO_SEPARATOR
from ..errors import EntrySyntaxError
from ..errors import FixableSyntaxError
from ..errors import HistoryInfoError
from ..errors import _PureCommentEntryError
from .time_field import TimestampFormat


_DISCARDED = 'DISCARDED'    # For entries marked via --discard
_HISTORY_INFO_SPACING = 12  # For the leftmost field in an entry


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
_FLOAT_RE = r'\d+(.\d+)?'    # Positive, without '+'
_MULTILINE = r'.*(?:\n.*)*'  # Non capturing to detect unknown lines
_SPACE_SEPARATED = rf'\d+({_SPACE}+\d+)*'

# Some special ones for the fields
# Notice the double escape of curly braces that we want to appear as
# literal braces: we .format these patterns later on.
# In principle, the format may be a bit more complex (see  bookkeeper           TODO: test
# _get_new_history_directory_name), but I (@michele-riva) think that
# the v1.0.0 modes make the alternative names obsolete.
_FOLDER_NAME_RE = r't\d+\.r\d+_(moved-)?\d{{6}}-\d{{6}}{job_name}'


def _make_entry_regex():
    """Return a compiled regular expression for a history.info entry."""
    # Regular expression for matching whole history entries:
    # Catch any character past the expected tags (notes: over more
    # lines). HistoryInfoEntry deals with the details of each field.
    fields_re = {field: rf'{tag}{_SPACE}*(?P<{field}>.*)'
                 for field, tag in _TAG.items()}
    fields_re['notes'] = fields_re['notes'].replace('.*', _MULTILINE)
    fields_re_keys = tuple(fields_re)

    def _except_next(ind=None):
        """Return a pattern for anything not starting with known keys."""
        all_next = (
            '|'.join(_TAG[item] for item in fields_re_keys) if ind is None
            else '|'.join(_TAG[item] for item in fields_re_keys[ind+1:])
            )
        return rf'(\n?(?!{all_next}).+)*'

    # Prepare pattern for a single entry:
    # Can start with anything, except the tag of the first field,
    # then have all fields in the right order, potentially with
    # stuff in between. All fields are optional and HistoryInfoEntry
    # deals with missing and extra ones.
    entry_re_pattern = _except_next()
    for i, field_re in enumerate(fields_re.values()):
        # Non capturing (?:...) is necessary to detect unknown lines
        entry_re_pattern += rf'(?:\n?{field_re})?'
        entry_re_pattern += rf'(\n?(?!{_except_next(i)}).+)*'
    return re.compile(entry_re_pattern)

_ENTRY_RE = _make_entry_regex()


# Sentinel for fields that do not exist at all.
_MISSING = '_tag_does_not_exist_in_entry__#oMBde$fsT&BLvIH6b2XEMjW$!iKp%lW'


# Some standard error messages
_MSG_NOT_UNDERSTOOD_PREFIX = f'{HISTORY_INFO_NAME}: Could not understand'
_MSG_EMPTY = 'Field is empty'
_MSG_NO_FLOAT = 'Field is not a floating-point number'
_MSG_NO_STRING = 'Field is not a string'
_MSG_MISSING = 'Mandatory field not found'


# TODO: the _fixup methods may be made even more lenient
# and accept some comment text following each field.
# TODO: this class is a bit messy. Perhaps one could have better code
# with an HistoryInfoEntryField base class (of which all fields are
# subclasses) that collects the base functionality, and subclasses
# for the different fields? Each field could then have a .format and
# a .check method that do what is done here now. It may a bit tricky
# to get it right, though.
# TODO: may be nicer to make _MISSING a MissingField class, then also
# have DuplicateField and ExtraField classes to handle various
# situations. An netry should then probably store its ._raw_contents
# as a list of fields. It would then make way more sense to have
# separate field classes altogether.
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
    r_ref: Optional[float] = None                                               # TODO: how to handle integer/fractional?
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
                    f'{_MSG_NOT_UNDERSTOOD_PREFIX} {tag} '
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
        if self.notes == _MISSING and self.discarded:
            # We fix-up automatically entries with missing
            # 'Notes:' that have been explicitly discarded
            txt += '\n' + self._format_notes()
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
        misses_stuff = any(getattr(self, f.name) == _MISSING
                           for f in fields(self)
                           if f.init and not is_optional_field(f))
        # Notes are marked as optional, but they really are not
        return misses_stuff or self.notes == _MISSING

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
        # Concerning the pylint disable: it would be better indeed
        # to use a named constant, but we'd then have to do it for
        # all. I'd then rather make _TAG an Enum.
        # pylint: disable-next=magic-value-comparison
        fix_time_fmt = 'timestamp' in kwargs
        kwargs.pop('timestamp', None)

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
        """Return an entry object from its string version.

        Parameters
        ----------
        entry_str : str
            The string version of this entry.

        Returns
        -------
        entry: HistoryInfoEntry or PureCommentEntry
            The entry. A PureCommentEntry is returned when entry_str
            does not contain any of the expected fields but only text.
        """
        cls._from_string_check_single_entry(entry_str)
        # Deal with DISCARDED later, as we need notes
        kwargs = {'discarded': False}
        try:
            cls._parse_fields_from_string(entry_str, kwargs)
        except _PureCommentEntryError:
            return PureCommentEntry(entry_str)
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
        set_frozen_attr(new, '_time_format', fmt)
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
            set_frozen_attr(self, attr, '')
            raise EntrySyntaxError(_MSG_EMPTY)

    def _check_missing(self, attr):
        """Complain if no value was provided at all for `attr`."""
        if getattr(self, attr) == _MISSING:
            raise EntrySyntaxError(_MSG_MISSING)

    def _check_notes_field(self):
        """Complain if the Notes field is missing."""
        self._check_missing('notes')

    def _check_r_factor_field(self, which_r):
        """Check whether one of the R-factor fields is fine."""
        value = getattr(self, which_r)
        if value is None:  # R factors are optional
            return
        if isinstance(value, str) and not value.strip():
            set_frozen_attr(self, which_r, '')
            raise EntrySyntaxError(_MSG_EMPTY)
        if isinstance(value, str):  # Non-empty
            value = self._fixup_float_string(value)
            set_frozen_attr(self, which_r, value)
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
        if value is None:  # OK, it's optional
            return
        if not isinstance(self.run_info, str):
            raise EntrySyntaxError(_MSG_NO_STRING)
        value = value.strip()
        if not value:
            set_frozen_attr(self, 'run_info', '')
            raise EntrySyntaxError(_MSG_EMPTY)
        if not re.fullmatch(_SPACE_SEPARATED, value):
            raise EntrySyntaxError('Not a space-separated list of integers.')

    def _check_tensor_nums_field(self):
        """Check whether the tensor_nums attribute is fine."""
        self._check_list_of_int_field('tensor_nums', accept_none=True)

    def _check_list_of_int_field(self, attr, accept_none=False):
        """Check whether the value of `attr` is a valid list of integers."""
        self._check_missing(attr)

        # Below we use set_frozen_attr to circumvent the frozen
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
            set_frozen_attr(self, attr, [])
            return
        if not accept_none and not value:
            set_frozen_attr(self, attr, [])
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
            set_frozen_attr(self, 'timestamp', '')
            raise EntrySyntaxError(_MSG_EMPTY)

        for fmt in TimestampFormat:
            try:
                set_frozen_attr(self,
                                'timestamp',
                                datetime.strptime(value, fmt.value))
            except ValueError:  # Not the right format
                continue
            # Conversion successful. Remember the format,
            # as long as it is acceptable for writing
            if fmt.writable:
                set_frozen_attr(self, '_time_format', fmt)
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
        set_frozen_attr(self, '_formatters', formatters)

    def _fixup_float_string(self, value_str):
        """Try converting `value_str` to a float value."""
        value_str = value_str.strip()
        if not re.fullmatch(_FLOAT_RE, value_str):
            raise EntrySyntaxError(_MSG_NO_FLOAT)
        return ast.literal_eval(value_str)  # Should never fail

    def _fixup_list_of_int_string(self, attr, value_str):
        """Try converting `value_str` to a list of int for `attr`."""
        value_str = value_str.strip()
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
            set_frozen_attr(self, attr, value)
            return
        reason = 'some space-separated items, rather than comma separated'
        self._needs_fixup[attr] = (reason, value)
        raise FixableSyntaxError(reason)

    def _format_field(self, field_name, value_str=None, fmt=''):
        """Return a formatted version of a field."""
        value = getattr(self, field_name)
        if value == _MISSING:
            value_str = ''
        elif value_str is None and value is None:
            return ''
        elif value_str is None and not value:
            value_str = ''
        elif value_str is None:
            value_str = value if isinstance(value, str) else f'{value:{fmt}}'
        return self._format_line(_TAG[field_name], value_str)

    def _format_folder_name(self):
        """Return a line with FOLDER."""
        return self._format_field('folder_name',
                                  value_str=str(self.folder_name).rstrip())

    def _format_job_nums(self):
        """Return a line with JOB ID."""
        value = self.job_nums
        if not isinstance(value, str) and not value:  # Empty
            value = ''
        return self._forma_list_of_ints('job_nums', value)

    @staticmethod
    def _format_line(tag, value):
        """Return a formatted line from a tag and values."""
        return f'{tag:<{_HISTORY_INFO_SPACING}}{value}'

    def _forma_list_of_ints(self, attr, value):
        """Return a line with a list of integers (if the value is correct)."""
        # Important: must check str first, as this is the best way
        # to retain as much as possible the value that users may
        # have removed
        if isinstance(value, str):
            value_str = value
        elif isinstance(value, Sequence):
            # Use list as it does not add trailing comma
            value_str = str(list(value))[1:-1]
        else:  # Not a sequence. Probably some problematic value.
            value_str = repr(value)
        return self._format_field(attr, value_str=value_str)

    def _format_notes(self):
        """Return a string with multi-line notes."""
        notes_str = self.notes.strip()
        if notes_str == _MISSING:
            notes_str = ''
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
        if not isinstance(value, str) and no_tensors:
            value = 'None'
        return self._forma_list_of_ints('tensor_nums', value)

    def _format_timestamp(self):
        """Return a line with TIME."""
        value = self.timestamp
        if isinstance(value, str):
            time_str = value
        elif isinstance(value, datetime):
            time_str = value.strftime(self.time_format.value)
        else:  # Something funny
            time_str = repr(value)
        return self._format_field('timestamp', value_str=time_str)

    @classmethod
    def _from_string_check_single_entry(cls, entry_str):
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
        n_entries = sum(cls._matched_known_fields(m)
                        for m in _ENTRY_RE.finditer(entry_str))
        if n_entries > 1:
            raise EntrySyntaxError(err_msg)

        # Now check the second option. Here we have to take off one
        # line at a time from the notes, and check if the rest matches
        match = _ENTRY_RE.match(entry_str)
        if not match:  # Will be caught in _parse_fields_from_string
            return
        notes = match['notes']
        if notes is None or not notes.strip():  # OK
            return
        if cls._matched_known_fields(_ENTRY_RE.match(notes.strip())):
            raise EntrySyntaxError(err_msg)

    @staticmethod
    def _matched_known_fields(entry_match):
        """Return whether a re.match object identified known entry fields."""
        return any(
            # We only have named groups for known fields
            v for v in entry_match.groupdict().values()
            )

    @classmethod
    def _parse_fields_from_string(cls, entry_str, kwargs):
        """Update `kwargs` by parsing `entry_str`."""
        match = _ENTRY_RE.match(entry_str.strip())
        if not match:  # Safeguard for something wrong with _ENTRY_RE
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
            elif not is_optional_field(field):  # No match at all
                kwargs[name] = _MISSING

        # Now check if there was any extra unknown field.
        known = {v for v in match.groupdict().values() if v}
        extras = [g
                  for g in match.groups()
                  if g and g not in known and g.strip()]
        if extras and not known:
            # Pure comment
            raise _PureCommentEntryError
        if extras:
            raise EntrySyntaxError(
                f'Found unexpected lines in entry\n{entry_str}\n'
                'Problematic lines:\n' + '\n-> '.join(extras)
                )

    @staticmethod
    def _process_string_notes(kwargs):
        """Update `kwargs` according to its notes."""
        notes = kwargs.get('notes', None)
        if notes is None:
            kwargs['notes'] = _MISSING
            return
        kwargs['discarded'] = notes.endswith(_DISCARDED)
        if kwargs['discarded']:
            kwargs['notes'] = notes[:-len(_DISCARDED)].rstrip()


@dataclass(frozen=True)
class PureCommentEntry:
    """An entry with only comments."""

    raw_comment: str = ''

    @property
    def can_be_removed(self):
        """Return that a pure-comment entry can never be removed."""
        return False
