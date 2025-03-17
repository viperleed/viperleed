"""Module time_field of viperleed.calc.bookkeeper.history.entry.

Defines the TimestampFormat and TimestampField classes for handling a
line of a history.info entry with the time stamp of the run.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-07-28'
__license__ = 'GPLv3+'

from datetime import datetime
from enum import Enum
from typing import Optional
from typing import Union

from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import replace_values
from viperleed.calc.lib.dataclass_utils import set_frozen_attr
from viperleed.calc.lib.time_utils import DateTimeFormat

from ..errors import EntrySyntaxError
from ..errors import FixableSyntaxError
from .enums import FieldTag
from .field import DefaultMessage
from .field import FieldBase
from .field import FixedFieldValue
from .field import MissingField

_OUTDATED = 'outdated format'


class TimestampFormat(Enum):
    """A collection of known formats for timestamps."""

    # Underscore names are used only for conversion purposes
    # while parsing strings, and never for formatting.
    GERMAN = '%d.%m.%y %H:%M:%S'              # < v0.13.0
    ISO = DateTimeFormat.ISO.value            # >= v0.13.0
    _CALC = DateTimeFormat.FILE_SUFFIX.value  # Format of log names
    DEFAULT = ISO  # Just an alias of one of the above

    @property
    def writable(self):
        """Return whether self is suitable for writing to history.info."""
        return not self.name.startswith('_')


@frozen
class TimestampField(FieldBase, tag=FieldTag.TIMESTAMP, mandatory=True):
    """A field containing a time stamp for the history run."""

    value: Union[datetime, str] = MissingField
    time_format: Optional[TimestampFormat] = None

    def as_fixed(self):
        """Return a version of this field with the default format."""
        # Overriding here is needed, as only the format needs to change
        return self.with_format(TimestampFormat.DEFAULT)

    def check_value(self):
        """Complain if a value is inappropriate."""
        error = None
        try:
            super().check_value()
        except (EntrySyntaxError, FixableSyntaxError) as exc:
            error = exc
        with self._register_errors():
            self._check_format_consistent()
        if error:
            raise error

    def with_format(self, fmt):
        """Return a version of this field with another date-time format."""
        if self.is_missing:
            raise EntrySyntaxError(DefaultMessage.MISSING)
        if not self.was_understood:
            raise EntrySyntaxError(self._not_understood)

        fmt = self._validate_time_format(fmt)
        if fmt is self.time_format:
            return self

        # Skip _value_str to recompute it the first time it is needed.
        instance = replace_values(self, time_format=fmt, skip='_value_str')
        fix_reason = (
            None if fmt is TimestampFormat.DEFAULT
            else FixedFieldValue(_OUTDATED, None)
            )
        set_frozen_attr(instance, '_needs_fix', fix_reason)
        return instance

    def _check_datetime_value(self):
        """Raise if this field is a datetime, but it has no format."""
        if self.time_format is None:
            raise EntrySyntaxError('Missing time-stamp format')

    def _check_format_consistent(self):
        """Raise if a time_format was given for an invalid value."""
        if self.was_understood or not self.time_format:
            return
        if self.is_empty:
            reason = DefaultMessage.EMPTY.value
        elif self.is_missing:
            reason = DefaultMessage.MISSING.value
        else:
            reason = 'Field has an invalid value'
        set_frozen_attr(self, 'time_format', None)
        set_frozen_attr(self, '_needs_fix', None)  # Can't fix this
        raise EntrySyntaxError(f'Cannot have a time_format defined: {reason}')

    def _check_str_value(self):
        """Raise if self.value is not a valid timestamp."""
        self._sanitize_string_value()
        super()._check_str_value()

        # See if we can parse the string with any known format
        value_str = self.value
        for fmt in TimestampFormat:
            if self._successfully_parsed(value_str, fmt):
                return

        # No valid format found. Make sure we're not storing any.
        set_frozen_attr(self, 'time_format', None)
        raise EntrySyntaxError('Field could not be parsed as a date-time')

    def _get_string_value(self):
        """Store and return a string value for this TimestampField."""
        if isinstance(self.value, datetime) and self.time_format:
            as_string = self.value.strftime(self.time_format.value)
            set_frozen_attr(self, '_value_str', as_string)
        return super()._get_string_value()

    def _sanitize_string_value(self):
        """Clean up a string value."""
        try:
            value_str = self.value.replace('moved-', '').strip()
        except AttributeError:
            # Not a string. Error pops up in super()._check_str_value()
            return
        # Replace self.value for the next check
        set_frozen_attr(self, 'value', value_str)

    def _successfully_parsed(self, value_str, fmt):
        """Attempt interpreting `value_str` with a date-time `fmt`."""
        try:
            parsed = datetime.strptime(value_str, fmt.value)
        except ValueError:  # Not the right format
            return False

        # Conversion successful. Remember value and format.
        set_frozen_attr(self, 'value', parsed)
        set_frozen_attr(self, 'time_format', fmt)
        if fmt is not TimestampFormat.DEFAULT:
            raise FixableSyntaxError(reason=_OUTDATED, fixed_value=None)
        return True

    @staticmethod
    def _validate_time_format(fmt):
        """Return a valid TimestampFormat from `fmt`."""
        known_fmts = ', '.join(repr(f.name) for f in TimestampFormat
                               if f.writable)
        if isinstance(fmt, str):
            try:
                return TimestampFormat[fmt.upper()]
            except KeyError:
                raise ValueError(f'Unknown timestamp format {fmt!r}. '
                                 f'Should be one of {known_fmts}') from None
        if not isinstance(fmt, TimestampFormat):
            raise TypeError(f'Expected {known_fmts}, or TimestampFormat. '
                            f'Found {type(fmt).__name__!r}')
        return fmt
