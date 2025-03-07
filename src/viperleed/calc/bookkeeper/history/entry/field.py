"""Module field of viperleed.calc.bookkeeper.history.entry.

Defines base classes for a single 'line' of a history entry.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-07-28'
__license__ = 'GPLv3+'

from collections import namedtuple
from contextlib import contextmanager
from dataclasses import fields as data_fields
from enum import Enum
import re
from typing import Any
from typing import ClassVar
from typing import Dict
from typing import Tuple
from typing import Union

from viperleed.calc.lib.dataclass_utils import check_types
from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import non_init_field
from viperleed.calc.lib.dataclass_utils import replace_values
from viperleed.calc.lib.dataclass_utils import set_frozen_attr
from viperleed.calc.lib.string_utils import to_snake_case
from viperleed.calc.lib.string_utils import strip_comments

from ..errors import EntrySyntaxError
from ..errors import FixableSyntaxError
from ..errors import HistoryInfoError
from .enums import FaultyLabel
from .enums import FieldTag


class CommonRegex(Enum):
    """Regular-expression patterns (and portions thereof) used for fields."""

    ANY = r'.*'
    SPACE = r'[ \t]'  # Exclude \n
    COMMA_SEPARATED_INTS = rf'\d+({SPACE}*,{SPACE}*\d+)*'
    COMMA_OR_SPACE_SEPARATED_INTS = rf'\d+({SPACE}*[,{SPACE}]?{SPACE}*\d+)*'
    FLOAT = r'\d+(.\d+)?'
    MULTILINE = r'(?:(?s:.)*)'  # Same as re.DOTALL
    ONE_LINE = ANY
    SPACE_SEPARATED_INTS = rf'\d+({SPACE}+\d+)*'


class DefaultMessage(Enum):
    """Some default error messages."""

    EMPTY = 'Field is empty'
    NOT_FLOAT = 'Field is not a floating-point number'
    NOT_STRING = 'Field is not a string'
    MISSING = 'Mandatory field not found'


class MissingField:  # pylint: disable=too-few-public-methods
    """A sentinel for a field that is completely absent."""


class EmptyField:    # pylint: disable=too-few-public-methods
    """A sentinel for a field that has no value."""


FixedFieldValue = namedtuple('FixedFieldValue', ('reason', 'value'))


@frozen
class FieldBase:
    """Base class for all the fields in a history.info entry."""

    # str is typically for not-understood values.
    # Subclasses define concrete supported types.
    value: Union[Any, str] = MissingField

    # _value_str can be set in subclasses to a string value
    # that is used instead of self.value for __str__.
    _value_str: str = non_init_field()
    _comment: str = non_init_field(default='')
    _needs_fix: FixedFieldValue = non_init_field()
    _not_understood: str = non_init_field()

    # Class attributes set in __init_subclass__
    is_mandatory: ClassVar[bool] = False
    tag: ClassVar[FieldTag] = None

    # Pattern used for string matching in .from_string:
    rgx_pattern: ClassVar[str] = CommonRegex.ONE_LINE.value

    # A cache of know sub-classes:
    _subclasses: ClassVar[Dict[FieldTag,'FieldBase']] = {}

    # Non-init attributes to be skipped when this field is being fixed:
    _skip_when_fixing: ClassVar[set] = {
        '_needs_fix',  # We're fixing the instance anyway
        '_value_str',  # Assume it will be recomputed for a new one
        }

    def __init_subclass__(cls, tag=None, mandatory=False):
        """Register subclasses with concrete tags."""
        already_handled = cls._subclasses.get(tag, None)
        if already_handled:
            raise ValueError('Cannot have two classes for the same tag. '
                             f'{already_handled} already takes care of {tag}')
        if tag is not None:
            cls._subclasses[tag] = cls
        # At this point, the class is not frozen yet.
        setattr(cls, 'tag', tag)
        setattr(cls, 'is_mandatory', bool(mandatory))

    # pylint: disable-next=unused-argument   # pylint bug #9843
    def __new__(cls, *args, **kwargs):
        """Complain when making an instance of a tag-less class."""
        if cls.tag is None:
            raise TypeError('Cannot instantiate a tag-less entry field.')
        return super().__new__(cls)

    def __post_init__(self):
        """Check initialization values for this field."""
        if isinstance(self.value, FieldBase):
            self._update_from_field(self.value)
            return
        if self.value is MissingField:  # Errors will pop up later
            return
        try:
            check_types(self, init_only=True)
        except TypeError as exc:
            set_frozen_attr(self, '_not_understood', str(exc))

    def __str__(self):
        """Return a string version of this field."""
        if self.is_missing:
            return ''
        return self._format_string_value(self._get_string_value())

    @property
    def is_empty(self):
        """Return whether this field has an empty value."""
        return self.value is EmptyField

    @property
    def is_missing(self):
        """Return whether this field is not defined at at all."""
        return self.value is MissingField

    @property
    def has_comments(self):
        """Return whether any comment was found at the right."""
        return bool(self._comment)

    @property
    def needs_fixing(self):
        """Return whether this field is wrongly formatted."""
        return self._needs_fix is not None

    @property
    def was_understood(self):
        """Return whether this field was recognized."""
        return self._not_understood is None

    def as_fixed(self):
        """Return a version of this field with its fixable value replaced."""
        if not self.needs_fixing:
            return self
        # About the disable: pylint can't know that _needs_fix
        # is set to a FixedFieldValue instance in this case
        # pylint: disable-next=no-member
        fixed_value = self._needs_fix.value
        if fixed_value is None:
            raise HistoryInfoError('No value for fixing field')
        assert self.was_understood
        return replace_values(self, value=fixed_value,
                              skip=self._skip_when_fixing)

    def check_value(self):
        """Raise if the value given at initialization is funny.

        Subclasses must not override this method, but rather the
        private implementation (_check_value) instead. Another
        option is implementing methods named
        _check_<snake_case_type_name>_value if they want to handle
        a value whose type is named SnakeCaseTypeName. These methods
        are called appropriately, based on the type of self.value.
        FieldBase implements a basic _check_str_value.

        Raises
        ------
        FixableSyntaxError
            If the value is understandable, but it is somewhat funny
            (e.g., invalid format).
        EntrySyntaxError
            If the value is not acceptable.
        """
        with self._register_errors():
            self._check_value()

    @classmethod
    def for_tag(cls, tag):
        """Return the subclass of FieldBase that handles tag."""
        try:
            return cls._subclasses[tag]
        except KeyError:
            raise  ValueError(f'No subclass to handle {tag}') from None

    def format_faulty(self, with_label=None):
        """Return a line indicating the problems in this field."""
        if self.is_missing and not self.is_mandatory:
            return ''
        prefix = with_label or FaultyLabel.for_field(self)
        contents = (self if not self.is_missing
                    else self._format_string_value('???'))
        return f'{prefix}{contents}'

    @classmethod
    def from_string(cls, field_string):
        """Return an instance of this field from its string contents.

        Parameters
        ----------
        field_string : str
            The string version of one or more lines of a
            history.info entry.

        Returns
        -------
        field : FieldBase
            A concrete FieldBase subclass that can understand
            `field_string`.

        Raises
        ------
        HistoryInfoError
            If there is no known subclass to interpret `field_string`.
        """
        # Figure out which subclass understands field_string.
        for tag in FieldTag:
            try:
                subclass = cls._subclasses[tag]
            except KeyError:
                continue
            pattern = (
                rf'{tag.value}'
                rf'(?P<field_value>{CommonRegex.SPACE.value}*'
                rf'{subclass.rgx_pattern})'
                )
            match_ = re.fullmatch(pattern, field_string)
            if match_ and match_['field_value'] is not None:
                break
        else:
            raise HistoryInfoError(f'No field can parse {field_string!r}')
        # Remove the tag and its trailing spaces. This way we
        # catch only the real contents of the field, and prevent
        # accumulation of leading white spaces that would end up
        # in match_['field_value'].
        raw_value = field_string.replace(str(tag), '', 1)
        if raw_value == field_string and raw_value.startswith(tag.value):
            # Probably the white space is missing
            raw_value = raw_value.replace(tag.value, '', 1)
        # pylint: disable-next=protected-access  # It's a subclass
        value_str, comments = subclass._strip_comments(raw_value)
        instance = subclass(value=value_str)
        set_frozen_attr(instance, '_comment', comments)
        return instance

    @contextmanager
    def _register_errors(self):
        """Register EntrySyntaxError and FixableSyntaxError when executing."""
        try:
            yield
        except FixableSyntaxError as exc:
            set_frozen_attr(self, '_needs_fix',
                            FixedFieldValue(exc.reason, exc.fixed_value))
            raise
        except EntrySyntaxError as exc:
            set_frozen_attr(self, '_not_understood', str(exc))
            raise

    @staticmethod
    def _strip_comments(string_value):
        """Return non-commented and commented parts of `string_value`."""
        without_comments = strip_comments(string_value,
                                          strip_whitespaces=False)
        comments = string_value.replace(without_comments, '', 1)
        return without_comments, comments

    def _check_not_missing(self):
        """Complain if no value was provided at all for `attr`."""
        if self.is_missing:
            raise EntrySyntaxError(DefaultMessage.MISSING)

    def _check_not_empty(self):
        """Complain if this field is empty."""
        value = self.value

        # About the disable: this seems to be a pylint inference bug.
        # It looks like it cannot realize that value.strip() is behind
        # a type check for string.
        # pylint: disable-next=no-member
        if isinstance(value, str) and not value.strip():
            set_frozen_attr(self, 'value', EmptyField)
        if self.is_empty:
            raise EntrySyntaxError(DefaultMessage.EMPTY)

    def _check_str_value(self):
        """Raise EntrySyntaxError if self.value is not a string."""
        if not isinstance(self.value, str):
            raise EntrySyntaxError(DefaultMessage.NOT_STRING)
        if not re.fullmatch(self.rgx_pattern, self.value):
            raise EntrySyntaxError('Unexpected format')
        self._cleanup_str_value()

    def _check_value(self):
        """Raise if the value given at initialization is funny.

        The base-class implementation complains if a field is empty
        or if a mandatory field is missing. The value of mandatory
        or non-missing fields is then checked according to its type,
        as long as a suitable checker method is found. Subclasses
        may extend this method to perform more checks.

        Raises
        ------
        EntrySyntaxError
            If a field is empty or a mandatory one is missing.
        EntrySyntaxError
            If a field value has an unexpected type (i.e., one
            for which no checker is found).
        EntrySyntaxError
            If a field value is not understood.
        FixableSyntaxError
            If a field value is understood, but does not follow
            the expected format.
        """
        if self.is_mandatory:
            self._check_not_missing()
        self._check_not_empty()
        if not self.is_mandatory and self.is_missing:
            return
        if self.is_empty:  # Complained already, or empty is acceptable
            return
        value_type = to_snake_case(type(self.value).__name__)
        try:
            checker = getattr(self, f'_check_{value_type}_value')
        except AttributeError:  # Some type we don't accept
            raise EntrySyntaxError(DefaultMessage.NOT_STRING) from None
        checker()

    def _cleanup_str_value(self):
        """Store a cleaned-up version of a string value."""
        # pylint: disable-next=no-member   # Pylint inference issue
        set_frozen_attr(self, 'value', self.value.rstrip())

    def _format_string_value(self, value_str):
        """Return a formatted and tagged version of `value_str`."""
        tag_str = str(self.tag)
        if not value_str.strip():
            tag_str = tag_str.rstrip()
        return f'{tag_str}{value_str}'

    def _get_string_value(self):
        """Return a string value from self."""
        if self.is_empty:
            return ''
        value_str = self._value_str
        if value_str is None and not self.value:
            return self._comment
        if value_str is None:
            return f'{self.value}  {self._comment}'.rstrip()
        return f'{value_str}  {self._comment}'.rstrip()

    def _update_from_field(self, other):
        """Return this field filled with the values of another one."""
        if not isinstance(other, type(self)):
            raise TypeError(
                f'{type(self).__name__}: Cannot fetch attributes from '
                f'incompatible field of type {type(other).__name__!r}'
                )
        for field in data_fields(self):
            attr = field.name
            # No need to worry about ClassVar, as it does not show up
            set_frozen_attr(self, attr, getattr(other, attr))


@frozen
class CommentLessField(FieldBase):
    """A field that stores its comments in the value."""

    @staticmethod
    def _strip_comments(string_value):
        """Return non-commented and commented parts of `string_value`."""
        # Include any comment in the overall value of this field
        return string_value, ''


@frozen
class MultiLineField(CommentLessField):
    """A field that spans multiple lines."""

    # Lines are stored as a tuple of strings. May also be
    # a multi-line string before .check_value is called.
    value: Union[Tuple[str], str, Any] = MissingField

    rgx_pattern: ClassVar[str] = CommonRegex.MULTILINE.value

    def format_faulty(self, with_label=None):
        """Return lines indicating the problems in this field."""
        single_string = super().format_faulty(with_label=with_label)
        lines = single_string.splitlines(True)  # Keeps \n
        indent = str(FaultyLabel.OK)
        indented_lines = [lines[0]]
        indented_lines.extend(indent + line for line in lines[1:])
        return ''.join(indented_lines)

    def _check_tuple_value(self):
        """Make sure all items are strings."""
        if not self.value:
            set_frozen_attr(self, 'value', EmptyField)
            return
        if isinstance(self.value, str):
            raise TypeError(f'{type(self).__name__}: Found unprocessed '
                            'string value. Call check_value() or '
                            '_cleanup_str_value() before')
        if not isinstance(self.value, tuple):
            raise EntrySyntaxError(
                f'{type(self).__name__} only accepts \'str\' or \'tuple\'. '
                f'Found {type(self.value).__name__!r} instead.'
                )
        # pylint: disable-next=not-an-iterable      # Can't infer tuple
        if not all(isinstance(i, str) for i in self.value):
            raise EntrySyntaxError(DefaultMessage.NOT_STRING)

    def _cleanup_str_value(self):
        """Store a string value as a tuple of lines without trailing spaces."""
        set_frozen_attr(self, 'value', self._convert_str_to_tuple(self.value))

    @staticmethod
    def _convert_str_to_tuple(value_str):
        """Return a tuple of lines without trailing spaces from `value_str`."""
        # This preserves the number of '\n' characters
        lines = [line.rstrip() for line in value_str.splitlines()]
        if value_str.endswith('\n'):  # Don't loose a final newline
            lines.append('')
        return tuple(lines)

    def _get_string_value(self):
        """Return a string version of the lines of this MultiLineField."""
        lines = self._prepare_lines_for_str()
        if lines is None:
            return super()._get_string_value()
        return '\n'.join(lines)

    def _prepare_lines_for_str(self):
        """Return a tuple of lines ready to '\n'-join. None otherwise."""
        if not self.is_empty and not self.is_missing:
            try:  # pylint: disable=too-many-try-statements
                # Here we make sure to realize whether there are issues
                # with the value. This may mean that .check_value was
                # not called appropriately before attempting to make a
                # string. We cannot place the "with" statement outside
                # the try...except block, as we would miss the
                # EntrySyntaxError that we want to register.
                with self._register_errors():
                    self._check_tuple_value()
            except (TypeError, EntrySyntaxError):
                # Something funny with the value
                return None
        # Notice that it is important to check emptiness now:
        # _check_tuple_value above may set the value to empty.
        if self.is_empty or self.is_missing:
            return None
        return self.value


class NoneIsEmptyField(FieldBase):
    """A field that considers a None value as an empty condition."""

    def _check_not_empty(self):
        """Complain if this field is empty."""
        super()._check_not_empty()
        if self.value is None:
            set_frozen_attr(self, 'value', EmptyField)
            raise EntrySyntaxError(DefaultMessage.EMPTY)


@frozen
class UnknownField(CommentLessField, tag=FieldTag.UNKNOWN):
    """A dummy field for lines that have no tag."""
