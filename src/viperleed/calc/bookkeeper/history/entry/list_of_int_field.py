"""Module list_of_int_field of viperleed.calc.bookkeeper.history.entry.

Defines the ListOfIntsField abstract base class, its PositiveIntsField,
CommaSeparatedIntsField, and SpaceSeparatedIntsField abstract
subclasses, and the concrete JobIdsField, RunInfoField, and
TensorNumsField classes.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-07-28'
__license__ = 'GPLv3+'

from abc import ABC
from abc import abstractmethod
import ast
from collections.abc import Sequence
import re
from typing import ClassVar
from typing import List
from typing import Tuple
from typing import Union

from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.dataclass_utils import set_frozen_attr
from viperleed.calc.sections.calc_section import CalcSection

from ..errors import EntrySyntaxError
from ..errors import FixableSyntaxError
from .enums import FieldTag
from .field import CommonRegex
from .field import DefaultMessage
from .field import EmptyField
from .field import FieldBase
from .field import MissingField


_INVALID_SECTION_ERR = 'Contains unknown section identifiers'
_NEGATIVE_ERR = 'Contains negative values'
_NOT_SEQUENCE_ERR = 'Not a sequence'
_NOT_INT_ERR = 'Not a list of integers'


def to_comma_separated(value_str):
    """Return a ', '-separated string of items from one with mixed spacing."""
    # Get rid of spaces before commas
    value_str = re.sub(r'\s+,', ',', value_str)
    # Then replace all multi-spaces not preceded by comma
    return re.sub(r'(?<!,)\s+', ',', value_str)


def to_space_separated(value_str):
    """Return a space-separated string of items from one with mixed spacing."""
    commas = to_comma_separated(value_str)
    return commas.replace(',', '')


@frozen
class ListOfIntsField(FieldBase, ABC):
    """A mandatory field with a comma-separated list of ints for value."""

    value: Union[int, Tuple[int], List[int], str] = MissingField
    rgx_loose: ClassVar[str] = None

    def _check_value(self):
        """Check whether this field's value is acceptable."""
        super()._check_value()
        if not self.is_mandatory and self.is_missing:
            return
        if not isinstance(self.value, str):
            set_frozen_attr(self, 'value', tuple(self.value))

    def _check_not_empty(self):
        """Extend the check for an empty field."""
        super()._check_not_empty()
        # Notice that this condition is looser than the one
        # of NoneIsEmptyFields that only checks for 'is None'.
        # This is the reason why we don't inherit that here.
        if not self.value:
            set_frozen_attr(self, 'value', EmptyField)
            raise EntrySyntaxError(DefaultMessage.EMPTY)

    def _check_is_sequence_of_int(self):
        """Raise if self.value is not a sequence of non-negative integers."""
        # Check that it is a sequence
        if not isinstance(self.value, Sequence):
            raise EntrySyntaxError(_NOT_SEQUENCE_ERR)

        # Make sure it is indeed a list of integers
        all_ints = all(
            isinstance(v, int) for v in self.value
            )
        if not all_ints:
            raise EntrySyntaxError(_NOT_INT_ERR)
        self._check_item_values()
        self._store_cleaned_up_string()

    _check_list_value = _check_is_sequence_of_int
    _check_tuple_value = _check_is_sequence_of_int
    _check_sequence_value = _check_is_sequence_of_int

    def _check_int_value(self):
        """Convert a single integer value to a tuple."""
        set_frozen_attr(self, 'value', (self.value,))
        self._check_item_values()
        self._store_cleaned_up_string()

    def _check_item_values(self):  # abstract
        """Raise unless the items of a sequence value are as expected."""

    def _check_str_value(self):
        """Try converting a string value to a tuple of integers."""
        super()._check_str_value()
        set_frozen_attr(self, 'value', self._clean_up_string())
        self._store_cleaned_up_string()
        self._check_item_values()

    def _clean_and_validate_string_loose(self, separator):
        """Raise if a string value does not match self.rgx_loose."""
        regex = self.rgx_loose
        if regex is None:
            raise NotImplementedError('Needs a rgx_loose class attribute!')
        # pylint: disable-next=no-member  # It's a string by now
        value_str = self.value.strip()
        matches = re.fullmatch(regex, value_str)
        if not matches:
            raise EntrySyntaxError(f'Not a {separator}-separated '
                                   'list of integers')
        return value_str

    @abstractmethod
    def _clean_up_string(self):
        """Return a tuple for self.value. Raise if syntax is invalid.

        Subclasses must override this method and NOT modify self.value
        to be a Sequence. The base class will take care of using the
        value returned by this method to update self.value. This method
        should raise EntrySyntaxError of FixableSyntaxError as
        appropriate.

        Returns
        -------
        value_as_tuple : tuple
            The value inferred from self.value (a string by the time
            this method is called). Will replace self.value after
            successful execution of this method.

        Raises
        ------
        EntrySyntaxError
            If self.value (a string) has an unacceptable format that
            made it not understandable.
        FixableSyntaxError
            If self.value (a string) has an unacceptable format that
            could be fixed.
        """

    def _store_cleaned_up_string(self):
        """Store a value_str into self._value_str after cleaning up."""
        # Use list since, compared to tuple, as it does
        # not add trailing comma for single elements
        value_str = str(list(self.value))[1:-1]
        set_frozen_attr(self, '_value_str', value_str)


@frozen
class CommaSeparatedIntsField(ListOfIntsField):
    """A comma-separated list of integers."""

    rgx_loose: ClassVar[str] = CommonRegex.COMMA_OR_SPACE_SEPARATED_INTS.value

    def _clean_up_string(self):
        """Raise if self.value is not a comma-separated list of integers."""
        value_str = self._clean_and_validate_string_loose('comma')
        comma_separated = re.fullmatch(CommonRegex.COMMA_SEPARATED_INTS.value,
                                       value_str)
        if not comma_separated:
            value_str = to_comma_separated(value_str)

        value = tuple(ast.literal_eval(value_str + ','))
        if comma_separated:
            return value
        reason = 'some space-separated items, rather than comma separated'
        raise FixableSyntaxError(reason=reason, fixed_value=value)


@frozen
class PositiveIntsField(ListOfIntsField):
    """A list of non-negative integers."""

    def _check_item_values(self):
        """Raise if some items are negative."""
        if any(v < 0 for v in self.value):
            raise EntrySyntaxError(_NEGATIVE_ERR)


@frozen
class SpaceSeparatedIntsField(ListOfIntsField):
    """A space-separated list of integers."""

    rgx_loose: ClassVar[str] = CommonRegex.COMMA_OR_SPACE_SEPARATED_INTS.value

    def _clean_up_string(self):
        """Raise if self.value is not a comma-separated list of integers."""
        value_str = self._clean_and_validate_string_loose('space')

        space_separated = re.fullmatch(CommonRegex.SPACE_SEPARATED_INTS.value,
                                       value_str)
        comma_separated = to_comma_separated(value_str)
        value = tuple(ast.literal_eval(comma_separated + ','))
        if space_separated:
            return value
        reason = 'some comma-separated items, rather than space separated'
        raise FixableSyntaxError(reason=reason, fixed_value=value)

    def _store_cleaned_up_string(self):
        """Store a value_str into self._value_str after cleaning up."""
        super()._store_cleaned_up_string()  # Gives commas
        space_separated = to_space_separated(self._value_str)
        set_frozen_attr(self, '_value_str', space_separated)


@frozen
class JobIdsField(PositiveIntsField, CommaSeparatedIntsField,
                  tag=FieldTag.JOB_NUMS, mandatory=True):
    """The field of the job numbers executed during this run."""


@frozen
class RunInfoField(SpaceSeparatedIntsField, tag=FieldTag.RUN_INFO):
    """An optional field containing information about run segments."""

    def _check_item_values(self):
        """Raise if some items are not section numbers."""
        section_ids = {s.value for s in CalcSection}
        if any(v not in section_ids for v in self.value):
            raise EntrySyntaxError(_INVALID_SECTION_ERR)


@frozen
class TensorNumsField(PositiveIntsField, CommaSeparatedIntsField,
                      tag=FieldTag.TENSOR_NUMS, mandatory=True):
    """The field of the tensor numbers used or generated during this run."""

    value: Union[int, Tuple[int], List[int], str, None] = MissingField

    @property
    def no_tensors(self):
        """Return whether no tensors were used/generated."""
        # About the disable: this seems to be a pylint inference bug.
        # It looks like it cannot realize that value.strip() is behind
        # a type check for string.
        # pylint: disable-next=no-member
        return isinstance(self.value, str) and self.value.strip() == str(None)

    def _check_str_value(self):
        """Accept 'None' as a valid string value."""
        if self.no_tensors:
            self._store_none()
            return
        super()._check_str_value()

    def _check_not_empty(self):
        """Extend the check for an empty field to include None-likes."""
        value = self.value
        is_akin_to_none = (  # Variations of None
            value in (None, [None], (None,))
            # pylint: disable-next=no-member  # pylint inference bug
            or (isinstance(value, str) and value.strip() == str(None))
            )
        is_initialization = value in (CalcSection.INITIALIZATION.value,
                                      (CalcSection.INITIALIZATION.value,),
                                      [CalcSection.INITIALIZATION.value])
        if is_akin_to_none or is_initialization:
            self._store_none()
            return
        # Notice that it is important to use super at the end, as
        # it checks for emptiness by using the logical value of
        # self.value, which would be False if the value is None
        # (that is instead acceptable).
        super()._check_not_empty()

    def _store_none(self):
        """Store a non-False value of 'None'."""
        set_frozen_attr(self, 'value', str(None))
        set_frozen_attr(self, '_value_str', self.value)
