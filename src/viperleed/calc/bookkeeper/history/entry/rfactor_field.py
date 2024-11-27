"""Module rfactor_field of viperleed.calc.bookkeeper.history.entry.

Defines the RFactorField class and its subclasses RRefField and
and RSuperField.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-28'
__license__ = 'GPLv3+'


import re
from typing import Union

from viperleed.calc.lib.dataclass_utils import set_frozen_attr
from viperleed.calc.lib.dataclass_utils import frozen

from ..errors import EntrySyntaxError
from .enums import FieldTag
from .field import CommonRegex
from .field import DefaultMessage
from .field import NoneIsEmptyField
from .field import MissingField

_EPS = 1e-5   # Tolerance on the limits for the R factor
_FMT = '.4f'  # Format for single-float value
_FLOAT_RE = CommonRegex.FLOAT.value
_R_FACTOR_RE = (  # Regular expression for parsing string value                 # TODO: consider whether we want to be more strict about spaces. May mean was edited.
    rf'(?P<r_tot>{_FLOAT_RE})\s*'
    rf'(?:\(\s*(?P<r_int>{_FLOAT_RE})\s*'
    rf'/\s*(?P<r_frac>{_FLOAT_RE})\s*\)\s*)?'
    )


@frozen
class RFactorField(NoneIsEmptyField):
    """Base class for all R-factor-containing fields."""

    # Either a string of the form "R_total (R_integer / R_fractional)"
    # or a single floating-point value. When a string does not contain
    # the part in parentheses, it is converted to a float.
    value: Union[str, float] = MissingField

    def _check_float_value(self):
        """Check that a float value is acceptable."""
        if not isinstance(self.value, float):
            raise EntrySyntaxError(DefaultMessage.NOT_FLOAT)
        self._check_in_range(self.value)
        set_frozen_attr(self, '_value_str', f'{self.value:{_FMT}}')

    @staticmethod
    def _check_in_range(value):
        """Raise if `value` is not an acceptable R factor float value."""
        in_range = -_EPS < value < 2 + _EPS
        if not in_range:
            raise EntrySyntaxError(f'{value} is not between zero and two')

    def _check_str_value(self):
        """Check that a string value is acceptable."""
        super()._check_str_value()
        # Store a left-stripped version of the value to
        # ensure we don't accumulate spaces at the left
        # pylint: disable-next=no-member  # Don't handle AttributeError
        set_frozen_attr(self, 'value', self.value.lstrip())
        value_str = self.value.strip()    # pylint: disable=no-member
        match = re.fullmatch(_R_FACTOR_RE, value_str)
        if not match:
            raise EntrySyntaxError(DefaultMessage.NOT_FLOAT)

        # Now check that the various bits make sense
        for r_value in match.groupdict().values():
            try:
                self._check_in_range(float(r_value))
            except TypeError:  # Did not match
                pass

        if not match['r_int']:  # Single float, store as such
            assert not match['r_frac']
            # The next line may potentially raise TypeError if
            # r_tot did not match (and is thus None). However,
            # this should never happen if there is a match.
            set_frozen_attr(self, 'value', float(match['r_tot']))

        set_frozen_attr(self, '_value_str', value_str)


@frozen
class RRefField(RFactorField, tag=FieldTag.R_REF):
    """The field of the R factor from a reference calculation."""


@frozen
class RSuperField(RFactorField, tag=FieldTag.R_SUPER):
    """The field of the R factor from a search/superpos calculation."""
