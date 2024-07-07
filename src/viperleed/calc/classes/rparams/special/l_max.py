"""Module l_max of viperleed.calc.classes.rparams.special.

Defines the LMax class, a convenience container for parameter LMAX.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-12-16'
__license__ = 'GPLv3+'

from dataclasses import dataclass

from ..defaults import NO_VALUE
from ..limits import PARAM_LIMITS
from .base import SpecialParameter


_MIN, _MAX = PARAM_LIMITS['LMAX']


@dataclass(repr=False)
class LMax(SpecialParameter, param='LMAX'):
    """A container for user-supplied LMAX values."""

    _min : int = NO_VALUE
    _max : int = NO_VALUE

    def __post_init__(self):
        """Process the input values to integers."""
        for attr_name in ('_min', '_max'):
            attr = getattr(self, attr_name)
            value = self._check_int_in_range(attr)
            setattr(self, attr_name, value)
        self._sort_values()

    def __bool__(self):
        """Return whether any of the LMax values is defined."""
        return self.has_max or self.has_min

    def __repr__(self):
        """Return a representation string of this LMax."""
        txt = 'LMax('
        if not self:
            return txt + ')'
        txt += str(self.min)  # == .max if only one value
        if not self.has_single_value:
            txt += f', {self.max}'
        return txt + ')'

    def __str__(self):
        """Return a string version of this LMax."""
        if not self:
            return ''
        if self.has_single_value:
            return str(self.max)
        return f'{self.min}-{self.max}'

    @staticmethod
    def _check_int_in_range(value):
        """Return an in-range integer version of value. Raise otherwise."""
        if value is NO_VALUE:
            return value
        try:
            f_value = float(value)
        except (TypeError, ValueError):
            raise TypeError('LMax must be numeric') from None
        i_value = int(f_value)
        if abs(i_value - f_value) > 1e-6:
            raise TypeError('LMax must be integer')
        if not _MIN <= i_value <= _MAX:
            raise ValueError(f'Invalid value={value} for LMax. Values '
                             f'must be between {_MIN} and {_MAX}')
        return i_value

    def _sort_values(self):
        """Swap min and max so that they're sorted as expected."""
        if not self or self.has_single_value:
            return
        if self._min > self._max:
            self._min, self._max = self._max, self._min

    @property
    def has_min(self):
        """Return whether there is a defined minimum value."""
        return self._min is not NO_VALUE

    @property
    def has_max(self):
        """Return whether there is a defined maximum value."""
        return self._max is not NO_VALUE

    @property
    def has_single_value(self):
        """Return whether this LMax does not correspond to a range."""
        return bool(self) and self.min == self.max

    @property
    def max(self):
        """Return the largest value of LMAX to be used for calculations."""
        return self._max if self.has_max else self._min

    @max.setter
    def max(self, new_max):
        """Assign a new maximum value for this LMax."""
        value = self._check_int_in_range(new_max)
        self._max = value
        self._sort_values()

    @property
    def min(self):
        """Return the smallest value of LMAX to be used for calculations."""
        return self._min if self.has_min else self._max

    @min.setter
    def min(self, new_min):
        """Assign a new minimum value of this LMax."""
        value = self._check_int_in_range(new_min)
        self._min = value
        self._sort_values()

    @classmethod
    def from_value(cls, value):
        """Return an LMax from a 2-tuple value."""
        return cls(*value)
