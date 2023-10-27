"""Module energy_range of viperleed.tleedmlib.classes.rparams.special.

Created on 2023-10-23

@author: Michele Riva (@michele-riva)

Defines the EnergyRange class, used as a base for parameters that
require a user input of the form 'start stop step'.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from math import isfinite
from numbers import Real

from ._base import SpecialParameter
from .._defaults import NO_VALUE


# Notice that we purposely do not pass a param keyword here, as
# this is exclusively a base class, and it is not associated with
# any parameter. Hence, it is not registered by SpecialParameter.

@dataclass
class EnergyRange(SpecialParameter):
    """A container for energies."""

    start: float = NO_VALUE
    stop: float = NO_VALUE
    step: float = NO_VALUE

    def __eq__(self, other):
        """Return whether this EnergyRange is identical to another."""
        if not isinstance(other, (EnergyRange, Sequence)):
            return NotImplemented
        equal = tuple(self) == tuple(other)
        return equal or NotImplemented

    def __iter__(self):
        """Yield items of self."""
        yield self.start
        yield self.stop
        yield self.step

    def __post_init__(self):
        """Check and process initialization values."""
        non_defaults = self._non_defaults
        if not all(isinstance(e, Real) for e in non_defaults):
            raise TypeError('Values must be real')
        if not all(isfinite(e) for e in non_defaults):
            raise ValueError('All values must be finite')
        self.check_consistency()

    @property
    def min(self):
        """Return the lower limit of this EnergyRange."""
        return self.start

    @property
    def max(self):
        """Return the upper limit of this EnergyRange."""
        return self.stop

    @property
    def _non_defaults(self):
        """Return the non-default values in self."""
        return [e for e in self if e is not NO_VALUE]

    @classmethod
    def from_value(cls, value):
        """Return an EnergyRange from value."""
        return cls(*value)

    def check_consistency(self):
        """Change inconsistent values or complain."""
        try:
            1 / self.step
        except ZeroDivisionError:
            raise ValueError('Step cannot be zero') from None
        except TypeError:  # NO_VALUE
            pass

        start, stop, step = self
        has_start_stop = all(v is not NO_VALUE for v in (start, stop))
        has_values = has_start_stop and step is not NO_VALUE

        if has_values and (stop - start) * step < 0:
            raise ValueError('Inconsistent step. Cannot shift from '
                             f'{start:.2f} to {stop:.2f} with {step=:.2f}')
        if has_start_stop and self._swap and stop < start:
            self._swap()

    def _swap(self):
        """Swap start and stop, change step."""
        self.start, self.stop = self.stop, self.start
        if self.step is not NO_VALUE:
            self.step *= -1
