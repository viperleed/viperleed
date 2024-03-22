"""Module energy_range of viperleed.calc.classes.rparams.special.

Classes
-------
EnergyRange
    Base class for parameters with 'start stop step' user input .
TheoEnergies(EnergyRange)
    Used as the Rparams attribute THEO_ENERGIES. The start -- stop
    intervalB is an integer multiple of step. All attributes are
    strictly positive.
IVShiftRange
    Used as the Rparams attribute IV_SHIFT_RANGE. Bounds are integer
    multiples of step. The step attribute is strictly positive.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-23'
__license__ = 'GPLv3+'

import ast
from collections.abc import Sequence
from dataclasses import dataclass
from decimal import Decimal
from math import ceil, floor, isfinite, remainder
from numbers import Real

from ._base import SpecialParameter
from .._defaults import NO_VALUE


EPS = 1e-8  # Tolerance for comparisons of floats

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
        tuple_other = tuple(other)
        if not self.has_step and len(tuple_other) == 2:
            tuple_other = *tuple_other, NO_VALUE
        if len(tuple(self)) != len(tuple_other):
            return NotImplemented
        scale = abs(self.step) if self.has_step else 1
        for v_self, v_other in zip(self, tuple_other):
            if (v_self, v_other) == (NO_VALUE, NO_VALUE):
                continue
            if v_self is NO_VALUE or v_other is NO_VALUE:
                return NotImplemented
            if abs(v_self - v_other) > EPS*scale:
                return NotImplemented
        return True

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
        self._check_consistency()

        # Mess with start/stop only if all the values are present
        if self.defined:
            self.adjust_to_fit_step()

    @property
    def min(self):
        """Return the lower limit of this EnergyRange."""
        return self.start

    @property
    def max(self):
        """Return the upper limit of this EnergyRange."""
        return self.stop

    @property
    def defined(self):
        """Return whether all of the attributes have a value set."""
        return not any(v is NO_VALUE for v in self)

    @property
    def has_bounds(self):
        """Return whether both start and stop are defined."""
        return not any(v is NO_VALUE for v in (self.start, self.stop))

    @property
    def has_step(self):
        """Return whether this EnergyRange has step defined."""
        return self.step is not NO_VALUE

    @property
    def is_adjusted(self):
        """Return True, or raise if this EnergyRange is undefined."""
        if not self.defined:
            raise RuntimeError(f'{self} has undefined items')
        return True

    @property
    def n_energies(self):
        """Return the number of energies in this TheoEnergies."""
        if not self.defined:
            raise RuntimeError(f'{self} has undefined items')
        # +1 because we include both start and stop
        return round((self.stop - self.start) / self.step) + 1

    @property
    def _non_defaults(self):
        """Return the non-default values in self."""
        return [e for e in self if e is not NO_VALUE]

    def adjust_to_fit_step(self):
        """Adjust bounds to fit constraints with steps.

        The base-class implementation complains if this EnergyRange
        is not defined. This method should be extended in subclasses
        to fix inconsistencies of the bounds with the step.

        Raises
        ------
        RuntimeError
            If this method is called before all the `start`, `stop`
            and `step` attributes have been assigned.
        """
        if not self.defined:
            raise RuntimeError(f'{self} has undefined items')

    def contains(self, other, ignore_step=False):
        """Return whether other is a subset of this EnergyRange."""
        if not isinstance(other, EnergyRange):
            raise TypeError
        if not self.has_bounds:
            raise RuntimeError('Cannot compare non-defined TheoEnergies')
        if not other.has_bounds:
            raise ValueError('Cannot compare non-defined TheoEnergies')
        if self.start > other.start or self.stop < other.stop:
            return False
        if ignore_step:
            return True
        if not self.defined:
            raise RuntimeError('Cannot compare non-defined TheoEnergies')
        if not other.defined:
            raise ValueError('Cannot compare non-defined TheoEnergies')
        if self.step != other.step:
            return False
        # Finally, make sure they're not shifted
        self_shift = remainder(self.start, self.step)
        other_shift = remainder(other.start, other.step)
        return abs(self_shift - other_shift) < self.step * EPS

    def copy(self):
        """Return a (deep)copy of this EnergyRange."""
        return self.__class__(*self)

    @classmethod
    def from_value(cls, value):
        """Return an EnergyRange from value."""
        return cls(*value)

    @classmethod
    def from_sorted_grid(cls, energy_grid):
        """Return an energy range from a sorted grid of energies."""
        n_energies = len(energy_grid)
        if n_energies < 2:
            raise ValueError('Not enough energy_grid values. Need '
                             f'at least 2, found {n_energies}')
        start, stop = energy_grid[0], energy_grid[-1]
        step = (stop - start) / (n_energies - 1)
        return cls(start, stop, step)

    def is_equivalent(self, other):
        """Return whether self is equal to other, including swapping."""
        if self == other:
            return True
        # Try converting to an EnergyRange. Notice that we do not
        # use self.__class__ but always EnergyRange, as the base
        # class only swaps the bound values but never adjusts them
        try:
            range_other = EnergyRange(*other)
        except (ValueError, TypeError):
            return False
        return self == range_other

    def intersected(self, other):
        """Return the subset of this range in common with another.

        Parameters
        ----------
        other : EnergyRange
            The other range to intersect this one with.

        Returns
        -------
        intersected : EnergyRange
            The range common to self and other. The step of intersected
            is the same as the one of self. Notice that intersected may
            be wider than the actual common range. This may be the case
            if type(self) needs its bounds to somehow fit its step.

        Raises
        ------
        TypeError
            If other is not an EnergyRange
        ValueError
            If self and other have no energies in common, or if other
            has no bounds defined.
        RuntimeError
            If self has no bounds defined.
        """
        if not isinstance(other, EnergyRange):
            raise TypeError
        if not self.has_bounds:
            raise RuntimeError(f'{self} has no bounds')
        if not other.has_bounds:
            raise ValueError(f'{other} has no bounds')
        new_start = max(self.min, other.min)
        new_stop = min(self.max, other.max)
        if new_stop < new_start:
            raise ValueError('No intersection')
        return self.__class__(new_start, new_stop, self.step)

    @staticmethod
    def parse_string_sequence(string_sequence):
        """Return floating-point or NO_VALUE from a string sequence."""
        # We will use ast to interpret the sequence as a tuple.
        # Notice the trailing comma, in case string_sequence is
        # only one-element-long
        string = ','.join(string_sequence) + ','

        # We have to replace '_' with something that AST can
        # handle. 'None' seems easy enough. We replace it
        # again further down when converting the rest to float
        string = string.replace('_', 'None')
        try:
            return [float(v) if v is not None else NO_VALUE
                    for v in ast.literal_eval(string)]
        except (SyntaxError, ValueError, TypeError,
                MemoryError, RecursionError) as exc:
            new_exc = TypeError if 'float' in exc.args[0] else ValueError
            raise new_exc(' '.join(string_sequence)) from exc

    def set_undefined_values(self, *new_values):
        """Assign undefined values from new_values, then adjust start."""
        if len(new_values) == 1:
            new_values = new_values[0]
        for attr, value in zip(('start', 'stop', 'step'), new_values):
            if getattr(self, attr) is NO_VALUE:
                setattr(self, attr, value)
        self.__post_init__()  # Check and process new values

    def _check_consistency(self):
        """Change inconsistent values or complain."""
        try:
            1 / self.step
        except ZeroDivisionError:
            raise ValueError('Step cannot be zero') from None
        except TypeError:  # NO_VALUE
            pass

        start, stop, step = self
        if self.defined and (stop - start) * step < 0:
            raise ValueError('Inconsistent step. Cannot shift from '
                             f'{start:.2f} to {stop:.2f} with step={step:.2f}')
        if self.has_bounds and self._swap and stop < start:
            self._swap()

    def _swap(self):
        """Swap start and stop, change step."""
        self.start, self.stop = self.stop, self.start
        if self.has_step:
            self.step *= -1


class TheoEnergies(EnergyRange, param='THEO_ENERGIES'):
    """Energy range used for calculating I(V) curves."""

    @property
    def is_adjusted(self):
        """Return whether (stop - start) is a multiple of step."""
        if not self.defined:
            raise RuntimeError(f'{self} has undefined items')
        start, stop, step = self
        return abs(remainder(stop - start,  step)) < EPS

    def adjust_to_fit_step(self):
        """Modify start so that (stop - start) is a multiple of step."""
        # The next one raises RuntimeError if we're not fully defined:
        # We should not mess with the start/stop till we know all.
        if self.is_adjusted:
            return

        start, stop, step = self
        # The next line could in principle also be done with
        # remainder, but there are some corner cases in which
        # it is complicated to get the same results as now.
        start -= step - (stop - start) % step
        # if start < -EPS:    # Testing reveals that this is never hit
            # start = start % step
        if abs(start) < EPS:
            start = step
        self.start = start

    def as_floats(self):
        """Return a list of float values, replacing NO_VALUE with -1."""
        return [-1 if e is NO_VALUE else e for e in self]

    def expanded_by(self, n_steps):
        """Return a copy of this range wider by n_steps left and right."""
        if not self.defined:
            raise RuntimeError(f'{self} has undefined items')
        delta = n_steps * self.step
        start = max(self.step, self.min - delta)
        stop = self.max + delta
        if stop < start:
            raise ValueError(f'Cannot expand by n_steps={n_steps}. '
                             f'{self} would become empty')
        return self.__class__(start, stop, self.step)

    def _check_consistency(self):
        """Change inconsistent values or complain."""
        super()._check_consistency()
        if not all(e > 0 for e in self._non_defaults):
            raise ValueError('Values must be positive')


class IVShiftRange(EnergyRange, param='IV_SHIFT_RANGE'):
    """EnergyRange for Rparams attribute IV_SHIFT_RANGE."""

    @property
    def is_adjusted(self):
        """Return whether bounds are integer multiples of the step."""
        if not self.defined:
            raise RuntimeError(f'{self} has undefined items')
        *bounds, step = self
        return all(abs(remainder(b, step)) < EPS for b in bounds)

    @property
    def is_fixed(self):
        """Return whether this IVShiftRange has equal bounds."""
        if not self.has_bounds:
            raise RuntimeError(f'{self} has no bounds')
        scale = self.step if self.has_step else 1
        return abs(self.stop - self.start) < EPS * scale

    def adjust_to_fit_step(self):
        """Expand bounds so that they are integer multiples of step."""
        # The next one raises RuntimeError if we're not fully defined:
        # We should not mess with the start/stop till we know all.
        if self.is_adjusted:
            return

        # Make sure that adjusting does not un-fix a fixed range
        was_fixed = self.is_fixed

        # Exact floor/ceil comes from the 'exact' solution using
        # Decimal in https://stackoverflow.com/questions/28425705/
        start, stop, step = (Decimal(str(v)) for v in self)
        self.start = float(floor(start / step) * step)
        self.stop = float(ceil(stop / step) * step)

        if was_fixed and not self.is_fixed:
            raise RuntimeError(
                f'step={step!s} cannot be used to fix IV_SHIFT_RANGE to '
                f'{start}: {start} is not an integer multiple of {step}'
                )

    @classmethod
    def fixed(cls, fixed_value):
        """Return an IVShiftRange with both bounds at the same value."""
        return cls(fixed_value, fixed_value, NO_VALUE)

    def set_undefined_step(self, new_step):
        """Assign a new_step if this range has it undefined."""
        super().set_undefined_values(NO_VALUE, NO_VALUE, new_step)
