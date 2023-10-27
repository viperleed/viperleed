"""Module theo_energies of viperleed.tleedmlib.classes.rparams.special.

Created on 2023-10-27

@author: Michele Riva (@michele-riva)

Defines the TheoEnergies subclass of EnergyRange, used as the Rparams
attribute THEO_ENERGIES.
"""

from dataclasses import dataclass
from math import remainder

from .energy_range import EnergyRange
from .._defaults import NO_VALUE


@dataclass(eq=False)  # __eq__ inherited from EnergyRange
class TheoEnergies(EnergyRange):
    """Energy range used for calculating I(V) curves."""

    def __post_init__(self):
        """Check and process initialization values."""
        super().__post_init__()
        if not all(e > 0 for e in self._non_defaults):
            raise ValueError('Values must be positive')
        if self.has_bounds and self.stop < self.start:
            raise ValueError('Maximum energy value should be at '
                             'least as large as the minimum')

        # Mess with start/stop only if all the values are present,
        # otherwise, leave it for when the others will be initialized
        # from experimental data. That's in Rparams.initTheoEnergies.
        if self.defined:
            self.adjust_to_fit_step()

    @property
    def is_adjusted(self):
        """Return whether (stop - start) is a multiple of step."""
        if not self.defined:
            raise RuntimeError(f'{self} has undefined items')
        start, stop, step = self
        return abs(remainder(stop - start,  step)) < 1e-6

    @property
    def n_energies(self):
        """Return the number of energies in this TheoEnergies."""
        if not self.defined:
            raise RuntimeError(f'{self} has undefined items')
        # +1 because we include both start and stop
        return round((self.stop - self.start) / self.step) + 1

    def adjust_to_fit_step(self):
        """Modify start so that (stop - start) is a multiple of step."""
        if self.is_adjusted:
            return

        start, stop, step = self
        # The next line could in principle also be done with
        # remainder, but there are some corner cases in which
        # it is complicated to get the same results as now.
        start -= step - (stop - start) % step
        # if start < -1e-6:    # Testing reveals that this is never hit
            # start = start % step
        if abs(start) < 1e-6:
            start = step
        self.start = start

    def as_floats(self):
        """Return a list of float values, replacing NO_VALUE with -1."""
        return [-1 if e is NO_VALUE else e for e in self]

    def contains(self, other):
        """Return whether other is a subset of this TheoEnergies."""
        if not isinstance(other, TheoEnergies):
            raise TypeError
        if not self.defined:
            raise RuntimeError('Cannot compare non-defined TheoEnergies')
        if not other.defined:
            raise ValueError('Cannot compare non-defined TheoEnergies')
        if self.start > other.start or self.stop < other.stop:
            return False
        if self.step != other.step:
            return False
        # Finally, make sure they're not shifted
        self_shift = remainder(self.start, self.step)
        other_shift = remainder(other.start, other.step)
        return abs(self_shift - other_shift) < self.step * 1e-6

    def set_undefined_values(self, new_values):
        """Assign undefined values from new_values, then adjust start."""
        for attr, value in zip(('start', 'stop', 'step'), new_values):
            if getattr(self, attr) is NO_VALUE:
                setattr(self, attr, value)
        self.adjust_to_fit_step()

    _swap = None  # Never swap a TheoEnergies. All items must be > 0
