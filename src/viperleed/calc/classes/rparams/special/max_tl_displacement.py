"""Module max_tl_displacement of viperleed.calc.classes.rparams.special.

Defines the MaxTLDisplacement class, a float with optional vib value.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-15'
__license__ = 'GPLv3+'

import numpy as np

from dataclasses import dataclass

from ._base import SpecialParameter
from .._defaults import NO_VALUE

@dataclass
class MaxTLDisplacement(SpecialParameter, param='MAX_TL_DISPLACEMENT'):
    """Maximum geometric and vibrational displacements relative to refcalc."""

    geo: float
    _vib: float = NO_VALUE

    def __post_init__(self):
        """Convert non-float inputs and check their range."""
        for attr_name in ('geo', '_vib'):
            attr = getattr(self, attr_name)
            if attr is NO_VALUE:
                continue
            try:
                attr_float = float(attr)
            except (TypeError, ValueError):
                raise TypeError('MAX_TL_DISPLACEMENT value must be '
                                'float') from None
            if attr_float <= 0:
                raise ValueError('MAX_TL_DISPLACEMENT value must be positive')
            setattr(self, attr_name, attr_float)

    @property
    def vib(self):
        """Return the maximum vibrational displacement."""
        return self.geo if self._vib is NO_VALUE else self._vib

    @classmethod
    def from_value(cls, value):
         """Return a MaxTLDisplacement from a 2-item tuple."""
         return cls(*value)

    def is_too_far(self, atom):
          """Return whether `atom` was displaced too much."""
          if atom.distance(atom.oriState) > self.geo:
              return True
          if any(np.abs(atom.site.vibamp[el] - atom.site.oriState.vibamp[el])
                 > self.vib for el in atom.site.vibamp):
              return True
          return False
