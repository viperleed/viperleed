"""Module symmetry_eps of viperleed.tleedmlib.classes.rparams.special.

Created on 2023-12-11

@author: Alexander Imre (@amimre)

Defines the class SymmetryEps, which is a float with optional z value.
"""
from dataclasses import dataclass
from typing import Optional

from ._base import SpecialParameter

@dataclass
class SymmetryEps(float, SpecialParameter):
    """SymmetryEps acts like a float but has an optional .Z value.

    Used for interpreting the parameter SYMMETRY_EPS. The user can specify a
    second float value to be used for symmetry comparisons in the z direction.
    If no second value is given, the first value is used by default."""

    value: float
    z: Optional[float] = None

    def __new__(cls, value, z=None):
        """Initialize instance."""
        return float.__new__(cls, value)

    def __init__(self, value, z=None):
        float.__init__(self)
        self.__z = z

    @property
    def Z(self):
        """Return z value of instance."""
        if self.__z is None:
            return float(self)
        return self.__z
