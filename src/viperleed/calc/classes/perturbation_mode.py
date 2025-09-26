"""Module perturbation_mode of viperleed.calc.classes."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-13'
__license__ = 'GPLv3+'

from enum import Enum


class PerturbationModeError(ValueError):
    """Error class for unknown perturbation type errors."""


class PerturbationMode(str, Enum):
    """Enum for perturbation types."""

    GEO = 'geo'
    VIB = 'vib'
    OCC = 'occ'
    DOM = 'dom'

    def __str__(self):
        """Return the string representation of the enum."""
        return self.value

    @classmethod
    def from_string(cls, mode_string: str):
        """Convert a string to a PerturbationType enum."""
        try:
            return cls(mode_string)
        except ValueError:
            msg = f'Unknown perturbation mode: {mode_string!r}. '
            raise PerturbationModeError(msg) from None
