"""Module perturbation_type of viperleed_calc.classes."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-05-13"
__license__ = "GPLv3+"

from enum import Enum
from viperleed.calc.lib.string_utils import harvard_commas


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
    def from_string(cls, s: str):
        """Convert a string to a PerturbationType enum."""
        try:
            return cls(s)
        except ValueError:
            msg = harvard_commas((repr(e.value) for e in cls), sep='or')
            raise PerturbationModeError(msg) from None
