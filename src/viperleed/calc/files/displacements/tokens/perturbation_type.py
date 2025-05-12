"""Module perturbation_type."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'


from enum import Enum


class PerturbationTypeError(ValueError):
    """Error class for unknown perturbation type errors."""

class PerturbationType(str, Enum):
    """Enum for perturbation types."""

    GEO = 'geo'
    VIB = 'vib'
    OCC = 'occ'
    DOM = 'dom'

    def __str__(self):
        """Return the string representation of the enum."""
        return self.value

    @classmethod
    def from_string(cls, s: str) -> 'PerturbationType':
        """Convert a string to a PerturbationType enum."""
        try:
            return cls(s)
        except ValueError:
            msg = (f'Unknown perturbation type tag: {s}. Tag must be one of '
                   '"geo", "vib", "occ", or "dom"."')
            raise PerturbationTypeError(msg) from None
