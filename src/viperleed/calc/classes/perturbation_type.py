"""Module perturbation_type."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'


from enum import Enum


class PerturbationType(str, Enum):
    """Enum for perturbation types."""

    GEO = 'geo'
    VIB = 'vib'
    OCC = 'occ'

    def __str__(self):
        """Return the string representation of the enum."""
        return self.value

    @classmethod
    def from_string(cls, s: str) -> 'PerturbationType':
        """Convert a string to a PerturbationType enum."""
        try:
            return cls(s)
        except ValueError:
            msg = f'Unknown offset type: {s}'
            raise ValueError(msg) from None
