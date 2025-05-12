"""Module perturbation_type."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'


from enum import Enum

from .base import DisplacementsFileToken, TokenParserError

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


class TypeTokenParserError(TokenParserError):
    """Class for errors during the TypeToken parsing."""

class TypeToken(DisplacementsFileToken):
    """Class for the TypeToken."""

    def __init__(self, type_str):
        self.type_str = type_str
        try:
            self.type = PerturbationType(type_str)
        except PerturbationTypeError as err:
            msg = f'Unable to parse perturbation type {self.type_str}.'
            raise TypeTokenParserError(msg) from err

    def __eq__(self, other):
        """Compare self to other TypeToken."""
        return self.type is other.type

    def __repr__(self):
        """Return representation of TypeToken."""
        return f'TypeToken({self.type})'
