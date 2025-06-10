"""Module for the <element> token in the DISPLACEMENTS file."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-05-13'

from viperleed.calc.lib import periodic_table

from .base import DisplacementsFileToken, TokenParserError


class ElementTokenParserError(TokenParserError):
    """Class for parsing Errors in the ElementToken."""


class ElementToken(DisplacementsFileToken):
    """Class to parse and represent the <element> token.

    The <element> token in the DISPLACEMENTS file is used to specify chemical
    (occupational) perturbations and ranges. It is used in the OCC_DELTA block
    and in lines relating to chemical parameters in the OFFSETS and CONSTRAIN
    blocks.
    The <element> token is specified by the user as a simple (case-insensitive)
    string corresponding to the Symbol of a element from the periodic table.

    Parameters
    ----------
    element_str : str
        The string containing the element to be parsed.
    """

    def __init__(self, element_str: str):
        """Construct an ElementToken from a string."""
        _element_str = element_str.strip()
        # first get atomic number
        try:
            self.atomic_number = periodic_table.get_atomic_number(_element_str)
        except ValueError as err:
            msg = f'Could not parse chemical element "{_element_str}".'
            raise ElementTokenParserError(msg) from err
        # also get symbol
        self.symbol = periodic_table.get_element_symbol(self.atomic_number)

    def __eq__(self, other):
        """Compare two RangeToken objects for equality."""
        if not isinstance(other, ElementToken):
            return False
        return other.atomic_number == self.atomic_number

    def __repr__(self):
        """Return a string representation of the ElementToken object."""
        return f'ElementToken(element={self.symbol})'

