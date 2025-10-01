"""Module for the <element> token in the DISPLACEMENTS file."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-07-18'
__license__ = 'GPLv3+'


from viperleed.calc.constants import DISPLACEMENTS_FILE_EPS

from .base import DisplacementsFileToken, TokenParserError


class TotalOccupationTokenParserError(TokenParserError):
    """Class for parsing Errors in the TotalOccupationToken."""


class TotalOccupationToken(DisplacementsFileToken):
    """Class to parse and represent the total occupation token.

    The total occupation token is used to specify a constant sum of the
    occupations of all elements on a site, allowing the individual occupations
    to vary freely. This can be used, for example, to vary the occupation in
    an alloy without introducing vacancies. The <total_occupation> token is
    given as a single floating-point number that must be in the range
    [0, 1].

    Parameters
    ----------
    total_occupation_str : str
        The string containing the total occupation to be parsed.
    """

    def __init__(self, total_occupation_str: str):
        """Construct a TotalOccupationToken from a string."""
        try:
            total_occupation_value = float(total_occupation_str.strip())
        except ValueError as err:
            msg = (
                'Non-numeric value in total occupation: '
                f'"{total_occupation_str}"'
            )
            raise TotalOccupationTokenParserError(msg) from err
        if not (0 <= total_occupation_value <= 1):
            msg = (
                f'Total occupation must be in the range [0, 1], '
                f'got {total_occupation_value}.'
            )
            raise TotalOccupationTokenParserError(msg)
        self.total_occupation = total_occupation_value

    def __eq__(self, other):
        """Compare two TotalOccupationToken objects for equality."""
        if not isinstance(other, TotalOccupationToken):
            return NotImplemented
        return (
            self.total_occupation - other.total_occupation
        ) < DISPLACEMENTS_FILE_EPS

    def __str__(self):
        """Return a string representation of the TotalOccupationToken object."""
        return (
            f'TotalOccupationToken(total_occupation={self.total_occupation})'
        )
