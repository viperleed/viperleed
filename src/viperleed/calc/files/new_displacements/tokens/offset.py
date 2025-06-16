"""Module for the <offset> token in the DISPLACEMENTS file."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-05-12"
__license__ = "GPLv3+"

from .base import DisplacementsFileToken, TokenParserError
from viperleed.calc.files.new_displacements import DISPLACEMENTS_FILE_EPS


class OffsetTokenParserError(TokenParserError):
    """Class raised by OffsetToken."""


class OffsetToken(DisplacementsFileToken):
    """Class to parse and represent offsets in the DISPLACEMENTS file.

    Ranges are provided by the user as strings containing a single float.

    Parameters
    ----------
    offset_str : str
        The offset string to parse.

    Raises
    ------
    OffsetTokenParserError
        If parsing fails.
    """

    def __init__(self, offset_str: str):
        """Construct a OffsetToken from a string."""
        parts = offset_str.strip().split()
        if len(parts) != 1:
            msg = (
                f'Invalid offset format: "{offset_str}". Expected format: '
                '"<offset>".'
            )
            raise OffsetTokenParserError(msg)

        try:
            offset = float(parts[0])
        except ValueError as err:
            msg = f'Non-numeric value in offset: "{offset_str}"'
            raise OffsetTokenParserError(msg) from err

        self.offset = offset

    @classmethod
    def from_floats(cls, offset: float) -> 'OffsetToken':
        """Alternate constructor using numeric values directly."""
        inst = cls.__new__(cls)
        try:
            inst.offset = offset = float(offset)
        except ValueError as err:
            msg = f'Non-numeric value in offset: "{offset}"'
            raise OffsetTokenParserError(msg) from err
        return inst

    def __eq__(self, other):
        """Compare two OffsetToken objects for equality."""
        if not isinstance(other, OffsetToken):
            return NotImplemented
        return abs(self.offset - other.offset) < DISPLACEMENTS_FILE_EPS

    def __str__(self):
        """Return a string representation of the OffsetToken object."""
        return f'OffsetToken(offset={self.offset})'
