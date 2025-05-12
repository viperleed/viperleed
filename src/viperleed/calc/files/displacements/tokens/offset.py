"""Offset token for DISPLACEMENT file parsing."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-05-12'

from .base import TokenParserError


class OffsetTokenParserError(TokenParserError):
    """Class raised by OffsetToken."""


class OffsetToken:
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

    _EPS = 1e-6

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
        if type(offset) is not float:
            raise OffsetTokenParserError('from_float require float argument.')
        inst = cls.__new__(cls)
        inst.offset = offset
        return inst

    def __eq__(self, other):
        """Compare two OffsetToken objects for equality."""
        if not isinstance(other, OffsetToken):
            return False
        return abs(self.offset - other.offset) < self._EPS

    def __repr__(self):
        """Return a string representation of the OffsetToken object."""
        return f'OffsetToken(offset={self.offset})'
