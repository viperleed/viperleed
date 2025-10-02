"""Module offset of viperleed.calc.files.new_displacements.tokens."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-12'
__license__ = 'GPLv3+'

import numpy as np

from viperleed.calc.constants import DISPLACEMENTS_FILE_EPS

from .base import DisplacementsFileToken, TokenParserError


class OffsetTokenParserError(TokenParserError):
    """Class raised by OffsetToken."""


class OffsetToken(DisplacementsFileToken):
    """Class to parse and represent offsets in the DISPLACEMENTS file.

    Ranges are provided by the user as strings containing 1 to 3 floats.

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
        """Construct an OffsetToken from a string."""
        parts = offset_str.strip().split()
        if not (1 <= len(parts) <= 3):
            msg = f'Invalid offset format (expect 1–3 floats): "{offset_str}".'
            raise OffsetTokenParserError(msg)

        try:
            offset_values = np.array([float(p) for p in parts], dtype=float)
        except ValueError as err:
            msg = f'Non-numeric value in offset: "{offset_str}".'
            raise OffsetTokenParserError(msg) from err

        self.offset = offset_values

    @property
    def dof(self):
        """Return the number of degrees of freedom (DOF) for the offset."""
        return self.offset.size

    @classmethod
    def from_floats(cls, *offsets: float) -> 'OffsetToken':
        """Alternate constructor using 1–3 numeric values directly."""
        if not (1 <= len(offsets) <= 3):
            msg = f'Expected 1 to 3 offset values, got {len(offsets)}.'
            raise OffsetTokenParserError(msg)

        inst = cls.__new__(cls)
        try:
            inst.offset = np.array([float(o) for o in offsets], dtype=float)
        except ValueError as err:
            msg = f'Non-numeric value in offset: "{offsets}"'
            raise OffsetTokenParserError(msg) from err
        return inst

    def __eq__(self, other):
        """Compare two OffsetToken objects for equality."""
        if not isinstance(other, OffsetToken):
            return NotImplemented
        if self.dof != other.dof:
            return False
        return abs(self.offset - other.offset) < DISPLACEMENTS_FILE_EPS

    def __str__(self):
        """Return a string representation of the OffsetToken object."""
        return f'OffsetToken(offset={self.offset})'
