"""Module offset of viperleed.calc.files.new_displacements.tokens."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
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

    def __init__(self, offset_str):
        """Construct an OffsetToken from a string.

        Parameters
        ----------
        offset_str : str
            One to three whitespace-separated numbers.
        """
        parts = offset_str.strip().split()
        self.offset = self._validate_offsets(parts, offset_str)

    @property
    def dof(self):
        """Return the number of degrees of freedom (DOF) for the offset."""
        return self.offset.size

    @classmethod
    def from_floats(cls, *offsets):
        """Alternate constructor using 1-3 numeric values directly.

        Parameters
        ----------
        *offsets : float
            The offset values.
        """
        inst = cls.__new__(cls)
        inst.offset = cls._validate_offsets(offsets, offsets)
        return inst

    @staticmethod
    def _validate_offsets(parts, description):
        """Return `parts` as a float array of 1-3 values, or raise.

        Parameters
        ----------
        parts : sequence
            Numbers, or strings representing numbers.
        description : object
            Identifies `parts` in error messages (a string, or the
            values themselves).

        Returns
        -------
        numpy.ndarray
            The offset values.

        Raises
        ------
        OffsetTokenParserError
            If there are not 1 to 3 values, or any is not numeric.
        """
        if not 1 <= len(parts) <= 3:
            msg = ('Invalid offset format (expect 1-3 floats): '
                   f'"{description}".')
            raise OffsetTokenParserError(msg)
        try:
            return np.array(parts, dtype=float)
        except ValueError as err:
            msg = f'Non-numeric value in offset: "{description}".'
            raise OffsetTokenParserError(msg) from err

    def __eq__(self, other):
        """Compare two OffsetToken objects for equality."""
        if not isinstance(other, OffsetToken):
            return NotImplemented
        if self.dof != other.dof:
            return False
        return np.allclose(
            self.offset, other.offset, atol=DISPLACEMENTS_FILE_EPS
        )

    def __str__(self):
        """Return a string representation of the OffsetToken object."""
        return f'OffsetToken(offset={self.offset})'
