"""Module linear_operation of viperleed.calc.files.new_displacements.tokens."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)')
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-05-14'
__license__ = 'GPLv3+'

import ast
import re

import numpy as np

from viperleed.calc.constants import DISPLACEMENTS_FILE_EPS

from .base import DisplacementsFileToken, TokenParserError


class LinearOperationTokenParserError(TokenParserError):
    """Class for parsing Errors in the LinearOperationToken."""


class LinearOperationToken(DisplacementsFileToken):
    """Class to parse and represent linear operations for CONSTRAIN blocks.

    Parameters
    ----------
    op_str : str
        The range string to parse.
    """

    def __init__(self, op_str):
        """Construct a LinearOperationToken from a string."""
        # some simple input sanitization
        cleaned = op_str.strip()

        # replace round brackets with square brackets
        cleaned = cleaned.replace('(', '[').replace(')', ']')

        # space separated syntax
        if ',' not in cleaned:
            # add commas between adjacent numbers:
            #  "<num>   <num>" -> "<num>, <num>"
            num_pattern = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
            cleaned = re.sub(rf'({num_pattern})\s+(?={num_pattern})', r'\1, ',
                             cleaned)
            # add commas between brackets
            cleaned = re.sub(r'\]\s+\[', '], [', cleaned)

        # parse to expression
        try:
            parsed = ast.literal_eval(cleaned)
        except (ValueError, SyntaxError) as err:
            msg = f'Could not parse linear operation "{op_str.strip()}".'
            raise LinearOperationTokenParserError(msg) from err

        self.arr = self._validate_array(parsed, op_str.strip())

    @staticmethod
    def _validate_array(values, description):
        """Return a square float array from `values`, or raise.

        Parameters
        ----------
        values : object
            Numeric data to convert. Must have 1, 4, or 9 elements.
        description : object
            Identifies `values` in error messages (a string, or the
            values themselves).

        Returns
        -------
        numpy.ndarray
            A 1x1, 2x2, or 3x3 array.

        Raises
        ------
        LinearOperationTokenParserError
            If `values` is not numeric or has the wrong size.
        """
        try:
            arr = np.array(values, dtype=float)
        except ValueError as err:
            msg = (
                f'Unable to convert linear operation "{description}" to '
                'a numeric array.'
            )
            raise LinearOperationTokenParserError(msg) from err

        # The linear operation must be a linear map between parameters of the
        # same degree of freedom (1, 2, or 3). I.e. it must be a square (1x1),
        # (2x2), or (3x3) matrix.
        if arr.size == 1:
            # single value is a 1x1 matrix
            return arr.reshape((1, 1))
        if arr.size == 4:
            # 2x2 matrix
            return arr.reshape((2, 2))
        if arr.size == 9:
            # 3x3 matrix
            return arr.reshape((3, 3))
        msg = (
            f'Invalid linear operation format: "{description}". '
            'Expected a 1x1, 2x2, or 3x3 matrix.'
        )
        raise LinearOperationTokenParserError(msg)

    @classmethod
    def from_array(cls, arr):
        """Alternate constructor using numeric values directly."""
        inst = cls.__new__(cls)
        inst.arr = cls._validate_array(arr, arr)
        return inst

    def __eq__(self, other):
        """Compare two LinearOperationToken objects for equality."""
        if not isinstance(other, LinearOperationToken):
            return NotImplemented
        return np.allclose(self.arr, other.arr, atol=DISPLACEMENTS_FILE_EPS)

    def __str__(self):
        """Return a string representation of the LinearOperationToken object."""
        return f'LinearOperationToken(arr={self.arr})'
