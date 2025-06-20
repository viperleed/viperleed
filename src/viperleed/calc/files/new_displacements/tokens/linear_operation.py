"""Module for the <linear_operation> token in the DISPLACEMENTS file."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-05-14"
__license__ = "GPLv3+"

import re
import ast

import numpy as np


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


    def __init__(self, op_str: str):
        """Construct a LinearOperationToken from a string."""
        # some simple input sanitization
        cleaned = op_str.strip()

        # space separated syntax
        if ',' not in cleaned:
            # add commas between numbers
            cleaned = re.sub(r'(?<=\d)\s+(?=\d)', ', ', cleaned)
            # add commas between brackets
            cleaned = re.sub(r'\]\s+\[', '], [', cleaned)
        # replace round brackets with square brackets
        cleaned = cleaned.replace('(', '[').replace(')', ']')

        # parse to expression
        try:
            parsed = ast.literal_eval(cleaned)
        except (ValueError, SyntaxError) as err:
            msg = f'Could not parse linear operation "{op_str.strip()}".'
            raise LinearOperationTokenParserError(msg) from err

        try:
            arr = np.array(parsed, dtype=float)
        except ValueError as err:
            msg = (f'Unable to convert linear operation "{op_str.strip()}" to '
                   'a numeric array.')
            raise LinearOperationTokenParserError(msg) from err

        # The linear operation must be a linear map between parameters of the
        # same degree of freedom (1, 2, or 3). I.e. it must be a square (1x1),
        # (2x2), or (3x3) matrix.
        if arr.size == 1:
            # single value is a 1x1 matrix
            arr = arr.reshape((1, 1))
        elif arr.size == 4:
            # 2x2 matrix
            arr = arr.reshape((2, 2))
        elif arr.size == 9:
            # 3x3 matrix
            arr = arr.reshape((3, 3))
        else:
            msg = (
                f'Invalid linear operation format: "{op_str.strip()}". '
                'Expected a 1x1, 2x2, or 3x3 matrix.'
            )
            raise LinearOperationTokenParserError(msg)
        self.arr = arr

    @classmethod
    def from_array(cls, arr) -> 'LinearOperationToken':
        """Alternate constructor using numeric values directly."""
        inst = cls.__new__(cls)
        inst.arr = np.array(arr, dtype=float)
        return inst


    def __eq__(self, other):
        """Compare two LinearOperationToken objects for equality."""
        if not isinstance(other, LinearOperationToken):
            return NotImplemented
        if self.arr.shape != other.arr.shape:
            return False
        return np.allclose(self.arr, other.arr)

    def __str__(self):
        """Return a string representation of the LinearOperationToken object."""
        return f'LinearOperationToken(arr={self.arr})'
