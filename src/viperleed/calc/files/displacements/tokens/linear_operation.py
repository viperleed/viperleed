"""Linear operation token."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-05-14'

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
    range_str : str
        The range string to parse.
    """

    _EPS = 1e-6

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

        # parse to expression
        try:
            parsed = ast.literal_eval(cleaned)
            arr = np.array(parsed, dtype=float)
        except (ValueError, SyntaxError) as err:
            msg = f'Could not parse linear operation "{op_str.strip()}".'
            raise LinearOperationTokenParserError(msg) from err
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
            return False
        if self.arr.shape != other.arr.shape:
            return False
        return np.allclose(self.arr, other.arr)

    def __repr__(self):
        """Return a string representation of the LinearOperationToken object."""
        return f'LinearOperationToken(arr={self.arr})'
