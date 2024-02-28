"""Test functionality of module woods of viperleed.guilib.leedsim.classes.

Created: 2024-02-20
Author: Michele Riva (@michele-riva)
"""

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.guilib.leedsim.classes import woods


_commensurate = {
    'integer floats': ([[1.0, 2.0], [3.0, 4.0]], True),
    'integers': ([[1, 2], [3, 4]], True),
    'singular': ([[1, 2], [2, 4]], False),
    'non int floats': ([[1.0, 2.0], [3.0, 4.01]], False),
    }

@parametrize('matrix,expect', _commensurate.values(), ids=_commensurate)
def test_is_commensurate(matrix, expect):
    """Check correct result of is_commensurate."""
    commensurate = woods.is_commensurate(matrix)
    assert commensurate is expect


_primes = {
    'int': (12, (2, 2, 3)),
    'int-like float': (42.0, (2, 3, 7)),
    'one': (1, ()),
    'prime': (97, (97,)),
    }

@parametrize('number,expect', _primes.values(), ids=_primes)
def test_prime_factors(number, expect):
    """Check correct splitting of a number into its prime factorization."""
    assert tuple(woods.prime_factors(number)) == expect

_squares = {
    1: (1, 1),
    12: (4, 3),
    16: (16, 1),
    }

@parametrize('number,expect', _squares.items(), ids=_squares)
def test_square_to_prod_of_squares(number, expect):
    """Check correct decomposition in a product of squares."""
    squares, remainder = woods.square_to_prod_of_squares(number)
    assert (squares, remainder) == expect

