"""Tests for module utils of viperleed.guilib.leedsim.classes.woods.

Created: 2024-02-20
Author: Michele Riva (@michele-riva)
"""

from pytest_cases import parametrize

from viperleed.guilib.leedsim.classes.woods.utils import (
    is_commensurate,
    prime_factors,
    square_to_prod_of_squares,
    )


_COMMENSURATE = {  # pylint: disable=R6101
    'integer floats': ([[1.0, 2.0], [3.0, 4.0]], True),
    'integers': ([[1, 2], [3, 4]], True),
    'singular': ([[1, 2], [2, 4]], False),
    'non int floats': ([[1.0, 2.0], [3.0, 4.01]], False),
    }

@parametrize('matrix,expect', _COMMENSURATE.values(), ids=_COMMENSURATE)
def test_is_commensurate(matrix, expect):
    """Check correct result of is_commensurate."""
    commensurate = is_commensurate(matrix)
    assert commensurate is expect


_PRIMES = {  # pylint: disable=R6101
    'int': (12, (2, 2, 3)),
    'int-like float': (42.0, (2, 3, 7)),
    'one': (1, ()),
    'prime': (97, (97,)),
    }

@parametrize('number,expect', _PRIMES.values(), ids=_PRIMES)
def test_prime_factors(number, expect):
    """Check correct splitting of a number into its prime factorization."""
    assert tuple(prime_factors(number)) == expect


_SQUARES = {  # pylint: disable=R6101
    1: (1, 1),
    12: (4, 3),
    16: (16, 1),
    }

@parametrize('number,expect', _SQUARES.items(), ids=_SQUARES)
def test_square_to_prod_of_squares(number, expect):
    """Check correct decomposition in a product of squares."""
    squares, remainder = square_to_prod_of_squares(number)
    assert (squares, remainder) == expect
