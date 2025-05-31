"""Tests for module utils of viperleed.gui.leedsim.classes.woods."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-02-20'
__license__ = 'GPLv3+'

from pytest_cases import parametrize

from viperleed.gui.leedsim.classes.woods.utils import is_commensurate
from viperleed.gui.leedsim.classes.woods.utils import prime_factors
from viperleed.gui.leedsim.classes.woods.utils import square_to_prod_of_squares


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
