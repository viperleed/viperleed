"""Test functionality of module woods of viperleed.guilib.leedsim.classes.

Created: 2024-02-20
Author: Michele Riva (@michele-riva)
"""

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.guilib.leedsim.classes.woods import is_commensurate
from viperleed.guilib.leedsim.classes.woods import prime_factors
from viperleed.guilib.leedsim.classes.woods import square_to_prod_of_squares
from viperleed.guilib.leedsim.classes.woods import Woods
from viperleed.guilib.leedsim.classes.woods import MatrixIncommensurateError
from viperleed.guilib.leedsim.classes.woods import WoodsNotRepresentableError
from viperleed.guilib.leedsim.classes.woods import WoodsSyntaxError


_commensurate = {
    'integer floats': ([[1.0, 2.0], [3.0, 4.0]], True),
    'integers': ([[1, 2], [3, 4]], True),
    'singular': ([[1, 2], [2, 4]], False),
    'non int floats': ([[1.0, 2.0], [3.0, 4.01]], False),
    }

@parametrize('matrix,expect', _commensurate.values(), ids=_commensurate)
def test_is_commensurate(matrix, expect):
    """Check correct result of is_commensurate."""
    commensurate = is_commensurate(matrix)
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
    assert tuple(prime_factors(number)) == expect

_squares = {
    1: (1, 1),
    12: (4, 3),
    16: (16, 1),
    }

@parametrize('number,expect', _squares.items(), ids=_squares)
def test_square_to_prod_of_squares(number, expect):
    """Check correct decomposition in a product of squares."""
    squares, remainder = square_to_prod_of_squares(number)
    assert (squares, remainder) == expect


class TestWoodsRaises:
    """Collection of tests for exceptions raised by the Woods class."""

    _init = {
        'invalid style': ({'style': 'invalid'}, ValueError),
        'non-string string': ({'string': 1}, TypeError),
        'not a wood notation': ({'string': '23'}, WoodsSyntaxError),
        'gamma syntax invalid character': (
            {'string': 'c(2, x 12)'},
            WoodsSyntaxError
            ),
        'gamma syntax unmatched': (
            {'string': 'c((2 x 12)'},
            WoodsSyntaxError
            ),
        'unsupported math': ({'string': 'cos(2)x3'}, WoodsSyntaxError),
        'gamma not int': ({'string': '1.3x8'}, ValueError),
        'missing bulk_basis': ({'matrix': np.eye(2)}, TypeError),
        'matrix shape': (
            {'matrix': np.eye(3), 'bulk_basis': np.eye(2)},
            ValueError
            ),
        'bulk_basis shape': (
            {'matrix': np.eye(2), 'bulk_basis': np.eye(3)},
            ValueError
            ),
        'inconsistent matrix and string': (
            {'matrix': np.eye(2), 'string': '2x2', 'bulk_basis': np.eye(2)},
            ValueError
            ),
        'incommensurate matrix': (
            {'matrix': np.eye(2)*1.2, 'bulk_basis': np.eye(2)},
            MatrixIncommensurateError
            ),
        'matrix not representable': (
            {'matrix': np.arange(4).reshape(2, 2), 'bulk_basis': np.eye(2)},
            WoodsNotRepresentableError
            )
        }

    @parametrize('kwargs,exc', _init.values(), ids=_init)
    def test_invalid_init_args(self, kwargs, exc):
        """Check that exc is raised with invalid kwargs at __init__."""
        with pytest.raises(exc):
            Woods(**kwargs)

