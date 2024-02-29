"""Test functionality of module woods of viperleed.guilib.leedsim.classes.

Created: 2024-02-20
Author: Michele Riva (@michele-riva)
"""

from collections import namedtuple

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


SQUARE = (1, 0), (0, 1)
HEX = (1, 0), (-0.5, 3**0.5/2)

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


_WoodsArgs = namedtuple('_WoodsArgs', ('string', 'bulk_basis', 'matrix'))


class TestFromAndToMatrix:
    """Collection of tests for initializing a Woods from matrix."""

    _matrix = {
        '1x1': _WoodsArgs('p(1×1)', SQUARE, np.eye(2)),
        'rt5 R+27': _WoodsArgs('p(√5×√5)R26.6°', SQUARE, ((2, 1), (-1, 2))),
        'rt5 R+63': _WoodsArgs('p(√5×√5)R63.4°', SQUARE, ((1, 2), (-2, 1))),
        'rt5 R-27': _WoodsArgs('p(√5×√5)R-26.6°', SQUARE, ((2, -1), (1, 2))),
        'rt5 R-117': _WoodsArgs('p(√5×√5)R-116.6°',
                                SQUARE, ((-1, -2), (2, -1))),
        'rt31 R9': _WoodsArgs('p(√31×√31)R8.9°', HEX, ((6, 1), (-1, 5))),
        'rt31 R-51': _WoodsArgs('p(√31×√31)R-51.1°', HEX, ((1, -5), (5, 6))),
        'rt67 R12': _WoodsArgs('p(√67×√67)R12.2°', HEX, ((9, 2), (-2, 7))),
        'c2x4': _WoodsArgs('c(2×4)', SQUARE, ((1, 2), (-1, 2))),
        }

    @parametrize(args=_matrix.values(), ids=_matrix)
    def test_from_matrix(self, args, subtests):
        """Check correct initialization from a matrix."""
        with subtests.test('matrix at __init__'):
            woods = Woods(matrix=args.matrix, bulk_basis=args.bulk_basis)
            assert woods.string == args.string
        with subtests.test('matrix at from_matrix call'):
            woods = Woods()
            woods.from_matrix(args.matrix, bulk_basis=args.bulk_basis)
            assert np.all(woods.bulk_basis == args.bulk_basis)

    @parametrize(args=_matrix.values(), ids=_matrix)
    def test_to_matrix_attribute(self, args):
        """Check correct identification of matrix given a Woods string."""
        woods = Woods(string=args.string, bulk_basis=args.bulk_basis)
        matrix = woods.matrix
        assert np.all(matrix == args.matrix)
        assert woods.matrix is matrix

    # Skip those cases that we know require to explicitly fix the
    # angle and would otherwise raise a MatrixIncommensurateError
    _to_matrix = {k: v for k, v in _matrix.items() if 'rt31' not in k}

    @parametrize(args=_to_matrix.values(), ids=_to_matrix)
    def test_to_matrix_method(self, args):
        """Check correct identification of matrix given a Woods string."""
        # Notice that the to_matrix method, compared to the .matrix
        # attribute, does not attempt to fix small angle mismatches
        woods = Woods()
        woods.to_matrix(args.string, bulk_basis=args.bulk_basis)
        assert np.all(woods.bulk_basis == args.bulk_basis)

    _cleared = {
        'bulk_basis': HEX,
        'string': 'c4x2',
        }

    @parametrize('attr,val', _cleared.items(), ids=_cleared)
    def test_matrix_cleared(self, attr, val):
        """Check that setting a new string value clears the stored matrix."""
        woods = Woods('1x1', matrix=np.eye(2), bulk_basis=SQUARE)
        assert woods._Woods__matrix is not None
        setattr(woods, attr, val)
        assert woods._Woods__matrix is None


class TestFromAndToString:
    """Collection of tests for initializing a Woods from string."""

    _parse = {
        '1x1': ('p', 1, 1, 0),
        'r 3xrt3R30': ('p', 3**.5, 3**.5, 30),
        'rt3xrt12R30': ('p', 3**.5, 12**.5, 30),
        'c4×3*4': ('c', 4, 12, 0),
        }
    _string = {
        '1x1': {'ascii': 'p(1x1)', 'unicode': 'p(1×1)'},
        'rt3xrt3R30': {'ascii': 'p(sqrt3 x sqrt3)R30.0',                        # TODO: not nice to have decimals with so round angles
                       'unicode': 'p(√3×√3)R30.0°'},
        'rt3xrt12R30': {'ascii': 'p(sqrt3 x 2sqrt3)R30.0',
                        'unicode': 'p(√3×2√3)R30.0°'},
        'c4×3*4': {'ascii': 'c(4x12)', 'unicode': 'c(4×12)'},
        }

    @parametrize('woods,expect', _parse.items(), ids=_parse)
    def test_parse(self, woods, expect):
        """Check expected outcome of Woods.parse."""
        parsed = Woods.parse(woods)
        assert parsed == expect

    @parametrize('string,style_and_result', _string.items(), ids=_string)
    def test_string(self, string, style_and_result, subtests):
        """Check expected result of Woods.string."""
        for style, expect in style_and_result.items():
            woods = Woods(string, style=style)
            with subtests.test(style):
                assert woods.string == expect


class TestRaises:
    """Collection of tests for exceptions raised by the Woods class."""

    _fmt = {
        'too_long': 'last_char_is_an_acceptable_u',
        'invalid': 'd',
        }

    @parametrize(fmt_spec=_fmt.values(), ids=_fmt)
    def test_format_invalid(self, fmt_spec):
        """Check complaints with an invalid format specification."""
        woods = Woods('p(1X1)')
        with pytest.raises(TypeError):
            format(woods, fmt_spec)

    def test_guess_rotation_invalid(self):
        """Check complaints for an invalid rotation angle."""
        woods = Woods('rt5xrt5R45', bulk_basis=SQUARE)
        with pytest.raises(MatrixIncommensurateError):
            woods.guess_correct_rotation()

    _init = {
        'invalid style': ({'style': 'invalid'}, ValueError),
        'invalid style type': ({'style': 3}, TypeError),
        'non-string string': ({'string': 1}, TypeError),
        'not a wood notation': ({'string': '23'}, WoodsSyntaxError),
        'gamma syntax invalid character': (
            {'string': 'c(2, x 12)'},
            WoodsSyntaxError
            ),
        'gamma syntax unmatched': ({'string': 'c((2 x 12)'}, WoodsSyntaxError),
        'too many gammas': ({'string': '2x3x4'}, WoodsSyntaxError),
        'unsupported math': ({'string': 'cos(2)x3'}, WoodsSyntaxError),
        'gamma not int': ({'string': '1.3x8'}, ValueError),
        'missing bulk_basis': ({'matrix': np.eye(2)}, TypeError),
        'matrix shape': (
            {'matrix': np.eye(3), 'bulk_basis': SQUARE},
            ValueError
            ),
        'bulk_basis shape': (
            {'matrix': np.eye(2), 'bulk_basis': np.eye(3)},
            ValueError
            ),
        'inconsistent matrix and string': (
            {'matrix': np.eye(2), 'string': '2x2', 'bulk_basis': SQUARE},
            ValueError
            ),
        'incommensurate matrix': (
            {'matrix': np.eye(2)*1.2, 'bulk_basis': SQUARE},
            MatrixIncommensurateError
            ),
        'matrix not representable': (
            {'matrix': np.arange(4).reshape(2, 2), 'bulk_basis': SQUARE},
            WoodsNotRepresentableError
            )
        }

    @parametrize('kwargs,exc', _init.values(), ids=_init)
    def test_invalid_init_args(self, kwargs, exc):
        """Check that exc is raised with invalid kwargs at __init__."""
        with pytest.raises(exc):
            Woods(**kwargs)

    def test_invalid_parse(self):
        """Check complaints with an invalid argument to Woods.parse."""
        # Most of the invalid cases are already tested in
        # test_invalid_init_args. We only miss a TypeError
        woods = Woods()
        with pytest.raises(TypeError):
            woods.parse(1)

    def test_from_matrix_no_basis(self):
        """Check complaints when assigning a Woods from matrix."""
        woods = Woods()
        with pytest.raises(ValueError):
            woods.from_matrix(np.eye(2))

    def test_to_matrix_empty_string(self):
        """Check complaints of to_matrix method with an empty woods."""
        woods = Woods()
        with pytest.raises(WoodsSyntaxError):
            woods.to_matrix()

    def test_to_matrix_invalid_basis(self):
        """Check complaints of to_matrix method without a bulk_basis."""
        woods = Woods(string='p(1x1)')
        with pytest.raises(ValueError):
            woods.to_matrix()
    
    def test_matrix_attr_not_commensurate(self):
        """Check complaints accessing .matrix when incommensurate."""
        woods = Woods('rt3xrt3R30', bulk_basis=SQUARE)
        with pytest.raises(MatrixIncommensurateError):
            woods.matrix


class TestStrReprFormat:
    """Collection of tests for string-ifying a Woods."""

    _str = {'unicode': 'p(1×1)', 'ascii': 'p(1x1)'}
    _repr = {
        'no basis, unicode': (_WoodsArgs('rt2xrt8R45', None, None), 'unicode',
                              "Woods('p(√2×2√2)R45.0°')"),
        'no basis, ascii': (_WoodsArgs('rt2xrt8R45', None, None), 'ascii',
                            "Woods('p(sqrt2 x 2sqrt2)R45.0', style='ascii')"),
        'basis': (_WoodsArgs('3x9', SQUARE, None), 'a',
                  "Woods('p(3x9)', bulk_basis=[[1,0], [0,1]], style='ascii')"),
        'matrix': (_WoodsArgs('2x7', SQUARE, np.diag((2, 7))), 'ascii',
                   "Woods('p(2x7)', bulk_basis=[[1,0], [0,1]], "
                   "matrix=[[2,0], [0,7]], style='ascii')"),
        }
    _fmt = {
        '5x4, ascii': ('5x4', 'a', 'p(5x4)'),
        '5x4, unicode not specified': ('5x4', '', 'p(5×4)'),
        'rt5, simple': ('(rt5 x rt5)R27', 's', '(√5×√5)R27.0°'),
        'rt3, simple ascii': ('(rt3 x rt3)R30', 'sa', '(sqrt3 x sqrt3)R30.0'),
        }

    @parametrize('style,expect', _str.items(), ids=_str)
    def test_str(self, style, expect):
        """Test __str__ method."""
        woods = Woods('1x1', style=style)
        assert str(woods) == expect

    @parametrize('args,style,expect', _repr.values(), ids=_repr)
    def test_repr(self, args, style, expect):
        """Test __repr__ method."""
        woods = Woods(*args, style=style)
        assert repr(woods) == expect

    @parametrize('string,spec,expect', _fmt.values(), ids=_fmt)
    def test_format(self, string, spec, expect):
        """Test __format__ method, via f-string."""
        woods = Woods(string)
        assert f'{woods:{spec}}' == expect

