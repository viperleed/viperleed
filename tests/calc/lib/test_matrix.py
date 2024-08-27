"""Tests for module matrix of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-27'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.matrix import NonIntegerMatrixError
from viperleed.calc.lib.matrix import SingularMatrixError
from viperleed.calc.lib.matrix import ensure_integer_matrix
from viperleed.calc.lib.matrix import rotation_matrix
from viperleed.calc.lib.matrix import rotation_matrix_order


class TestEnsureIntegerMatrix:
    """Tests for the ensure_integer_matrix function."""

    _valid = {
        'int': ([[1, 2], [3, 4]], [[1, 2], [3, 4]]),
        'almost eps away': (
            [[1 + 9.999e-7, 3 - 9.999e-7], [3 + 9.999e-7, 4]],
            [[1, 3], [3, 4]],
            ),
        'int as float': ([[1.0, 2.0], [3.0, 4.0]], [[1, 2], [3, 4]]),
        'less than eps': (
            [[1 + 5e-7, 2 - 5e-7], [3 + 5e-7, 4 + 8e-7]],
            [[1, 2], [3, 4]],
            ),
        'negative': (
            [[-1 - 1e-7, -2 + 1e-7], [-3 - 9.999e-7, -5 + 7e-7]],
            [[-1, -2], [-3, -5]],
            ),
        'some zero': (
            [[1e-7, 1 - 1e-7], [1e-7, -1 + 1e-7]],
            [[0, 1], [0, -1]],
            ),
        'zero': (np.zeros((3, 3)), np.zeros((3, 3))),
        }

    @parametrize('matrix,expect', _valid.values(), ids=_valid)
    def test_valid(self, matrix, expect):
        """Check expected outcome for acceptable arguments."""
        result = ensure_integer_matrix(matrix)
        np.testing.assert_array_equal(result, expect)

    def test_large_matrix(self):
        """Check correct result for a large matrix."""
        rounded_matrix = np.random.rand(1000, 1000).round()
        result = ensure_integer_matrix(rounded_matrix + 1e-6)
        np.testing.assert_array_equal(result, rounded_matrix)

    _invalid = {
        'large': (np.random.rand(1000, 1000),),
        'mixed inf nan': ([[np.nan, np.inf], [-np.inf, 4]],),
        'not integer': ([[1.1, 2.5], [3.7, 4.2]],),
        'nan': ([[1, 2], [np.nan, 4]],),
        'nan, small eps': ([[1, 2], [np.nan, 4]], 1e-12),
        '+inf': ([[1, 2], [np.inf, 4]],),
        '+inf, small eps': ([[1, 2], [np.inf, 4]],),
        '-inf': ([[1, 2], [-np.inf, 4]],),
        'large eps': ([[1.1, 2.2], [3.3, 4.4]], 0.3),
        'small eps': ([[1 + 1e-6, 2 - 1e-7], [3, 4]], 1e-12),
        }

    @parametrize(args=_invalid.values(), ids=_invalid)
    def test_invalid(self, args):
        """Check complaints for non-integer matrices."""
        with pytest.raises(NonIntegerMatrixError):
            ensure_integer_matrix(*args)


class TestRotationMatrix:
    """Tests for the rotation_matrix function."""

    _valid = {
        '2d': ((np.pi / 4,), [[2**-0.5, -2**-0.5], [2**-0.5, 2**-0.5]]),
        '3d': ((np.pi / 2, 3), [[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
        '4d': ((np.pi / 3, 4), [
            [0.5, -np.sqrt(3)/2, 0, 0],
            [np.sqrt(3)/2, 0.5, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
            ]),
        'large, negative': ((-15 * np.pi - np.pi / 6,),
                            [[-np.sqrt(3)/2, -0.5], [0.5, -np.sqrt(3)/2]]),
        'large, positive': ((20 * np.pi + np.pi / 4,),
                            [[2**-0.5, -2**-0.5], [2**-0.5, 2**-0.5]]),
        'negative': ((-np.pi / 4,), [[2**-0.5, 2**-0.5], [-2**-0.5, 2**-0.5]]),
        'zero': ((0, 3), np.eye(3)),
        'zero, large': ((20 * np.pi,), np.eye(2)),
        'zero, large, negative': ((-10 * np.pi,), np.eye(2)),
        }

    @parametrize('args,expect', _valid.values(), ids=_valid)
    def test_valid(self, args, expect):
        """Check expected outcome for acceptable arguments."""
        matrix = rotation_matrix(*args)
        np.testing.assert_allclose(matrix, expect, atol=1e-10)

    @parametrize(ndims=(-3, -1, 0, 1))
    def test_raises(self, ndims):
        """Check complaints for invalid number of dimensions."""
        with pytest.raises(ValueError):
            rotation_matrix(np.pi / 4, dim=ndims)


class TestRotationMatrixOrder:
    """Tests for the rotation_matrix_order function."""

    _valid = {
        '2d': ((4,), [[0, -1], [1, 0]]),
        '3d': ((4, 3), [[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
        'identity': ((1, 3), np.eye(3)),
        'large': ((1_000_000,),
                  [[1, -np.pi/500_000], [np.pi/500_000, 1]]),
        'large, neg': ((-1_000_000,),
                       [[1, np.pi/500_000], [-np.pi/500_000, 1]]),
        'negative': ((-4, 2), [[0, 1],[-1, 0]]),
        }

    @parametrize('args,expect', _valid.values(), ids=_valid)
    def test_valid(self, args, expect):
        """Check expected outcome for acceptable arguments."""
        matrix = rotation_matrix_order(*args)
        np.testing.assert_allclose(matrix, expect, atol=1e-10)

    _float_order = {
        'small': (2.5, 2),
        'large': (1000.2, 3),
        }

    @parametrize('order,ndims', _float_order.values(), ids=_float_order)
    def test_float_order(self, order, ndims):
        """Check correct behavior for a non-integer rotation order."""
        matrix = rotation_matrix_order(order, dim=ndims)
        expect = rotation_matrix(2 * np.pi / order, dim=ndims)
        np.testing.assert_allclose(matrix, expect, atol=1e-10)

    @parametrize(ndims=(-3, -1, 0, 1))
    def test_raises(self, ndims):
        """Check complaints for invalid number of dimensions."""
        with pytest.raises(ValueError):
            rotation_matrix_order(4, dim=ndims)
