"""Tests for module math_utils of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-26'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.math_utils import angle
from viperleed.calc.lib.math_utils import cosvec
from viperleed.calc.lib.math_utils import floor_eps
from viperleed.calc.lib.math_utils import lcm


class TestAngle:
    """Tests for the angle function."""

    _valid = {
        '45deg': (([1, 0], [1, 1]), np.pi / 4),
        'opposite': (([1, 0], [-1, 0]), np.pi),
        'orthogonal': (([1, 0], [0, 1]), np.pi / 2),
        'parallel': (([1.3, 2.6], [4.2, 8.4]), 0),
        'zero': (([0, 0], [1, 7]), 0),
        'multiple': (
            ([[1, 0], [ 1, 0], [1, 0], [1.3, 2.6], [0, 0]],
             [[1, 1], [-1, 0], [0, 1], [4.2, 8.4], [1, 7]]),
            (np.pi/4, np.pi, np.pi/2, 0, 0),
            ),
        }

    @parametrize('vectors,expect', _valid.values(), ids=_valid)
    def test_valid(self, vectors, expect):
        """Check expected outcome for acceptable arguments."""
        result = angle(*vectors)
        assert result == pytest.approx(expect)

    _nargs = {
        'none': (),
        'one': (1,),
        'three': (1, 1, 1),
        }

    @parametrize(args=_nargs.values(), ids=_nargs)
    def test_raises_nargs(self, args):
        """Check complaints for too few/many arguments."""
        with pytest.raises(TypeError):
            angle(*args)


class TestCosvec:
    """Tests for the cosvec function."""

    _valid = {
        'opposite': (([1, 0], [-1, 0]), -1),
        'orthogonal': (([1, 0], [0, 1]), 0),
        'parallel': (([1.3, 2.6], [4.2, 8.4]), 1),
        'multiple': (
            ([[1, 0], [1, 0], [1.3, 2.6]],
             [[-1, 0], [0, 1], [4.2, 8.4]]),
            (-1, 0, 1),
            ),
        }
    _valid_3d = {
        'arbitrary': (([1, 2, 3], [4, 5, 6]), pytest.approx(0.97463185)),
        'identical': (([3, 4, 5], [3, 4, 5]), 1),
        'mixed components': (([1, -2, 3], [-3, 2, -1]),
                             pytest.approx(-0.7142857)),
        'near zero angle': (([2, 1, 1e-10], [2, 1, -1e-10]),
                            pytest.approx(1)),
        'near zero components': (
            ([1e-10, 1e-10, 1e-10], [1e-15, 1e-15, 1e-15]),
            pytest.approx(1)
            ),
        'negative': (([-1, -2, -3], [3, 2, 1]), pytest.approx(-0.7142857)),
        'not unit vecs': (([3, 4, 0], [0, 8, 6]), 0.64),
        'opposite': (([-1, -2, -3], [1, 2, 3]), -1),
        'parallel': (([1, 4, 7], [2, 8, 14]), 1),
        'xy': (([1, 0, 0], [0, 1, 0]), 0),
        'xz': (([1, 0, 0], [0, 0, 1]), 0),
        'yz': (([0, 1, 0], [0, 0, 1]), 0),
        }

    @parametrize('vectors,expect', _valid.values(), ids=_valid)
    def test_valid(self, vectors, expect):
        """Check expected outcome for acceptable arguments."""
        result = cosvec(*vectors)
        assert result == pytest.approx(expect)

    @parametrize('vectors,expect', _valid_3d.values(), ids=_valid_3d)
    def test_valid_3d(self, vectors, expect):
        """Check expected outcome for 3D vectors."""
        self.test_valid(vectors, expect)

    _raises = {
        'zero norm 2d': (([0, 0], [1, 1]), ValueError),
        'zero norm 3d': (([0, 0, 0], [1, 1, 1]), ValueError),
        'multiple, one zero': ((([1, 0, 0], [0, 1, 0]),
                               ([1, 0, 0], [0, 0, 0])),
                               ValueError)
        }

    @parametrize('vectors,exc', _raises.values(), ids=_raises)
    def test_raises(self, vectors, exc):
        """Check complaints for unacceptable arguments."""
        with pytest.raises(exc):
            cosvec(*vectors)


class TestFloorEps:
    """Tests for the floor_eps function."""

    def _check_floor_eps(self, eps, values, expect):
        """Check expected result of floor_eps. Return it."""
        floor_func = floor_eps(eps)
        result = floor_func(values)
        np.testing.assert_equal(result, expect)
        return result

    _valid = {  # eps, values, expected
        'eps small': (1e-8,
                      np.array([1 - 1e-8, 1 + 1e-8, 3 - 1e-8]),
                      [1, 1, 3]),
        'eps zero': (0.0, np.array([0.7, 1.5, 2.2]), [0, 1, 2]),
        'inf': (0.1, np.array([np.inf, -np.inf]), [np.inf, -np.inf]),
        'mixed': (0.4, np.array([1.5, -0.5, 2.3, -2.7]), [1, -1, 2, -3]),
        'nan': (0.1, np.array([np.nan, 1.2, np.nan]), [np.nan, 1.0, np.nan]),
        'one value': (0.7, 4.1, 4),
        'one value, array': (0.7, np.array([4.1]), [4]),
        'negative': (0.3, np.array([-1.7, -2.2, -3.1]), [-2, -2, -3]),
        'no values': (0.5, np.array([]), []),
        'positive': (0.5, np.array([1.2, 2.8, 3.3]), [1, 3, 3]),
        }

    @parametrize('eps,values,expect', _valid.values(), ids=_valid)
    def test_valid(self, eps, values, expect):
        """Check expected outcome for acceptable arguments."""
        self._check_floor_eps(eps, values, expect)

    # Highlight the difference between floor_eps and integer division
    _div = {  # eps, values, expect_floor, expect_div
        'edge': (1e-3, np.array([1 - 1e-3, -1 + 1e-3]), [1, -1], [0, -1]),
        'eps small': (1e-4, np.array([3.00005, 7.9999]), [3, 8], [3, 7]),
        'near int': (1e-4, np.array([6 - 1e-4, 4 + 1e-4]), [6, 4], [5, 4]),
        'negative': (0.05, np.array([-4.95, -5.05, -6.95]),
                     [-5, -5, -7], [-5, -6, -7]),
        'positive': (0.1, np.array([4.1, 5.9, 6.5]), [4, 6, 6], [4, 5, 6]),
        }

    @parametrize('eps,values,expect_floor,expect_div', _div.values(), ids=_div)
    def test_diff_with_integer_division(self, eps, values,
                                        expect_floor, expect_div):
        """Highlight differences between floor_eps and //."""
        floor_eps_result = self._check_floor_eps(eps, values, expect_floor)
        floor_div_result = values // 1
        np.testing.assert_array_equal(floor_div_result, expect_div)
        assert not np.allclose(floor_eps_result, floor_div_result)


class TestLCM:
    """Tests for the lcm function."""

    _valid = {
        'primes': ((5, 7), 35),
        'same': ((7, 7), 7),
        'simple': ((4, 6), 12),
        'zero first': ((0, 10), 0),
        'zero second': ((10, 0), 0),
        'multiple': (
            ((5, 7, 4, 0, 10),
             (7, 7, 6, 10, 0)),
            (35, 7, 12, 0, 0),
            )
        }

    @parametrize('values,expect', _valid.values(), ids=_valid)
    def test_valid(self, values, expect):
        """Check expected outcome for acceptable arguments."""
        result = lcm(*values)
        assert result == pytest.approx(expect)
