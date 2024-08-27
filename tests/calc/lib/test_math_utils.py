"""Tests for module math_utils of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-26'
__license__ = 'GPLv3+'

import numpy as np
from pytest_cases import parametrize

from viperleed.calc.lib.math_utils import floor_eps


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
