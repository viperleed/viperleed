"""Tests for module spline_interpolation of viperleed.calc.lib."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from viperleed.calc.lib.spline_interpolation import CachedSpline


class TestCachedSpline:
    """Tests for the CachedSpline proxy."""

    @pytest.fixture(name='mock_spline')
    def fixture_mock_spline(self, mocker):
        """Return a mock wrapping a cubic spline."""
        x_data = np.linspace(0, 10, 11)
        return mocker.MagicMock(wraps=CubicSpline(x_data, np.sin(x_data)))

    @pytest.fixture(name='spline')
    def fixture_spline(self, mock_spline):
        """Return a CachedSpline wrapping `mock_spline`."""
        return CachedSpline(mock_spline)

    def test_same_array_cached(self, spline, mock_spline):
        """Check that evaluating the same array twice is cached."""
        grid = np.linspace(1, 9, 50)
        first = spline(grid)
        second = spline(grid)
        assert second is first
        mock_spline.assert_called_once()

    @pytest.mark.parametrize('make_other', [
        pytest.param(lambda grid: grid + 0.5, id='different'),
        pytest.param(np.copy, id='equal copy'),
        ])
    def test_other_array_recomputed(self, spline, mock_spline, make_other):
        """Check that arrays other than the cached one are evaluated."""
        grid = np.linspace(1, 9, 50)
        other_grid = make_other(grid)
        spline(grid)
        result = spline(other_grid)
        assert mock_spline.call_args.args[0] is other_grid
        np.testing.assert_array_equal(result, mock_spline(other_grid))

    def test_derivative_cached(self, spline):
        """Check that derivatives are created only once per order."""
        assert spline.derivative(nu=1) is spline.derivative(nu=1)
        assert spline.derivative(nu=1) is not spline.derivative(nu=2)
