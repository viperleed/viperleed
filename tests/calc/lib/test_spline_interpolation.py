"""Tests for module spline_interpolation of viperleed.calc.lib."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from viperleed.calc.lib.spline_interpolation import CachedSpline


class TestCachedSpline:
    """Tests for the CachedSpline proxy."""

    # pylint: disable=protected-access                # OK in tests

    @pytest.fixture(name='spline')
    def fixture_spline(self, mocker):
        """Return a CachedSpline wrapping a mock of a cubic spline."""
        x_data = np.linspace(0, 10, 11)
        return CachedSpline(
            mocker.MagicMock(wraps=CubicSpline(x_data, np.sin(x_data)))
            )

    def test_same_array_cached(self, spline):
        """Check that evaluating the same array twice is cached."""
        grid = np.linspace(1, 9, 50)
        first = spline(grid)
        second = spline(grid)
        assert second is first
        assert spline._spline.call_count == 1

    def test_different_array_recomputed(self, spline):
        """Check that a different array is evaluated anew."""
        grid = np.linspace(1, 9, 50)
        other_grid = grid + 0.5
        first = spline(grid)
        second = spline(other_grid)
        # pylint: disable-next=magic-value-comparison
        assert spline._spline.call_count == 2
        assert not np.allclose(second, first)

    def test_equal_copy_recomputed(self, spline):
        """Check that an equal but distinct array is evaluated anew."""
        grid = np.linspace(1, 9, 50)
        first = spline(grid)
        second = spline(grid.copy())
        # pylint: disable-next=magic-value-comparison
        assert spline._spline.call_count == 2
        np.testing.assert_array_equal(second, first)

    def test_derivative_cached(self, spline):
        """Check that derivatives are created only once per order."""
        assert spline.derivative(nu=1) is spline.derivative(nu=1)
        assert spline.derivative(nu=1) is not spline.derivative(nu=2)
