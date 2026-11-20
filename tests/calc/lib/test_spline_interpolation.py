"""Tests for module spline_interpolation of viperleed.calc.lib."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from viperleed.calc.lib import spline_interpolation

# Indexing the results of calls to the spline with numpy integers or
# arrays is flagged: pylint cannot infer what the spline returns
# pylint: disable=invalid-sequence-index


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
        return spline_interpolation.CachedSpline(mock_spline)

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


@pytest.fixture(name='ragged_data')
def fixture_ragged_data():
    """Return a grid and beams that start/end at different energies."""
    grid = np.arange(10.0, 20.0)
    beams = np.full((len(grid), 3), np.nan)
    first_last = ((0, 9), (0, 5), (3, 7))  # Inclusive indices of the data
    for col, (first, last) in enumerate(first_last):
        beams[first:last + 1, col] = np.sin(grid[first:last + 1] / 3 + col)
    return grid, beams


class TestInterpolateRaggedArray:
    """Tests for interpolate_ragged_array."""

    def test_last_knot_is_data_point(self, ragged_data):
        """Check that each beam's last data point is not lost."""
        grid, beams = ragged_data
        spline = spline_interpolation.interpolate_ragged_array(grid, beams)
        values = spline(grid)
        for col in range(beams.shape[1]):
            last = np.flatnonzero(~np.isnan(beams[:, col]))[-1]
            assert values[last, col] == pytest.approx(beams[last, col])

    def test_nan_outside_beam_range(self, ragged_data):
        """Check that beams are NaN outside their own data range."""
        grid, beams = ragged_data
        spline = spline_interpolation.interpolate_ragged_array(grid, beams)
        eps = 0.01
        for col in range(beams.shape[1]):
            energies = grid[~np.isnan(beams[:, col])]
            below, above = energies[0] - eps, energies[-1] + eps
            for deriv in (spline, spline.derivative(nu=1)):
                assert np.isnan(deriv(np.array([below, above]))[:, col]).all()
                inside = deriv(np.array([energies[0], energies[-1]]))
                assert np.isfinite(inside[:, col]).all()

    @pytest.mark.parametrize('order', [1, 2])
    def test_derivatives_at_last_knot(self, ragged_data, order):
        """Check derivatives at the last knot match a single spline."""
        grid, beams = ragged_data
        spline = spline_interpolation.interpolate_ragged_array(grid, beams)
        derivatives = spline.derivative(nu=order)
        for col in range(beams.shape[1]):
            has_data = ~np.isnan(beams[:, col])
            single = CubicSpline(grid[has_data], beams[has_data, col],
                                 bc_type='not-a-knot')
            last = grid[has_data][-1:]
            expected = single(last, nu=order)[0]
            assert derivatives(last)[0, col] == pytest.approx(expected)


_ENERGY_TOL = 1e-6


@pytest.mark.parametrize('step', [0.1, 0.2, 1/3])
@pytest.mark.parametrize('shift', [0, 0.3, -1.7])
def test_computed_energies_keep_beam_ends(step, shift):
    """Check beam ends are kept for computed (rounded) energies."""
    knots = np.round(27.3 + np.arange(400) * step, 10)
    beams = np.full((len(knots), 2), np.nan)
    beams[:, 0] = np.sin(knots / 20)
    beams[5:-7, 1] = np.cos(knots[5:-7] / 20)
    spline = spline_interpolation.interpolate_ragged_array(knots, beams)
    # Energies computed as in the V0r scan: arange grid, minus a shift
    grid = np.arange(knots[0] + shift - 1, knots[-1] + shift + 1, step)
    energies = grid - shift
    for col, (first, last) in enumerate(((0, -1), (5, -8))):
        for knot in (first, last):
            idx = np.argmin(np.abs(energies - knots[knot]))
            assert abs(energies[idx] - knots[knot]) < _ENERGY_TOL
            value = spline(energies[idx:idx + 1])[0, col]
            assert value == pytest.approx(beams[knot, col])
            deriv = spline.derivative(nu=1)(energies[idx:idx + 1])[0, col]
            assert np.isfinite(deriv)
