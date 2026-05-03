"""Module pendry of viperleed.calc.lib.rfactor."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib import dynamic_numerical_lib as dnl

from .groups import group_rfactors
from .utils import (
    nansum_trapezoid,
    shift_theo_intensity_non_negative,
    coerce_to_energy_grid,
)


def r_pendry(
    v0_imag,
    energy_step,
    energy_grid,
    data_spline_1=None,
    data_and_derivatives_1=None,
    data_spline_2=None,
    data_and_derivatives_2=None,
    shift_2nd_spline=0.0,  # only available if passed as spline
    **kwargs,
):
    """Calculate Pendry's R factor between I(V) curves.

    See https://doi.org/10.1088/0022-3719/13/5/024
    Uses two sets of beam data, either passed as splines or as already
    calculated arrays of intensities and derivatives. When comparing
    experimental and theoretical data, the theoretical data should be the
    second dataset, as this is also shifted to correct for non-positive values.
    For each dataset, pass either data_spline or data_and_derivatives, but not
    both. Using splines for one and pre-computed data for the other is allowed.
    Note that data_and_derivatives are expected to be 3-tuples for uniform call
    signatures with R-factors that also use the 2nd derivative; however, the
    2nd derivative is not used and can be left empty.

    Parameters
    ----------
    v0_imag : float
        The imaginary part of the inner potential.
    energy_step : float
        The energy step of the data grid.
    energy_grid : array
        The grid on which data should be evaluated. All evaluations will ignore
        regions in which either of the datasets is nan.
    data_spline_1, data_spline_2 : scipy.interpolate.CubicSpline or similar
        Splines of the data, which will be evaluated on energy_grid.
        Ignored if the respective data_and_derivatives are passed.
    data_and_derivatives_1, data_and_derivatives_2 : 3-tuples of float arrays
        Pre-evaluated tuples (intensity, 1st derivative, 2nd derivative) at the
        energy_grid points. This avoids re-calculation of splines and is more
        JAX-friendly. The 2nd derivative is not actually used, so passing
        (intensity, 1st derivative, None) would be acceptable.
    shift_2nd_spline : float, optional
        Evaluate the 2nd spline on a shifted energy grid. This is meant for
        testing V0r variations. Note that this is NOT available when the 2nd
        dataset is passed as data_and_derivatives.
    **kwargs
        Additional keyword arguments are passed to group_rfactors, which
        determines the grouping of beams for the final R-factor
        calculation.

    Raises
    ------
    TypeError
        If neither data_spline nor data_and_derivatives is passed for
        either dataset, or if both are passed for the same dataset.
    """
    # Get data either as splines or as pre-computed arrays (mainly for JAX)
    data_1_intensity, data_1_derivative = coerce_to_energy_grid(
        data_and_derivatives_1, data_spline_1, energy_grid, derivs=1
    )

    # 2nd data uses shifted grid if spline
    shifted_grid = energy_grid - shift_2nd_spline
    data_2_intensity, data_2_derivative = coerce_to_energy_grid(
        data_and_derivatives_2, data_spline_2, shifted_grid, derivs=1
    )

    y_1 = y_pendry(data_1_intensity, data_1_derivative, v0_imag)
    y_2 = y_pendry(data_2_intensity, data_2_derivative, v0_imag)

    return r_pendry_from_y(y_1, y_2, energy_step, **kwargs)


def r_pendry_from_y(y_1, y_2, energy_step, **kwargs):
    """Calculate Pendry's R factor from pre-computed Y function values.

    Parameters
    ----------
    y_1, y_2 : ndarray, shape (n_beams, n_points)
        The Y function values for the two datasets, evaluated on the
        same energy grid. n_beams is the number of diffraction beams,
        and n_points is the number of energy points.
        Arrays may contain NaN values marking points that
        should be ignored in the calculation.
    energy_step : float
        The energy step of the data grid.
    **kwargs
        Additional keyword arguments are passed to group_rfactors.

    Returns
    -------
    r_factors : ndarray, shape (n_groups,)
        Overall or per-group R factors, depending on the grouping
        requested.
    """
    # mask out NaNs for this calculation
    y_1_mask = dnl.xp.isnan(y_1)
    y_2_mask = dnl.xp.isnan(y_2)
    mask = dnl.xp.logical_or(y_1_mask, y_2_mask)

    y_1 = dnl.xp.where(mask, 0, y_1)
    y_2 = dnl.xp.where(mask, 0, y_2)

    numerators = nansum_trapezoid((y_1 - y_2) ** 2, dx=energy_step, axis=0)
    denominators = nansum_trapezoid((y_1**2 + y_2**2), dx=energy_step, axis=0)

    # calculate R-factor with requested grouping
    return group_rfactors(numerators, denominators, **kwargs)


def y_pendry(intensity, intensity_derivative, v0_imag):
    """Calculate the Y function used in Pendry's R factor.

    Parameters
    ----------
    intensity : array
        The intensity I(V) values of the beam.
    intensity_derivative : array
        The first derivative dI/dV of the intensity values.
    v0_imag : float
        The imaginary part of the inner potential V0.

    Returns
    -------
    array
        The Y function values, calculated as
        I * (dI/dV) / (I^2 + (V0_imag * dI/dV)^2).
    """
    numerator = intensity * intensity_derivative
    denominator = intensity**2 + v0_imag**2 * intensity_derivative**2
    return numerator / denominator
