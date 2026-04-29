"""viperleed.calc.lib.rfactor.zannazi_jona."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib import dynamic_numerical_lib as dnl

from .utils import nansum_trapezoid
from .groups import group_rfactors


def R_zj(
    v0_imag,
    energy_step,
    energy_grid,
    data_spline_1=None,
    data_and_derivatives_1=None,
    data_spline_2=None,
    data_and_derivatives_2=None,
    shift_2nd_spline=0.0,   # only available if passed as spline
    **kwargs,
    ):
    """
    Zannazi-Jona R-factor, see
    https://www.sciencedirect.com/science/article/pii/0039602877904289
    Uses two sets of beam data, either passed as splines or as already
    calculated arrays of intensities and derivatives.

    For each dataset, pass either data_spline or data_and_derivatives, but not
    both. Using splines for one and pre-computed data for the other is allowed.

    Parameters
    ----------
    v0_imag : float
        The imaginary part of the inner potential.
    energy_step : float
        The energy step of the data grid.
    energy_grid : array
        The grid on which data should be evaluated. All evaluations will ignore
        regions in which either of the datasets is nan.
    data_spline_1, data_spline_2 : arrays of splines
        Splines of the data, which will be evaluated on energy_grid. Need to be
        at least cubic for calculating derivatives. Ignored if the respective
        data_and_derivatives are passed.
    data_and_derivatives_1, data_and_derivatives_2 : 3-tuples of float arrays
        Pre-evaluated tuples (intensity, 1st derivative, 2nd derivative) at the
        energy_grid points. This avoids re-calculation of splines and is more
        JAX-friendly.
    shift_2nd_spline : float, optional
        Evaluate the 2nd spline on a shifted energy grid. This is meant for
        testing V0r variations. Note that this is NOT available when the 2nd
        dataset is passed as data_and_derivatives.
    """
    # Get data either as splines or as pre-computed arrays (mainly for JAX)
    if data_and_derivatives_1 is None:
        if data_spline_1 is None:
            raise TypeError('R_zj requires either data splines or pre-computed'
                            'data_and_derivatives arrays.')
        # when using splines, this can be sped up via CashedSplines
        data_1_deriv_1_spline = data_spline_1.derivative()
        data_1_deriv_2_spline = data_1_deriv_1_spline.derivative()
        data_1_intensity = data_spline_1(energy_grid)
        data_1_derivative_1 = data_1_deriv_1_spline(energy_grid)
        data_1_derivative_2 = data_1_deriv_2_spline(energy_grid)
    else:
        data_1_intensity, data_1_derivative_1, data_1_derivative_2 = (
            data_and_derivatives_1
            )
    # 2nd data also uses shifted grid if spline
    shifted_grid = energy_grid - shift_2nd_spline
    if data_and_derivatives_2 is None:
        if data_spline_2 is None:
            raise TypeError('R_zj requires either data splines or pre-computed'
                            'data_and_derivatives arrays.')
        # when using splines, this can be sped up via CashedSplines
        data_2_deriv_1_spline = data_spline_2.derivative()
        data_2_deriv_2_spline = data_2_deriv_1_spline.derivative()
        # evaluate on shifted grid, allowing continuous shifts
        data_2_intensity = data_spline_2(shifted_grid)
        data_2_derivative_1 = data_2_deriv_1_spline(shifted_grid)
        data_2_derivative_2 = data_2_deriv_2_spline(shifted_grid)
    else:
        data_2_intensity, data_2_derivative_1, data_2_derivative_2 = (
            data_and_derivatives_2
            )
    mask_1 = dnl.xp.isnan(data_1_intensity)
    data_1_intensity = dnl.xp.where(mask_1, 0.0, data_1_intensity)
    data_1_derivative_1 = dnl.xp.where(mask_1, 0.0, data_1_derivative_1)
    data_1_derivative_2 = dnl.xp.where(mask_1, 0.0, data_1_derivative_2)

    mask_2 = dnl.xp.isnan(data_2_intensity)
    mask = dnl.xp.logical_or(mask_1, mask_2)
    data_2_intensity = dnl.xp.where(mask, 0.0, data_2_intensity)
    data_2_derivative_1 = dnl.xp.where(mask, 0.0, data_2_derivative_1)
    data_2_derivative_2 = dnl.xp.where(mask, 0.0, data_2_derivative_2)

    # calculate experimental energy ranges (without NaNs)
    exp_energy_ranges = dnl.xp.logical_not(mask).sum(axis=0) * energy_step

    # Factor 0.027 for random correlation, Zannazi & Jona 1977
    prefactors = (
        1 / nansum_trapezoid(data_1_intensity, energy_step, axis=0) / 0.027
    )

    # calculate normalization for each beam
    beam_normalization = nansum_trapezoid(
        data_1_intensity, dx=energy_step, axis=0
    ) / nansum_trapezoid(data_2_intensity, dx=energy_step, axis=0)

    numerators = abs(
        beam_normalization * data_2_derivative_2 - data_1_derivative_2
    ) * abs(beam_normalization * data_2_derivative_1 - data_1_derivative_1)
    numerators = dnl.xp.where(mask, 0.0, numerators)

    denominators = abs(data_1_derivative_1) + dnl.xp.nanmax(
        data_1_derivative_1, axis=0
    )
    denominators = dnl.xp.clip(denominators, 1e-12, None)
    denominators = dnl.xp.where(mask, 1.0, denominators)

    quotient = numerators / denominators

    # calculate R-factor with requested grouping
    return group_rfactors(
        prefactors
        * nansum_trapezoid(quotient, axis=0, dx=energy_step)
        * exp_energy_ranges,
        exp_energy_ranges,
        **kwargs,
    )
