"""viperleed.calc.lib.rfactor.smooth."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib.dynamic_numerical_lib import xp
from .utils import shift_theo_intensity_non_negative, nansum_trapezoid
from .groups import group_rfactors


def R_zj(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, groups=None
):
    # Experimental data
    exp_deriv_1_spline = exp_spline.derivative()
    exp_deriv_2_spline = exp_deriv_1_spline.derivative()

    exp_intensity = exp_spline(energy_grid)
    exp_derivative_1 = exp_deriv_1_spline(energy_grid)
    exp_derivative_2 = exp_deriv_2_spline(energy_grid)

    exp_mask = xp.isnan(exp_intensity)
    exp_intensity = xp.where(exp_mask, 0.0, exp_intensity)
    exp_derivative_1 = xp.where(exp_mask, 0.0, exp_derivative_1)
    exp_derivative_2 = xp.where(exp_mask, 0.0, exp_derivative_2)

    # Theory data
    theo_deriv_1_spline = theo_spline.derivative()
    theo_deriv_2_spline = theo_deriv_1_spline.derivative()

    theo_intensity = theo_spline(energy_grid)
    theo_mask = xp.isnan(theo_intensity)
    mask = xp.logical_or(exp_mask, theo_mask)

    theo_derivative_1 = theo_deriv_1_spline(energy_grid)
    theo_derivative_2 = theo_deriv_2_spline(energy_grid)

    # apply mask to theory as well
    theo_intensity = xp.where(mask, 0.0, theo_intensity)
    theo_derivative_1 = xp.where(mask, 0.0, theo_derivative_1)
    theo_derivative_2 = xp.where(mask, 0.0, theo_derivative_2)

    # calculate experimental energy ranges (without NaNs)
    exp_energy_ranges = xp.logical_not(mask).sum(axis=0) * energy_step

    # Factor 0.027 for random correlation, Zannazi & Jona 1977
    prefactors = (
        1 / nansum_trapezoid(exp_intensity, energy_step, axis=0) / 0.027
    )

    # # calculate normalization for each beam
    beam_normalization = nansum_trapezoid(
        exp_intensity, dx=energy_step, axis=0
    ) / nansum_trapezoid(theo_intensity, dx=energy_step, axis=0)

    numerators = abs(
        beam_normalization * theo_derivative_2 - exp_derivative_2
    ) * abs(beam_normalization * theo_derivative_1 - exp_derivative_1)
    numerators = xp.where(mask, 0.0, numerators)

    denominators = abs(exp_derivative_1) + xp.nanmax(exp_derivative_1, axis=0)
    denominators = xp.clip(denominators, a_min=1e-12)
    denominators = xp.where(mask, 1.0, denominators)

    quotient = numerators / denominators

    # calculate R-factor with requested grouping
    return group_rfactors(
        prefactors
        * nansum_trapezoid(quotient, axis=0, dx=energy_step)
        * exp_energy_ranges,
        exp_energy_ranges,
        groups=groups,
    )
