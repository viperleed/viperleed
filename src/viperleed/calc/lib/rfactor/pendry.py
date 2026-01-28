"""viperleed.calc.lib.rfactor.pendry."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib import dynamic_numerical_lib as dnl

from .groups import group_rfactors
from .utils import nansum_trapezoid, shift_theo_intensity_non_negative


def R_pendry(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline,
    theo_shift=0.0, **kwargs
):
    """Calculate the R-factor for two beams."""
    # Experimental data
    exp_deriv_spline = exp_spline.derivative()

    exp_intensity = exp_spline(energy_grid)
    exp_derivative = exp_deriv_spline(energy_grid)

    # Theory data
    theo_deriv_spline = theo_spline.derivative()

    shifted_grid = energy_grid - theo_shift
    theo_intensity = theo_spline(shifted_grid)
    # shift theo_intensity to be non-negative
    theo_intensity = shift_theo_intensity_non_negative(
        theo_intensity, exp_intensity
    )

    theo_derivative = theo_deriv_spline(shifted_grid)

    exp_y = y_pendry(exp_intensity, exp_derivative, v0_imag)
    theo_y = y_pendry(theo_intensity, theo_derivative, v0_imag)

    return R_pendry_from_y(exp_y, theo_y, energy_step, **kwargs)


def R_pendry_from_y(y_1, y_2, energy_step, **kwargs):
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
    intens_deriv_ratio = intensity / intensity_derivative
    return intens_deriv_ratio / (intens_deriv_ratio**2 + v0_imag**2)
