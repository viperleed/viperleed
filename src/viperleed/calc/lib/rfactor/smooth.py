"""viperleed.calc.lib.rfactor.smooth."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib import dynamic_numerical_lib as dnl

from .pendry import R_pendry_from_y
from .utils import shift_theo_intensity_non_negative


def R_ms(theo_spline, v0_imag, energy_step, energy_grid, exp_spline, **kwargs):
    # Experimental data
    exp_deriv_1_spline = exp_spline.derivative()
    exp_deriv_2_spline = exp_deriv_1_spline.derivative()

    exp_intensity = exp_spline(energy_grid)
    exp_derivative_1 = exp_deriv_1_spline(energy_grid)
    exp_derivative_2 = exp_deriv_2_spline(energy_grid)

    # Theory data
    theo_deriv_1_spline = theo_spline.derivative()
    theo_deriv_2_spline = theo_deriv_1_spline.derivative()

    theo_intensity = theo_spline(energy_grid)
    # shift theo_intensity to be non-negative
    theo_intensity = shift_theo_intensity_non_negative(
        theo_intensity, exp_intensity
    )

    theo_derivative_1 = theo_deriv_1_spline(energy_grid)
    theo_derivative_2 = theo_deriv_2_spline(energy_grid)

    y_exp = y_ms(
        exp_intensity, exp_derivative_1, exp_derivative_2, v0_imag, energy_step
    )
    y_theo = y_ms(
        theo_intensity,
        theo_derivative_1,
        theo_derivative_2,
        v0_imag,
        energy_step,
    )

    return R_pendry_from_y(y_exp, y_theo, energy_step, **kwargs)


def y_ms(intensity, first_derivative, second_derivative, v0_imag, e_step):
    numerator = first_derivative
    condition = second_derivative > 0
    denominator = intensity**2 + 0.5 * (first_derivative * v0_imag) ** 2
    denominator += condition * 0.1 * (second_derivative * v0_imag**2) ** 2
    denominator = dnl.xp.sqrt(denominator)
    return numerator / denominator


def R_s(
    theo_spline,
    v0_imag,
    energy_step,
    energy_grid,
    exp_spline,
    alpha=4.0,
    beta=0.15,
    theo_shift=0.0,
    **kwargs,
):
    # Experimental data (cacheable via CachedSpline)
    exp_deriv_1_spline = exp_spline.derivative()
    exp_deriv_2_spline = exp_deriv_1_spline.derivative()

    exp_intensity = exp_spline(energy_grid)
    exp_derivative_1 = exp_deriv_1_spline(energy_grid)
    exp_derivative_2 = exp_deriv_2_spline(energy_grid)

    y_exp = y_s(
        exp_intensity, exp_derivative_1, exp_derivative_2,
        v0_imag, energy_step, alpha=alpha, beta=beta
    )

    # Theory data
    theo_deriv_1_spline = theo_spline.derivative()
    theo_deriv_2_spline = theo_deriv_1_spline.derivative()

    shifted_grid = energy_grid - theo_shift
    theo_intensity = theo_spline(shifted_grid)
    theo_intensity = shift_theo_intensity_non_negative(theo_intensity, exp_intensity)
    theo_derivative_1 = theo_deriv_1_spline(shifted_grid)
    theo_derivative_2 = theo_deriv_2_spline(shifted_grid)

    y_theo = y_s(
        theo_intensity, theo_derivative_1, theo_derivative_2,
        v0_imag, energy_step, alpha=alpha, beta=beta
    )

    return R_pendry_from_y(y_exp, y_theo, energy_step, **kwargs)


def y_s(
    intensity,
    first_derivative,
    second_derivative,
    v0_imag,
    e_step,
    alpha=4.0,
    beta=0.15,
):
    intermediate_1 = (
        intensity / second_derivative
        - 0.5 * first_derivative**2 / second_derivative**2
    )
    intermediate_1 = (alpha / v0_imag**2) * intermediate_1 + beta
    intermediate_2 = intermediate_1 / dnl.xp.sqrt(1 + intermediate_1**2)

    numerator = first_derivative
    condition = dnl.xp.logical_and(second_derivative > 0, intermediate_1 > 0)

    denominator = intensity**2 + 4 * v0_imag**2 * first_derivative**2
    conditional_denominator = (
        second_derivative**2 * v0_imag**4 * intermediate_2**2
    )
    denominator += condition * conditional_denominator
    denominator = dnl.xp.sqrt(denominator)
    return numerator / denominator
