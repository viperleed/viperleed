"""Module R-factor."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-02-21'
__license__ = 'GPLv3+'


import numpy as _np
xp = _np


def pendry_R(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, groups=None
):
    """Calculate the R-factor for two beams."""
    # Experimental data
    exp_deriv_spline = exp_spline.derivative()

    exp_intensity = exp_spline(energy_grid)
    exp_derivative = exp_deriv_spline(energy_grid)

    # Theory data
    theo_deriv_spline = theo_spline.derivative()

    theo_intensity = theo_spline(energy_grid)
    # shift theo_intensity to be non-negative
    theo_intensity = _shift_theo_intensity_non_negative(
        theo_intensity, exp_intensity
    )

    theo_derivative = theo_deriv_spline(energy_grid)

    exp_y = pendry_y(exp_intensity, exp_derivative, v0_imag)
    theo_y = pendry_y(theo_intensity, theo_derivative, v0_imag)

    return pendry_R_from_y(exp_y, theo_y, energy_step, groups=groups)



def pendry_R_from_y(y_1, y_2, energy_step, groups=None):
    # mask out NaNs for this calculation
    y_1_mask = xp.isnan(y_1)
    y_2_mask = xp.isnan(y_2)
    mask = xp.logical_or(y_1_mask, y_2_mask)

    y_1 = xp.where(mask, 0, y_1)
    y_2 = xp.where(mask, 0, y_2)

    # TODO?: potentially, one could do these integrals analytically based on the spline coefficients
    numerators = nansum_trapezoid((y_1 - y_2) ** 2, dx=energy_step, axis=0)
    denominators = nansum_trapezoid((y_1**2 + y_2**2), dx=energy_step, axis=0)

    # calculate R-factor with requested grouping
    return group_rfactors(numerators, denominators, groups=groups)


def pendry_y(intensity, intensity_derivative, v0_imag):
    intens_deriv_ratio = intensity / intensity_derivative
    return intens_deriv_ratio / (intens_deriv_ratio**2 + v0_imag**2)




### R2 ###


def R_2(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, groups=None
):
    # calculate interpolation only – no derivatives needed for R2

    # Experimental data
    exp_intensity = exp_spline(energy_grid)

    # Theory data
    theo_intensity = theo_spline(energy_grid)

    # calculate normalization for each beam
    beam_normalization = nansum_trapezoid(
        exp_intensity, energy_step, axis=0
    ) / nansum_trapezoid(theo_intensity, energy_step, axis=0)

    numerators = nansum_trapezoid(
        (exp_intensity - beam_normalization * theo_intensity) ** 2,
        energy_step,
        axis=0,
    )
    denominators = nansum_trapezoid(exp_intensity**2, energy_step, axis=0)

    # calculate R-factor with requested grouping
    return group_rfactors(numerators, denominators, groups=groups)


def R_1(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, groups=None
):
    # Experimental data
    exp_deriv_spline = exp_spline.derivative()

    exp_intensity = exp_spline(energy_grid)

    # Theory data
    theo_deriv_spline = theo_spline.derivative()

    theo_intensity = theo_spline(energy_grid)

    # calculate normalization for each beam
    beam_normalization = nansum_trapezoid(
        exp_intensity, energy_step, axis=0
    ) / nansum_trapezoid(theo_intensity, energy_step, axis=0)

    numerators = nansum_trapezoid(
        abs(exp_intensity - beam_normalization * theo_intensity),
        energy_step,
        axis=0,
    )
    denominators = nansum_trapezoid(exp_intensity, energy_step, axis=0)

    # calculate R-factor with requested grouping
    return group_rfactors(numerators, denominators, groups=groups)


### RMS ###


def R_ms(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, groups=None
):
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
    theo_intensity = _shift_theo_intensity_non_negative(
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

    return pendry_R_from_y(y_exp, y_theo, energy_step, groups=groups)


def y_ms(intensity, first_derivative, second_derivative, v0_imag, e_step):
    numerator = first_derivative
    condition = second_derivative > 0
    denominator = intensity**2 + 0.5 * (first_derivative * v0_imag) ** 2
    denominator += condition * 0.1 * (second_derivative * v0_imag**2) ** 2
    denominator = xp.sqrt(denominator)
    return numerator / denominator


def R_s(
    theo_spline,
    v0_imag,
    energy_step,
    energy_grid,
    exp_spline,
    groups=None,
    alpha=4.0,
    beta=0.15,
):
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
    theo_intensity = _shift_theo_intensity_non_negative(
        theo_intensity, exp_intensity
    )

    theo_derivative_1 = theo_deriv_1_spline(energy_grid)
    theo_derivative_2 = theo_deriv_2_spline(energy_grid)

    y_exp = y_s(
        exp_intensity, exp_derivative_1, exp_derivative_2, v0_imag, energy_step
    )
    y_theo = y_s(
        theo_intensity,
        theo_derivative_1,
        theo_derivative_2,
        v0_imag,
        energy_step,
        alpha=alpha,
        beta=beta,
    )

    return pendry_R_from_y(y_exp, y_theo, energy_step, groups=groups)


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
    intermediate_2 = intermediate_1 / xp.sqrt(1 + intermediate_1**2)

    numerator = first_derivative
    condition = xp.logical_and(second_derivative > 0, intermediate_1 > 0)

    denominator = intensity**2 + 4 * v0_imag**2 * first_derivative**2
    conditional_denominator = (
        second_derivative**2 * v0_imag**4 * intermediate_2**2
    )
    denominator += condition * conditional_denominator
    denominator = xp.sqrt(denominator)
    return numerator / denominator


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


