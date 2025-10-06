"""Module R-factor."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-02-21'

import jax
from jax import numpy as jnp


def pendry_R(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, per_beam=False
):
    """Calculate the R-factor for two beams."""
    # Experimental data
    exp_deriv_spline = exp_spline.derivative()

    exp_intensity = exp_spline(energy_grid)
    exp_derivative = exp_deriv_spline(energy_grid)

    # Theory data
    theo_deriv_spline = theo_spline.derivative()

    theo_intensity = theo_spline(energy_grid)
    theo_derivative = theo_deriv_spline(energy_grid)

    exp_y = pendry_y(exp_intensity, exp_derivative, v0_imag)
    theo_y = pendry_y(theo_intensity, theo_derivative, v0_imag)

    return pendry_R_from_y(exp_y, theo_y, energy_step, per_beam=per_beam)


def pendry_R_from_intensity_and_derivative(
    intens_deriv_1,
    intens_deriv_2,
    v0_real_steps,
    v0_imag,
    energy_step,
    per_beam=False,
):
    intens_1, deriv_1 = intens_deriv_1
    intens_2, deriv_2 = intens_deriv_2

    y_1 = pendry_y(intens_1, deriv_1, v0_imag)
    y_2 = pendry_y(intens_2, deriv_2, v0_imag)

    # shift y_1 by v0_real_steps
    y_1 = integer_shift_v0r(y_1, v0_real_steps)

    return pendry_R_from_y(y_1, y_2, energy_step, per_beam=per_beam)


def pendry_R_from_y(y_1, y_2, energy_step, per_beam=False):
    # mask out NaNs for this calculation
    y_1_mask = jnp.isnan(y_1)
    y_2_mask = jnp.isnan(y_2)
    mask = jnp.logical_or(y_1_mask, y_2_mask)

    y_1 = jnp.where(mask, 0, y_1)
    y_2 = jnp.where(mask, 0, y_2)

    # TODO?: potentially, one could do these integrals analytically based on the spline coefficients
    numerators = nansum_trapezoid((y_1 - y_2) ** 2, dx=energy_step, axis=0)
    denominators = nansum_trapezoid((y_1**2 + y_2**2), dx=energy_step, axis=0)
    if per_beam:
        # return R factor for each beam
        return numerators / denominators
    # R factor for all beams
    return jnp.sum(numerators) / jnp.sum(denominators)


def pendry_y(intensity, intensity_derivative, v0_imag):
    intens_deriv_ratio = intensity / intensity_derivative
    return intens_deriv_ratio / (intens_deriv_ratio**2 + v0_imag**2)


def nansum_trapezoid(y, dx, axis=-1):
    y_arr = jnp.moveaxis(y, axis, -1)
    # select the axis to integrate over
    return jnp.nansum(y_arr[..., 1:] + y_arr[..., :-1], axis=-1) * dx * 0.5


def integer_shift_v0r(array, n_steps):
    """Applies a v0r shift to the array by shifting the values n_steps up or
    down the first axis (energy) and padding with NaNs.
    """
    # NB, TODO: This only allows for integer shifts (multiples of the set
    # energy step). This is a limitation of the current implementation.
    # In principle, we could implement a more general shift and allow real
    # numbers by doing this earlier and changing the knot values in the
    # interpolator.
    n_energies, n_beams = array.shape[0], array.shape[1]

    rolled_array = jnp.roll(array, n_steps, axis=0)
    row_ids = jnp.arange(n_energies).reshape(-1, 1)
    row_ids_tiled = jnp.tile(row_ids, (1, n_beams))
    mask = jnp.logical_or(
        row_ids_tiled < n_steps, row_ids >= n_energies + n_steps
    )
    return jnp.where(mask, jnp.nan, rolled_array)


### R2 ###


def R_2(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, per_beam=False
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
    if per_beam:
        # return R factor for each beam
        return numerators / denominators
    # R factor for all beams
    return jnp.sum(numerators) / jnp.sum(denominators)


def R_1(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, per_beam=False
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
    if per_beam:
        # return R factor for each beam
        return numerators / denominators
    # R factor for all beams
    return jnp.sum(numerators) / jnp.sum(denominators)


### RMS ###


def R_ms(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, per_beam=False
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

    return pendry_R_from_y(y_exp, y_theo, energy_step, per_beam=per_beam)


def y_ms(intensity, first_derivative, second_derivative, v0_imag, e_step):
    numerator = first_derivative
    condition = second_derivative > 0
    denominator = intensity**2 + 0.5 * (first_derivative * v0_imag) ** 2
    denominator += condition * 0.1 * (second_derivative * v0_imag**2) ** 2
    denominator = jnp.sqrt(denominator)
    return numerator / denominator


def R_s(
    theo_spline,
    v0_imag,
    energy_step,
    energy_grid,
    exp_spline,
    per_beam=False,
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

    return pendry_R_from_y(y_exp, y_theo, energy_step, per_beam=per_beam)


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
    intermediate_2 = intermediate_1 / jnp.sqrt(1 + intermediate_1**2)

    numerator = first_derivative
    condition = jnp.logical_and(second_derivative > 0, intermediate_1 > 0)

    denominator = intensity**2 + 4 * v0_imag**2 * first_derivative**2
    conditional_denominator = (
        second_derivative**2 * v0_imag**4 * intermediate_2**2
    )
    denominator += condition * conditional_denominator
    denominator = jnp.sqrt(denominator)
    return numerator / denominator


def R_zj(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, per_beam=False
):
    # Experimental data
    exp_deriv_1_spline = exp_spline.derivative()
    exp_deriv_2_spline = exp_deriv_1_spline.derivative()

    exp_intensity = exp_spline(energy_grid)
    exp_derivative_1 = exp_deriv_1_spline(energy_grid)
    exp_derivative_2 = exp_deriv_2_spline(energy_grid)

    exp_mask = jnp.isnan(exp_intensity)
    exp_mask = jax.lax.stop_gradient(exp_mask)
    exp_intensity = jnp.where(exp_mask, 0.0, exp_intensity)
    exp_derivative_1 = jnp.where(exp_mask, 0.0, exp_derivative_1)
    exp_derivative_2 = jnp.where(exp_mask, 0.0, exp_derivative_2)

    # Theory data
    theo_deriv_1_spline = theo_spline.derivative()
    theo_deriv_2_spline = theo_deriv_1_spline.derivative()

    theo_intensity = theo_spline(energy_grid)
    theo_mask = jnp.isnan(theo_intensity)
    theo_mask = jax.lax.stop_gradient(theo_mask)
    mask = jnp.logical_or(exp_mask, theo_mask)

    theo_derivative_1 = theo_deriv_1_spline(energy_grid)
    theo_derivative_2 = theo_deriv_2_spline(energy_grid)

    # apply mask to theory as well
    theo_intensity = jnp.where(mask, 0.0, theo_intensity)
    theo_derivative_1 = jnp.where(mask, 0.0, theo_derivative_1)
    theo_derivative_2 = jnp.where(mask, 0.0, theo_derivative_2)

    # calculate experimental energy ranges (without NaNs)
    exp_energy_ranges = jnp.logical_not(mask).sum(axis=0) * energy_step

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
    numerators = jnp.where(mask, 0.0, numerators)

    denominators = abs(exp_derivative_1) + jnp.nanmax(exp_derivative_1, axis=0)
    denominators = jnp.clip(denominators, a_min=1e-12)
    denominators = jnp.where(mask, 1.0, denominators)

    quotient = numerators / denominators

    r_beams = prefactors * nansum_trapezoid(quotient, axis=0, dx=energy_step)
    if per_beam:
        return r_beams

    return jnp.nansum(r_beams * exp_energy_ranges) / jnp.nansum(
        exp_energy_ranges
    )


R_FACTOR_SYNONYMS = {
    pendry_R: ('pendry', 'r_p', 'rp', 'pendry r-factor'),
    R_1: ('r1', 'r_1', 'r1 factor'),
    R_2: ('r2', 'r_2', 'r2 factor'),
    R_s: ('s', 'rs', 'r_s', 'smooth'),
    R_zj: (
        'zj',
        'zj factor',
        'zannazi',
        'zannazi jona',
        'zannazi-jona',
    ),
}
