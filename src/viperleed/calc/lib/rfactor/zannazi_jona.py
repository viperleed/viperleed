"""viperleed.calc.lib.rfactor.zannazi_jona."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib import dynamic_numerical_lib as dnl

from .groups import group_rfactors, overlap_length_terms
from .utils import coerce_to_energy_grid, overlap_trapezoid


def r_zjj(
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
    """
    Zannazi-Jona R-factor, see
    https://doi.org/10.1016/0039-6028(77)90428-9.
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
    data_1_intensity, data_1_derivative_1, data_1_derivative_2 = (
        coerce_to_energy_grid(
            data_and_derivatives_1, data_spline_1, energy_grid, derivs=2
        )
    )
    # 2nd data also uses shifted grid if spline
    shifted_grid = energy_grid - shift_2nd_spline
    data_2_intensity, data_2_derivative_1, data_2_derivative_2 = (
        coerce_to_energy_grid(
            data_and_derivatives_2, data_spline_2, shifted_grid, derivs=2
        )
    )
    nan_mask_data_1 = dnl.xp.isnan(data_1_intensity)
    nan_mask_data_2 = dnl.xp.isnan(data_2_intensity)
    no_overlap_mask = dnl.xp.logical_or(nan_mask_data_1, nan_mask_data_2)

    # set all values to 0 in regions where either dataset is nan
    # to avoid issues in the calculation

    data_1_intensity = dnl.xp.where(no_overlap_mask, 0.0, data_1_intensity)
    data_1_derivative_1 = dnl.xp.where(
        no_overlap_mask, 0.0, data_1_derivative_1
    )
    data_1_derivative_2 = dnl.xp.where(
        no_overlap_mask, 0.0, data_1_derivative_2
    )
    data_2_intensity = dnl.xp.where(no_overlap_mask, 0.0, data_2_intensity)
    data_2_derivative_1 = dnl.xp.where(
        no_overlap_mask, 0.0, data_2_derivative_1
    )
    data_2_derivative_2 = dnl.xp.where(
        no_overlap_mask, 0.0, data_2_derivative_2
    )

    overlap_mask = dnl.xp.logical_not(no_overlap_mask)

    # Factor 0.027 for random correlation, Zannazi & Jona 1977
    exp_integral = overlap_trapezoid(
        data_1_intensity, energy_step, overlap_mask, axis=0
    )
    # Array comparisons: bool(array) would be ambiguous
    # pylint: disable=compare-to-zero
    # (nan for beams without overlap)
    prefactors = 1 / dnl.xp.where(
        exp_integral != 0, exp_integral, dnl.xp.nan
    ) / 0.027

    # calculate normalization for each beam
    theo_integral = overlap_trapezoid(
        data_2_intensity, energy_step, overlap_mask, axis=0
    )
    beam_normalization = exp_integral / dnl.xp.where(
        theo_integral != 0, theo_integral, dnl.xp.nan
    )

    numerators = abs(
        beam_normalization * data_2_derivative_2 - data_1_derivative_2
    ) * abs(beam_normalization * data_2_derivative_1 - data_1_derivative_1)
    numerators = dnl.xp.where(no_overlap_mask, 0.0, numerators)

    # epsilon: largest experimental |I'| in the overlap
    # (data_1_derivative_1 is zero outside it)
    denominators = abs(data_1_derivative_1) + dnl.xp.max(
        abs(data_1_derivative_1), axis=0
    )
    denominators = dnl.xp.clip(denominators, 1e-12, None)
    denominators = dnl.xp.where(no_overlap_mask, 1.0, denominators)

    quotient = numerators / denominators

    # Per-beam values are the numerators; groups and the overall
    # value are the means over beams weighted by the length of the
    # overlap.
    beam_rfactors = prefactors * overlap_trapezoid(
        quotient, energy_step, overlap_mask, axis=0
    )
    numerators, denominators = overlap_length_terms(
        beam_rfactors, overlap_mask, energy_step
    )
    return group_rfactors(numerators, denominators, **kwargs)
