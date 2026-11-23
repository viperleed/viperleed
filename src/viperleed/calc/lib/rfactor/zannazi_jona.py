"""viperleed.calc.lib.rfactor.zannazi_jona."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib import dynamic_numerical_lib as dnl

from . import utils
from .groups import group_beam_rfactors


def _zj_quotient(beam_normalization, exp_derivs, theo_derivs, no_overlap):
    """Return the ZJ integrand |num| / eps at each point.

    Parameters
    ----------
    beam_normalization : dnl.xp.ndarray
        Factor scaling the theoretical to the experimental intensity,
        one value per beam.
    exp_derivs, theo_derivs : tuple of dnl.xp.ndarray
        First and second derivatives of the experimental and
        theoretical intensities. Must be zero outside the overlap.
    no_overlap : dnl.xp.ndarray
        Boolean array, True for the points outside the overlap.

    Returns
    -------
    dnl.xp.ndarray
        The integrand, zero outside the overlap.
    """
    exp_1, exp_2 = exp_derivs
    theo_1, theo_2 = theo_derivs
    numerators = abs(beam_normalization * theo_2 - exp_2) * abs(
        beam_normalization * theo_1 - exp_1
    )
    numerators = dnl.xp.where(no_overlap, 0.0, numerators)

    # epsilon: largest experimental |I'| in the overlap
    epsilon = abs(exp_1) + dnl.xp.max(abs(exp_1), axis=0)
    epsilon = dnl.xp.clip(epsilon, 1e-12, None)
    epsilon = dnl.xp.where(no_overlap, 1.0, epsilon)
    return numerators / epsilon


def _zj_beam_rfactors(data_1, data_2, overlap, energy_step):
    """Return the ZJ R factor of each beam.

    Parameters
    ----------
    data_1, data_2 : tuple of dnl.xp.ndarray
        Intensity, first and second derivative of the experimental and
        theoretical data, with shape (n_energies, n_beams). May be NaN
        outside the overlap.
    overlap : dnl.xp.ndarray
        Boolean array, True where neither dataset is NaN.
    energy_step : float
        The energy step of the data grid.

    Returns
    -------
    dnl.xp.ndarray
        One value per beam. NaN for beams without overlap or with
        vanishing integrated intensity.
    """
    # set all values to 0 in regions where either dataset is nan
    # to avoid issues in the calculation
    no_overlap = dnl.xp.logical_not(overlap)
    exp_intensity, *exp_derivs = (
        dnl.xp.where(no_overlap, 0.0, data) for data in data_1
    )
    theo_intensity, *theo_derivs = (
        dnl.xp.where(no_overlap, 0.0, data) for data in data_2
    )

    exp_integral = utils.overlap_trapezoid(
        exp_intensity, energy_step, overlap, axis=0
    )
    # calculate normalization for each beam (nan without overlap)
    beam_normalization = exp_integral / utils.nonzero_or_nan(
        utils.overlap_trapezoid(theo_intensity, energy_step, overlap, axis=0)
    )
    quotient = _zj_quotient(
        beam_normalization, exp_derivs, theo_derivs, no_overlap
    )
    # Factor 0.027 for random correlation, Zannazi & Jona 1977
    prefactors = 1 / utils.nonzero_or_nan(exp_integral) / 0.027
    return prefactors * utils.overlap_trapezoid(
        quotient, energy_step, overlap, axis=0
    )


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

    The overall value, and the value of each group of beams, is the
    mean of the per-beam values weighted with the length of their
    overlap. Beams without overlap, and beams whose R is undefined
    (zero denominator, e.g., I_e = 0 on the overlap), are excluded
    from the overall and group values.

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
    data_1 = utils.coerce_to_energy_grid(
        data_and_derivatives_1, data_spline_1, energy_grid, derivs=2
    )
    # 2nd data also uses shifted grid if spline
    data_2 = utils.coerce_to_energy_grid(
        data_and_derivatives_2,
        data_spline_2,
        energy_grid - shift_2nd_spline,
        derivs=2,
    )
    overlap = utils.overlap_mask(data_1[0], data_2[0])
    beam_rfactors = _zj_beam_rfactors(data_1, data_2, overlap, energy_step)

    # Groups and the overall value are the means over beams weighted
    # by the length of the overlap.
    return group_beam_rfactors(
        beam_rfactors, overlap, energy_step, **kwargs
    )
