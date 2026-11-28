"""Module r_1 of viperleed.calc.lib.rfactor."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from . import utils
from .groups import group_rfactors


def r_1(
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
    Calculate the R1 R factor.

    The R1 R factor is kept for legacy reasons only. It cannot be
    selected via the R_FACTOR_TYPE parameter.

    Uses two sets of beam data, either passed as splines or as
    already calculated arrays of intensities.

    For each dataset, pass either data_spline or data_and_derivatives, but not
    both. Using splines for one and pre-computed data for the other is allowed.
    Note that data_and_derivatives are expected to be 3-tuples for uniform call
    signatures with R-factors that also use the derivatives; however, the
    derivatives are not used and can be left empty.

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
        Splines of the data, which will be evaluated on energy_grid. Ignored if
        the respective data_and_derivatives are passed.
    data_and_derivatives_1, data_and_derivatives_2 : 3-tuples of float arrays
        Pre-evaluated tuples (intensity, 1st derivative, 2nd derivative) at the
        energy_grid points. This avoids re-calculation of splines and is more
        JAX-friendly. The derivatives are not actually used, so passing
        (intensity, None, None) would be acceptable.
    shift_2nd_spline : float, optional
        Evaluate the 2nd spline on a shifted energy grid. This is meant for
        testing V0r variations. Note that this is NOT available when the 2nd
        dataset is passed as data_and_derivatives.

    Raises
    ------
    TypeError
        If neither data_spline nor data_and_derivatives is passed for
        either dataset, or if both are passed for the same dataset.
    """
    # Get data either as splines or as pre-computed arrays (mainly for JAX)
    data_1_intensity, *_ = utils.coerce_to_energy_grid(
        data_and_derivatives_1, data_spline_1, energy_grid, derivs=0
    )
    # 2nd data also shifted grid if spline
    data_2_intensity, *_ = utils.coerce_to_energy_grid(
        data_and_derivatives_2,
        data_spline_2,
        energy_grid - shift_2nd_spline,
        derivs=0,
    )

    # all integrals are restricted to the overlap of the two datasets
    overlap = utils.overlap_mask(data_1_intensity, data_2_intensity)

    # calculate normalization for each beam
    beam_normalization = utils.overlap_trapezoid(
        data_1_intensity, energy_step, overlap, axis=0
    ) / utils.overlap_trapezoid(
        data_2_intensity, energy_step, overlap, axis=0
    )

    numerators = utils.overlap_trapezoid(
        abs(data_1_intensity - beam_normalization * data_2_intensity),
        energy_step,
        overlap,
        axis=0,
    )
    denominators = utils.overlap_trapezoid(
        data_1_intensity**2, energy_step, overlap, axis=0
    )

    # calculate R-factor with requested grouping
    return group_rfactors(numerators, denominators, **kwargs)
