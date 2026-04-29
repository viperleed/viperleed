"""viperleed.calc.lib.rfactor.r_1."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from .groups import group_rfactors
from .utils import nansum_trapezoid


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
    Calculate the R2 R factor.

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
    if data_and_derivatives_1 is None:
        if data_spline_1 is None:
            raise TypeError(
                'r_2 requires either data splines or '
                'pre-computed data_and_derivatives arrays.'
            )
        data_1_intensity = data_spline_1(energy_grid)
    else:
        data_1_intensity, _, _ = data_and_derivatives_1
    # 2nd data also uses shifted grid if spline
    shifted_grid = energy_grid - shift_2nd_spline
    if data_and_derivatives_2 is None:
        if data_spline_2 is None:
            raise TypeError(
                'r_2 requires either data splines or '
                'pre-computed data_and_derivatives arrays.'
            )
        # evaluate on shifted grid, allowing continuous shifts
        data_2_intensity = data_spline_2(shifted_grid)
    else:
        data_2_intensity, _, _ = data_and_derivatives_2

    # calculate normalization for each beam
    beam_normalization = nansum_trapezoid(
        data_1_intensity, energy_step, axis=0
    ) / nansum_trapezoid(data_2_intensity, energy_step, axis=0)

    numerators = nansum_trapezoid(
        (data_1_intensity - beam_normalization * data_2_intensity),
        energy_step,
        axis=0,
    )
    denominators = nansum_trapezoid(data_1_intensity**2, energy_step, axis=0)

    # calculate R-factor with requested grouping
    return group_rfactors(numerators, denominators, **kwargs)
