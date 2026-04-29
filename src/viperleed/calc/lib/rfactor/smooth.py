"""viperleed.calc.lib.rfactor.smooth."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib import dynamic_numerical_lib as dnl

from .pendry import r_pendry_from_y
from .utils import shift_theo_intensity_non_negative

DEFAULT_ALPHA = 4.0
DEFAULT_BETA = 0.15

def r_s(
    v0_imag,
    energy_step,
    energy_grid,
    data_spline_1=None,
    data_and_derivatives_1=None,
    data_spline_2=None,
    data_and_derivatives_2=None,
    shift_2nd_spline=0.0,  # only available if passed as spline
    alpha=DEFAULT_ALPHA,
    beta=DEFAULT_BETA,
    **kwargs,
):
    """
    Calculate the Smooth R-factor.

    See https://doi.org/10.1088/1361-648X/ae4af8.
    Uses two sets of beam data, either passed as splines or as already
    calculated arrays of intensities and derivatives. When comparing
    experimental and theoretical data, the theoretical data should be the
    second dataset, as this is also shifted to correct for non-negative values.

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
    alpha : float, optional
        Tuning parameter, default 4.0.
        Determines the influence of the intensity offset of the minimum, see
        https://doi.org/10.1088/1361-648X/ae4af8.
    beta : float, optional
        Tuning parameter, default 0.15.
        Determines the behaviour at a minimum reaching zero intensity, see
        https://doi.org/10.1088/1361-648X/ae4af8.

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
                'r_s requires either data splines or pre-computed '
                'data_and_derivatives arrays.'
            )
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
            raise TypeError(
                'r_s requires either data splines or pre-computed '
                'data_and_derivatives arrays.'
            )
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
    # shift data_2 to be non-negative in overlapping regions. This should fix
    #  any potential issues caused by undershooting splines.
    data_2_intensity = shift_theo_intensity_non_negative(data_2_intensity,
                                                         data_1_intensity)

    y_1 = y_s(
        data_1_intensity, data_1_derivative_1, data_1_derivative_2,
        v0_imag, alpha=alpha, beta=beta
    )
    y_2 = y_s(
        data_2_intensity, data_2_derivative_1, data_2_derivative_2,
        v0_imag, alpha=alpha, beta=beta
    )

    return r_pendry_from_y(y_1, y_2, energy_step, **kwargs)


def y_s(
    intensity,
    first_derivative,
    second_derivative,
    v0_imag,
    alpha=DEFAULT_ALPHA,
    beta=DEFAULT_BETA,
):
    """Calculate the Y function for the Smooth R-factor.

    Parameters
    ----------
    intensity : array
        The intensity values of the beam, evaluated on the energy grid.
    first_derivative : array
        The first derivative of the intensity with respect to energy
        evaluated on the same energy grid.
    second_derivative : array
        The second derivative of the intensity with respect to energy
        evaluated on the same energy grid.
    v0_imag : float
        The imaginary part of the inner potential.
    alpha : float, optional
        Tuning parameter, default 4.0 (module global,
        see also https://doi.org/10.1088/1361-648X/ae4af8).
    beta : float, optional
        Tuning parameter, default 0.15 (module global,
        see also https://doi.org/10.1088/1361-648X/ae4af8).
    """
    positive_d2 = second_derivative > 0
    # Note the value of 1.0 is arbitrary and does not affect the result.
    # The second derivative is only used in the case of a positive
    # curvature. However, this masking is required to avoid NaNs
    # in the JAX implementation.
    safe_d2 = dnl.xp.where(positive_d2, second_derivative, 1.0)
    intermediate_1 = (
        intensity / safe_d2
        - 0.5 * first_derivative**2 / safe_d2**2
    )
    intermediate_1 = (alpha / v0_imag**2) * intermediate_1 + beta
    intermediate_2 = intermediate_1 / dnl.xp.sqrt(1 + intermediate_1**2)

    numerator = first_derivative
    condition = dnl.xp.logical_and(positive_d2, intermediate_1 > 0)

    denominator = intensity**2 + 4 * v0_imag**2 * first_derivative**2
    conditional_denominator = (
        second_derivative**2 * v0_imag**4 * intermediate_2**2
    )
    denominator += condition * conditional_denominator
    denominator = dnl.xp.sqrt(denominator)
    return numerator / denominator
