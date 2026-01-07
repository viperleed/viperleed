"""viperleed.calc.lib.rfactor.utils."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib import dynamic_numerical_lib as dnl


def nansum_trapezoid(y, dx, axis=-1):
    y_arr = dnl.xp.moveaxis(y, axis, -1)
    # select the axis to integrate over
    return dnl.xp.nansum(y_arr[..., 1:] + y_arr[..., :-1], axis=-1) * dx * 0.5


def shift_theo_intensity_non_negative(theo_intensity, exp_intensity):
    """Shift the theoretical intensity so that it is non-negative.

    The array containing the theoretical intensity is calculated from a spline
    interpolation, which can lead to small negative values due to undershoots
    near the minima. This function shifts the theoretical intensities
    into the non-negative (physical) range by adding a constant value.
    The constant is chosen on a per-beam basis, and only considers the region of
    overlap between the theoretical and experimental intensities.

    Parameters
    ----------
    theo_intensity : dnl.xp.ndarray
        Interpolated theoretical intensity array of shape (n_energies, n_beams)
    exp_intensity : dnl.xp.ndarray
        Interpolated experimental intensity array of shape (n_energies, n_beams)

    Returns
    -------
    dnl.xp.ndarray
        Shifted theoretical intensity array of shape (n_energies, n_beams)
    """
    # Determine the mask for valid (non-NaN) intensities
    valid_exp_mask = ~dnl.xp.isnan(exp_intensity)
    valid_theo_mask = ~dnl.xp.isnan(theo_intensity)
    overlap_mask = dnl.xp.logical_and(valid_exp_mask, valid_theo_mask)

    # Calculate the minimum theoretical intensity in the valid regions
    masked_theo_intensity = dnl.xp.where(
        overlap_mask, theo_intensity, dnl.xp.nan
    )
    masked_theo_mins = dnl.xp.nanmin(masked_theo_intensity, axis=0)
    # Ensure that we do not consider NaNs in the minimum calculation
    min_theo_intensity = dnl.xp.where(
        dnl.xp.isnan(masked_theo_mins), 0.0, masked_theo_mins
    )
    # only shift if minimum is negative
    shifts = dnl.xp.where(min_theo_intensity < 0, -min_theo_intensity, 0.0)
    # stop gradient on shifts to avoid affecting optimization

    # broadcast shifts and add to theo_intensity
    return theo_intensity + shifts
