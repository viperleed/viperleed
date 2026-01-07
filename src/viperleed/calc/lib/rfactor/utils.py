"""viperleed.calc.lib.rfactor.utils."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib.dynamic_numerical_lib import xp


def nansum_trapezoid(y, dx, axis=-1):
    y_arr = xp.moveaxis(y, axis, -1)
    # select the axis to integrate over
    return xp.nansum(y_arr[..., 1:] + y_arr[..., :-1], axis=-1) * dx * 0.5


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

    rolled_array = xp.roll(array, n_steps, axis=0)
    row_ids = xp.arange(n_energies).reshape(-1, 1)
    row_ids_tiled = xp.tile(row_ids, (1, n_beams))
    mask = xp.logical_or(
        row_ids_tiled < n_steps, row_ids >= n_energies + n_steps
    )
    return xp.where(mask, xp.nan, rolled_array)


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
    theo_intensity : xp.ndarray
        Interpolated theoretical intensity array of shape (n_energies, n_beams)
    exp_intensity : xp.ndarray
        Interpolated experimental intensity array of shape (n_energies, n_beams)

    Returns
    -------
    xp.ndarray
        Shifted theoretical intensity array of shape (n_energies, n_beams)
    """
    # Determine the mask for valid (non-NaN) intensities
    valid_exp_mask = ~xp.isnan(exp_intensity)
    valid_theo_mask = ~xp.isnan(theo_intensity)
    overlap_mask = xp.logical_and(valid_exp_mask, valid_theo_mask)

    # Calculate the minimum theoretical intensity in the valid regions
    masked_theo_intensity = xp.where(overlap_mask, theo_intensity, xp.nan)
    masked_theo_mins = xp.nanmin(masked_theo_intensity, axis=0)
    # Ensure that we do not consider NaNs in the minimum calculation
    min_theo_intensity = xp.where(
        xp.isnan(masked_theo_mins), 0.0, masked_theo_mins
    )
    # only shift if minimum is negative
    shifts = xp.where(min_theo_intensity < 0, -min_theo_intensity, 0.0)
    # stop gradient on shifts to avoid affecting optimization

    # broadcast shifts and add to theo_intensity
    return theo_intensity + shifts
