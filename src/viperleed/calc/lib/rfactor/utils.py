"""viperleed.calc.lib.rfactor.utils."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

import numpy as np

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
    shifts = dnl.stop_gradient(shifts)

    # broadcast shifts and add to theo_intensity
    return theo_intensity + shifts


def average_beam_array(beam_array, beam_correspondence):
    """Average the beam array over the beam correspondence.

    Parameters
    ----------
    beam_array : array_like
        The beam array to average, shape (n_energies, n_beams).
    beam_correspondence : tuple
        A tuple containing the beam correspondence, which maps the
        experimental beams to the theoretical beams. It should be a 1D array
        of integers with shape (n_beams,).

    Returns
    -------
    array_like
        The averaged beam array, shape (n_energies, n_averaged_beams).

    Raises
    ------
    ValueError
        If the number of beams in the beam array does not match the length of
        the beam correspondence.
    """
    # convert beam correspondence to numpy array
    beam_corr = np.array(beam_correspondence, dtype=np.int32)  # force numpy

    # check if beam_array and beam_corr have the same number of beams
    if beam_array.shape[1] != beam_corr.shape[0]:
        raise ValueError(
            'Beam array and beam correspondence must have the same number of beams.'
        )

    # beam_corr may contain -1 for beams that have no correspondence; these
    # need to be filtered out
    valid_beam_mask = beam_corr != -1
    beam_corr = beam_corr[valid_beam_mask]
    filtered_beam_array = beam_array[:, valid_beam_mask]
    # determine the number of averaged beams after filtering
    n_averaged_beams = np.unique(beam_corr).size  # force numpy

    # get weights for averaging
    ones = dnl.xp.ones_like(beam_corr, dtype=float)
    summed = dnl.segment_sum(ones, beam_corr, num_segments=n_averaged_beams)
    averaged_beam_weights = dnl.xp.reciprocal(summed)

    # sum beams according to the beam correspondence
    mix_beams_vmap = dnl.vmap(dnl.segment_sum, in_axes=(0, None, None))
    averaged = mix_beams_vmap(filtered_beam_array, beam_corr, n_averaged_beams)

    # apply the averaged weights
    return averaged * averaged_beam_weights[dnl.xp.newaxis, :]
