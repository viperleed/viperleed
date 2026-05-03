"""viperleed.calc.lib.rfactor.utils."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

import numpy as np

from viperleed.calc.lib import dynamic_numerical_lib as dnl


def nansum_trapezoid(y, dx, axis=-1):
    """Calculate the trapezoidal integral, ignoring NaNs."""
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
        two sets of beams to each other (usually experimental and
        theoretical). Should be a 1D array of integers with shape
        (n_beams,). Any beams assigned a value of -1 are removed from
        the beam array.

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
    # force numpy to ensure it can be used as an index array
    beam_corr = np.array(beam_correspondence, dtype=np.int32)

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


def coerce_to_energy_grid(
    data_and_derivatives, data_spline, energy_grid, derivs
):
    """Coerce data or splines to the energy grid.

    Utitlity function to allow flexible input of either pre-computed
    data and derivatives or splines that can be evaluated on the energy grid.
    If both are passed, data_and_derivatives will be used and data_spline
    will be ignored.

    Parameters
    ----------
    data_and_derivatives : tuple of arrays or None
        Pre-computed data and derivatives, expected to be a tuple of
        arrays (intensity, 1st derivative, ..., nth derivative) on the
        energy_grid. It is not checked (and not possible to check)
        whether the data is actually on the energy_grid.
        If None, data_spline will be used to evaluate the data.
    data_spline : scipy.interpolate.CubicSpline or similar or None
        Spline of the data, which will be evaluated on energy_grid.
        Ignored if data_and_derivatives is not None.
    energy_grid : array
        The grid on which data should be evaluated. Only used if
        data_and_derivatives is None.
    derivs : int
        The number of derivatives to compute. Only used if
        data_and_derivatives is None.

    Returns
    -------
    tuple of arrays
        Tuple of arrays (intensity, 1st derivative, ..., nth derivative)
        on the energy_grid, either taken from data_and_derivatives or
        evaluated from data_spline.

    Raises
    ------
    TypeError
        If both data_and_derivatives and data_spline are None.
    """
    if data_and_derivatives is None and data_spline is None:
        raise TypeError(
            'R_pendry requires either data splines or '
            'pre-computed data_and_derivatives arrays.'
        )
    if data_and_derivatives is None:
        # evaluate splines on energy grid to get data and derivatives
        sampled_data = data_spline(energy_grid)
        sampled_derivs = [
            data_spline.derivative(nu=deriv)(energy_grid)
            for deriv in range(derivs)
        ]
        return (sampled_data, *sampled_derivs)

    # otherwise just return data and derivatives truncated to the
    # requested number of derivatives
    return data_and_derivatives[: derivs + 1]
