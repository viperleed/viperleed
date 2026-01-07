"""Module spline_interpolation of viperleed.calc.lib.

Helper functions for the interpolation of ragged arrays.
"""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-02-19'
__license__ = 'GPLv3+'


from viperleed.calc.lib import dynamic_numerical_lib as dnl

def make_1d_ragged_cubic_spline(
    x, y, axis=0, bc_type='not-a-knot', extrapolate=False
):
    """Construct a piecewise cubic spline interpolator with ragged edges.

    The interpolator uses a cubic spline to interpolate data.
    """
    if x.ndim > 1 or y.ndim > 1:
        raise ValueError('x and y must be 1-dimensional arrays.')
    y_mask = dnl.xp.isnan(y)
    x_subarray, y_subarray = x[~y_mask], y[~y_mask]
    start_index = dnl.xp.where(~y_mask)[0][0]
    subarray_spline = dnl.CubicSpline(
        x_subarray, y_subarray, axis, bc_type, extrapolate
    )

    return subarray_spline, start_index


def interpolate_ragged_array(
    x, y, axis=0, bc_type='not-a-knot', extrapolate=False
):
    all_coeffs = dnl.xp.full(
        (4, y.shape[0], y.shape[1]), fill_value=dnl.xp.nan
    )
    for dim in range(y.shape[1]):
        _y = y[:, dim]
        all_nans = dnl.xp.all(dnl.xp.isnan(_y))
        if all_nans:
            continue
        spline, start_id = make_1d_ragged_cubic_spline(
            x, y[:, dim], axis=0, bc_type=bc_type, extrapolate=None
        )
        all_coeffs[:, start_id : start_id + spline.c.shape[1], dim] = spline.c
    return dnl.PPoly.construct_fast(all_coeffs, x)
