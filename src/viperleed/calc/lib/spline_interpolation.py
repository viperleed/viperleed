"""Module spline_interpolation of viperleed.calc.lib.

Helper functions for the interpolation of ragged arrays.
"""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-02-19'
__license__ = 'GPLv3+'


import numpy as np

from viperleed.calc.lib import dynamic_numerical_lib as dnl


class CachedSpline:
    """
    Proxy that caches derivative splines and last evaluation.
    Assumes spline is pure: same x -> same y.
    """

    __slots__ = ('_spline', '_deriv_cache', '_last_x', '_last_y')

    def __init__(self, spline):
        self._spline = spline
        self._deriv_cache = {}  # key: (args, frozenset(kwargs.items()))
        self._last_x = None
        self._last_y = None

    def __call__(self, x):
        if x is self._last_x:
            return self._last_y
        y = self._spline(x)
        self._last_x = x
        self._last_y = y
        return y

    def derivative(self, *args, **kwargs):
        key = (args, tuple(sorted(kwargs.items())))
        if key not in self._deriv_cache:
            self._deriv_cache[key] = CachedSpline(
                self._spline.derivative(*args, **kwargs)
            )
        return self._deriv_cache[key]

    def __getattr__(self, name):
        # delegate anything else (e.g., extrapolate, c, x, etc.)
        return getattr(self._spline, name)


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
    # NB: this is intentionally using a Numpy array, even if the backend is
    # swapped out
    all_coeffs = np.full(
        (4, y.shape[0] - 1, y.shape[1]), fill_value=dnl.xp.nan
    )
    for dim in range(y.shape[1]):
        _y = y[:, dim]
        all_nans = dnl.xp.all(dnl.xp.isnan(_y))
        if all_nans:
            continue
        spline, start_id = make_1d_ragged_cubic_spline(
            x, _y, axis=0, bc_type=bc_type, extrapolate=None
        )
        all_coeffs[:, start_id : start_id + spline.c.shape[1], dim] = spline.c
    return dnl.PPoly.construct_fast(all_coeffs, x)
