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
    Assumes spline is pure: same x -> same y. The last evaluation is
    reused only if called again with the very same x object, which
    must not have been modified in place in the meantime.
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
        # this is only invoced if the attribute is not defined on the
        # CachedSpline itself
        return getattr(self._spline, name)


def make_1d_ragged_cubic_spline(
    x, y, axis=0, bc_type='not-a-knot', extrapolate=False
):
    """Construct a piecewise cubic spline with ragged edges.

    The interpolator uses a cubic spline to interpolate data.

    Parameters
    ----------
    x : array-like
        1D array of x-coordinates of the data points.
    y : array-like
        1D array of y-coordinates of the data points, which may contain
        NaNs to indicate missing data. NaNs are allowed only at the
        beginning and end of the array, and the non-NaN region must be
        contiguous.
    axis : int, optional
        Axis along which to interpolate. Default is 0.
    bc_type : str, optional
        Boundary condition type for the spline. Default is 'not-a-knot'.
    extrapolate : bool or str, optional
        Whether to extrapolate beyond the data range. Default is False.
        If True, extrapolation is done according to the specified
        boundary conditions. See scipy.interpolate for more details.

    Returns
    -------
    dnl.CubicSpline
        A cubic spline object from the selected numerical library.
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


def _shift_cubic(coeffs, shift):
    """Return the coefficients of a cubic re-expanded around x0 + shift.

    Parameters
    ----------
    coeffs : numpy.ndarray
        Coefficients (c3, c2, c1, c0) of a cubic in powers of (x - x0),
        in the order used by scipy.interpolate.PPoly.
    shift : float
        Displacement of the new expansion point from x0.

    Returns
    -------
    numpy.ndarray
        Coefficients of the same cubic in powers of (x - x0 - shift).
    """
    c_3, c_2, c_1, c_0 = coeffs
    return np.array([
        c_3,
        c_2 + 3 * c_3 * shift,
        c_1 + 2 * c_2 * shift + 3 * c_3 * shift**2,
        c_0 + c_1 * shift + c_2 * shift**2 + c_3 * shift**3,
    ])


def interpolate_ragged_array(
    x, y, axis=0, bc_type='not-a-knot', extrapolate=False
):
    """Interpolate a ragged array along a specified axis.

    This function constructs a piecewise cubic spline interpolator for each
    column of the input array `y`, allowing for NaN values that indicate
    missing data. The NaN values are allowed only at the beginning and end of
    the array, and the non-NaN region must be contiguous. The function returns
    a single piecewise polynomial object that can be evaluated on the original
    x values.

    Parameters
    ----------
    x : array-like
        1D array of x-coordinates of the data points.
    y : array-like
        2D array of y-coordinates of the data points, which may contain
        NaNs to indicate missing data. NaNs are allowed only at the
        beginning and end of each column, and the non-NaN region must be
        contiguous.
    axis : int, optional
        Axis along which to interpolate. Default is 0.
    bc_type : str, optional
        Boundary condition type for the spline. Default is 'not-a-knot'.
    extrapolate : bool or str, optional
        Whether to extrapolate beyond the data range. Default is False.
        If True, extrapolation is done according to the specified
        boundary conditions. See scipy.interpolate for more details.

    Returns
    -------
    dnl.PPoly
        A piecewise polynomial object from the selected numerical
        library that can be evaluated on a 1D array of x values
        to return the interpolated y values.
    """
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
            x, _y, axis=axis, bc_type=bc_type, extrapolate=extrapolate
        )
        n_intervals = spline.c.shape[1]
        end_id = start_id + n_intervals  # Index of the last knot of the beam
        all_coeffs[:, start_id:end_id, dim] = spline.c
        if end_id < len(x) - 1:
            # The interval starting at the last knot has no data of its
            # own. Continue the last cubic there, so that the beam's
            # last data point is not lost when evaluating at the knot.
            all_coeffs[:, end_id, dim] = _shift_cubic(
                np.asarray(spline.c[:, -1]), x[end_id] - x[end_id - 1]
            )
    return dnl.PPoly.construct_fast(all_coeffs, x)
