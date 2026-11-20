"""Module dynamic_numerical_lib of viperleed.calc.lib.

This module provides dynamic numerical library support, allowing for
switching between NumPy and JAX for some array computations.

To swap the numerical library in use, set the corresponding variable
from this module.

For example, to use JAX instead of NumPy, you can do:

    from viperleed.calc.lib import dynamic_numerical_lib as dnl
    import jax
    dnl.xp = jax.numpy

This module further provides Numpy-compatible function patches for
JAX-specific functions such as vmap and stop_gradient.

Note that not all functionality is guaranteed to be fully compatible.
Scipy splines may be replaced with their JAX-compatible counterparts
from the interpax package.
"""

__authors__ = ('Alexandra M. Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

import numpy as _np  # NumPy is the default numerical library
from scipy.interpolate import (
    CubicSpline as _SciCubicSpline,
)
from scipy.interpolate import (
    PPoly as _SciPPoly,
)

xp = _np
CubicSpline = _SciCubicSpline
PPoly = _SciPPoly


def stop_gradient(array):
    """Stop gradient calculation of given array.

    This is a no-op when using NumPy as the numerical library: unlike
    jax.lax.stop_gradient, it returns `array` itself (not a copy), for
    any type of `array`.
    When using the R-factor calculation with a numerical library that
    supports automatic differentiation (such as JAX in viperleed-jax),
    this utility function should be patched with a method that stops
    gradient propagation. In JAX this can be done as
    ```
    from viperleed.calc.lib import dynamic_numerical_lib as dnl

    dnl.xp = jax.numpy
    dnl.stop_gradient = jax.lax.stop_gradient
    ```
    """
    return array


def bincount(x, weights, length):
    """Return the sum of `weights` for each value in `x`.

    Wrapper around numpy.bincount compatible with JAX argument names.
    Deviations from jax.numpy.bincount:
    - `weights` and `length` are required, and the third positional
      argument is `length`, not `minlength`. There is no `minlength`.
    - `length` is only a minimum length: values of `x` larger than
      or equal to `length` extend the output instead of being dropped.
    - Negative values of `x` raise ValueError instead of being
      counted as zero.

    Parameters
    ----------
    x : array_like of int
        Non-negative values to count.
    weights : array_like
        Weight of each value in `x`, same shape as `x`.
    length : int
        Minimum number of bins of the output.

    Returns
    -------
    numpy.ndarray
        Sum of `weights` for each value in `x`, of length
        max(`length`, max(`x`) + 1).
    """
    return _np.bincount(x, weights=weights, minlength=length)


def vmap(func, in_axes=0, out_axes=0):
    """Substitute for JAX vmap using NumPy.

    `func` is called once per element along the mapped axis, in a
    Python loop, and the outputs are stacked. Deviations from jax.vmap:
    - An int `in_axes` maps only the first positional argument; the
      others are passed unchanged. (JAX maps all arguments.)
    - Only a tuple, an int or None is accepted as `in_axes`: no lists
      or nested containers (pytrees) of axes.
    - Keyword arguments are never mapped: they are passed unchanged
      to each call. (JAX maps them along their first axis.)
    - `out_axes` must be an int. The outputs of `func` are stacked
      into a single array with numpy.stack, so `func` must return
      one array (or array-like) of the same shape on each call.
    - `axis_name`, `axis_size`, `spmd_axis_name` and `sum_match` are
      not implemented.

    Parameters
    ----------
    func : callable
        Function to map.
    in_axes : int or None or tuple, optional
        Axis or axes to map over. If an int, the axis of the first
        positional argument; the others are not mapped. If a tuple,
        one entry per positional argument: the axis to map over, or
        None if that argument is not mapped. At least one argument
        must be mapped, and all mapped axes must have the same
        length. Default is 0.
    out_axes : int, optional
        Axis of the output along which the results of each call are
        stacked. Default is 0.

    Returns
    -------
    callable
        Wrapped function that maps `func` over the specified axes.
        Calling it raises ValueError if `in_axes` does not match the
        positional arguments, if no argument is mapped, or if the
        mapped axes have different lengths.
    """

    def wrapped(*args, **kwargs):
        # Normalize in_axes
        if isinstance(in_axes, tuple):
            axes = in_axes
        else:
            axes = (in_axes,) + (None,) * (len(args) - 1)

        if len(axes) != len(args):
            raise ValueError(
                'in_axes must match number of positional arguments.'
            )

        # Find mapped axis length
        map_lens = []
        for a, ax in zip(args, axes):
            if ax is not None:
                a = _np.asarray(a)
                map_lens.append(a.shape[ax])
        if not map_lens:
            raise ValueError(
                'At least one argument must have in_axes != None.'
            )
        n = map_lens[0]
        if any(L != n for L in map_lens):
            raise ValueError('All mapped axes must have the same length.')

        outs = []
        for i in range(n):
            sliced_args = []
            for a, ax in zip(args, axes):
                if ax is None:
                    sliced_args.append(a)
                else:
                    a_arr = _np.asarray(a)
                    sliced_args.append(_np.take(a_arr, i, axis=ax))
            outs.append(func(*sliced_args, **kwargs))

        return _np.stack(outs, axis=out_axes)

    return wrapped


def segment_sum(values, segment_ids, num_segments):
    """NumPy equivalent of jax.ops.segment_sum.

    Deviations from jax.ops.segment_sum:
    - Only 1-D `values` are supported.
    - `num_segments` is required.
    - `indices_are_sorted`, `unique_indices`, `bucket_size` and
      `mode` are not implemented.
    - Segment ids outside [-`num_segments`, `num_segments`) raise
      IndexError, and negative ones count from the end, instead of
      being dropped.

    Parameters
    ----------
    values : array_like, shape (N,)
        Values to sum.
    segment_ids : array_like of int, shape (N,)
        Segment index for each value. Must satisfy
        0 <= segment_ids[i] < num_segments.
    num_segments : int
        Total number of segments.

    Returns
    -------
    np.ndarray
        Array of shape (num_segments,) where entry g is the sum of
        values with segment_ids == g.

    Raises
    ------
    ValueError
        If `values` is not 1-D, or if `segment_ids` has a different
        shape.
    """
    values = _np.asarray(values)
    segment_ids = _np.asarray(segment_ids, dtype=_np.int64)

    if values.ndim != 1:
        raise ValueError('values must be 1D')
    if segment_ids.shape != values.shape:
        raise ValueError('segment_ids must have the same shape as values')

    out = _np.zeros(num_segments, dtype=values.dtype)
    _np.add.at(out, segment_ids, values)
    return out
