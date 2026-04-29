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

    This is a no-op when using NumPy as the numerical library.
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
    """Bincount implementation compatible with dynamic numerical libraries.

    This is a wrapper around `dnl.xp.bincount` to provide a consistent
    interface across different numerical libraries. Some libraries may
    use `minlength` instead of `length` as a parameter name.
    """
    return _np.bincount(x, weights=weights, minlength=length)


def vmap(func, in_axes=0, out_axes=0):
    """Substitute for JAX vmap using NumPy.

    Parameters
    ----------
    func : callable
        Function to map.
    in_axes : int or None or tuple
        Axis or axes to map over for each argument. If an int, the same
        axis is used for all arguments. If None, the argument is not
        mapped over. If a tuple, it should have the same length as the
        number of arguments, and specify the axis for each argument.
        See the JAX documentation for more details.
    out_axes : int
        Axis along which to stack outputs.

    Returns
    -------
    callable
        Wrapped function that maps `func` over the specified axis.
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
    """NumPy equivalent of jax.lax.segment_sum.

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
