"""viperleed.calc.lib.dynamic_numerical_lib.

This module provides dynamic numerical library support, allowing for
switching between NumPy and JAX for some array computations.

To swap the numerical library in use, set the corresponding variable
from this module.

For example, to use JAX instead of NumPy, you can do:

    from viperleed.calc.lib import dynamic_numerical_lib as dnl
    import jax
    dnl.xp = jax.numpy
"""

__authors__ = ('Alexander M. Imre (@amimre)',)
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
