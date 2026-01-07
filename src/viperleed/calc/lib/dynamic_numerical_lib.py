"""viperleed.calc.lib.dynamic_numerical_lib.

This module provides dynamic numerical library support, allowing for
switching between NumPy and JAX for some array computations.

To swap the numerical library in use, set the corresponding variable
from this module.

For example, to use JAX instead of NumPy, you can do:

    from viperleed.calc.lib import dynamic_numerical_lib as dnl
    import jax.numpy as jnp
    dnl.xp = jnp
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
