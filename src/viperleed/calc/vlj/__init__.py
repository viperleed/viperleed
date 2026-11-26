"""ViperLEED JAX plugin environment.
"""

__authors__ = (
    "Alexander M. Imre (@amimre)",
)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-07-01"
__license__ = "GPLv3+"

from importlib.util import find_spec

# Do not import viperleed_jax here: importing it (or any of its
# submodules) replaces the numerical backend of the R-factor library
# in viperleed.calc.lib with JAX, which would affect every user of
# this package, whether or not the JAX backend is requested.
VLJ_AVAILABLE = find_spec('viperleed_jax') is not None
