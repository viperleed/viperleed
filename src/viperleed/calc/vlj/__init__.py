"""ViperLEED JAX plugin environment.
"""

__authors__ = (
    "Alexander M. Imre (@amimre)",
)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-07-01"
__license__ = "GPLv3+"

try:
    import viperleed_jax
except (ModuleNotFoundError, ImportError):
    VLJ_AVAILABLE = False
else:
    VLJ_AVAILABLE = True
