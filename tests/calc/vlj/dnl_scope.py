"""Scoping of the numerical backend swapped by the viperleed_jax plugin.

Importing viperleed_jax.tensor_calculator replaces attributes of
viperleed.calc.lib.dynamic_numerical_lib with JAX equivalents, for
as long as the process lives. This module confines the replacement to
the tests that need it.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-26'
__license__ = 'GPLv3+'

from contextlib import contextmanager
import sys

from viperleed.calc.lib import dynamic_numerical_lib as dnl

PATCHED_NAMES = (
    'xp', 'stop_gradient', 'bincount', 'vmap', 'segment_sum',
    'CubicSpline', 'PPoly',
    )
_PLUGIN_MODULE = 'viperleed_jax.tensor_calculator'
_plugin_values = {}  # What the plugin sets, once it has been imported


def _snapshot():
    """Return the current values of the attributes the plugin patches."""
    return {name: getattr(dnl, name) for name in PATCHED_NAMES}


def _apply(values):
    """Set the attributes of dnl to values."""
    for name, value in values.items():
        setattr(dnl, name, value)


@contextmanager
def scoped_dnl_patch():
    """Undo, on exit, any swap of the numerical backend of dnl.

    If the plugin was imported before, its swap is re-applied on entry,
    as importing it again does not repeat it.
    """
    originals = _snapshot()
    _apply(_plugin_values)
    try:
        yield
    finally:
        if not _plugin_values and _PLUGIN_MODULE in sys.modules:
            _plugin_values.update(_snapshot())
        _apply(originals)
