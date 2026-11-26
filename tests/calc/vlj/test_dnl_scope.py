"""Tests for the scoping of the plugin's numerical backend."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-26'
__license__ = 'GPLv3+'

from importlib import import_module
from importlib.util import find_spec

import pytest

from viperleed.calc.lib import dynamic_numerical_lib as dnl

from .dnl_scope import PATCHED_NAMES
from .dnl_scope import scoped_dnl_patch

pytestmark = pytest.mark.skipif(find_spec('viperleed_jax') is None,
                                reason='viperleed_jax plugin not installed')


@pytest.fixture(name='fixture_scoped_dnl_patch')
def fixture_no_outer_scope():
    """Replace the autouse scope: this test checks the scope itself.

    The autouse fixture re-applies the plugin's values on entry if any
    earlier test imported the plugin. The test would then start from
    the JAX backend, whatever the test order.
    """


def _snapshot():
    """Return the current values of the attributes the plugin patches."""
    return {name: getattr(dnl, name) for name in PATCHED_NAMES}


def test_patch_confined_to_scope():
    """Check that the plugin's swap is undone, and re-applied on entry."""
    before = _snapshot()
    jax_numpy = import_module('jax.numpy')
    assert before['xp'] is not jax_numpy
    for _ in range(2):  # Second time the plugin is already imported
        with scoped_dnl_patch():
            import_module('viperleed_jax.tensor_calculator')
            assert dnl.xp is jax_numpy
        assert _snapshot() == before
