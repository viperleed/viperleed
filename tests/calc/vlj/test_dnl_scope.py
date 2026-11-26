"""Tests for the scoping of the plugin's numerical backend."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-26'
__license__ = 'GPLv3+'

from importlib import import_module
from importlib.util import find_spec

import numpy as np
import pytest

from viperleed.calc.lib import dynamic_numerical_lib as dnl

from .dnl_scope import scoped_dnl_patch

pytestmark = pytest.mark.skipif(find_spec('viperleed_jax') is None,
                                reason='viperleed_jax plugin not installed')


def test_patch_confined_to_scope():
    """Check that the plugin's swap is undone, and re-applied on entry."""
    for _ in range(2):  # Second time the plugin is already imported
        with scoped_dnl_patch():
            import_module('viperleed_jax.tensor_calculator')
            assert dnl.xp is not np
        assert dnl.xp is np
