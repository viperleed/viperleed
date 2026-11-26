"""Tests for module helpers of viperleed.tests."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-28'
__license__ = 'GPLv3+'

import random

from viperleed.calc import psgen

from . import helpers


class TestPrivateRandom:
    """Tests for the private_random context manager."""

    def test_global_state_untouched(self):
        """Check that module draws leave the global state alone."""
        state = random.getstate()
        with helpers.private_random(psgen, seed=1) as rng:
            assert psgen.random is rng
            psgen.random.choice(range(100))
        assert random.getstate() == state

    def test_module_restored(self):
        """Check that the module gets its `random` back."""
        with helpers.private_random(psgen, seed=1):
            pass
        assert psgen.random is random

    def test_reproducible(self):
        """Check that the same seed gives the same draws."""
        draws = []
        for _ in range(2):
            with helpers.private_random(psgen, seed=1):
                draws.append([psgen.random.random() for _ in range(5)])
        assert draws[0] == draws[1]
