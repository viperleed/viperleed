"""Tests for module symmetry_eps of viperleed.tleedmlib.classes.rparams.special.

Created on 2023-12-11

@author: Michele Riva (@michele-riva)
"""

import pytest
from pytest_cases import parametrize

from viperleed.tleedmlib.classes.rparams import SymmetryEps


class TestSymmetryEps:
    """Collection of tests for SymmetryEps objects."""

    valid = {
        'no z': ((0.1,), (0.1, 0.1)),
        'z None': ((0.1, None), (0.1, 0.1)),
        'with z': ((0.1, 0.2), (0.1, 0.2)),
        }

    invalid = {
        'float': (('invalid',), TypeError),
        'float z': ((0.1, 'invalid'), TypeError),
        }

    @parametrize('values,expected', valid.values(), ids=valid)
    def test_valid(self, values, expected):
        """Check correct interpretation of valid values."""
        eps = SymmetryEps(*values)
        assert (eps, eps.z) == expected

    @parametrize('values,exc', invalid.values(), ids=invalid)
    def test_invalid(self, values, exc):
        """Ensure exceptions are raised for invalid inputs."""
        with pytest.raises(exc):
            SymmetryEps(*values)
