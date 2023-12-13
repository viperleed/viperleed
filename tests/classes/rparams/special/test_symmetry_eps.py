"""Tests for module symmetry_eps of viperleed.tleedmlib.classes.rparams.special.

Created on 2023-12-11

@author: Michele Riva (@michele-riva)
"""

import operator

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
        assert (float(eps), eps.z) == expected

    @parametrize('values,exc', invalid.values(), ids=invalid)
    def test_invalid(self, values, exc):
        """Ensure exceptions are raised for invalid inputs."""
        with pytest.raises(exc):
            SymmetryEps(*values)

    ordering = {  # self, other, ordering_op
        'identical, no z': ((0.1,), SymmetryEps(0.1), operator.eq),
        'identical, z': ((0.1, 0.2), SymmetryEps(0.1, 0.2), operator.eq),
        'different, no z': ((0.1,), SymmetryEps(0.2), operator.ne),
        'different, z': ((0.1, 0.2), SymmetryEps(0.2, 0.1), operator.ne),
        'smaller, no z': ((0.1,), SymmetryEps(0.2), operator.lt),
        'equal, float': ((0.1,), 0.1, operator.eq),
        'different, z, float': ((0.1, 0.2), 0.1, operator.ne),
        'different, no z, float': ((0.1,), 0.2, operator.ne),
        'smaller, no z, float': ((0.1,), 0.2, operator.lt),
        'larger, z, float': ((0.1, 0.1), 0.05, operator.gt),
        'smaller, z, float': ((0.1, 0.1), 0.5, operator.lt),
        'less equal, z, float': ((0.1, 0.05), 0.1, operator.le),
        'not comparable': ((0.1,), 'abc', operator.ne),
        }
    not_ordered = {  # self, other, ordering_op to be negated
        'not comparable lt': ((0.1,), 'abc', operator.lt),
        'not comparable eq': ((0.1,), 'abc', operator.eq),
        }

    @parametrize('args,other,compare', ordering.values(), ids=ordering)
    def test_ordering(self, args, other, compare):
        """Check the expected outcome of the comparison with other."""
        eps = SymmetryEps(*args)
        assert compare(eps, other)

    @parametrize('args,other,compare', not_ordered.values(), ids=not_ordered)
    def test_not_ordered(self, args, other, compare):
        """Check that comparison with other is False."""
        eps = SymmetryEps(*args)
        try:
            result = compare(eps, other)
        except TypeError:  # Not comparable
            return
        assert not result

    hash_ = {  # self, other
        'identical, no z': ((0.1,), SymmetryEps(0.1)),
        'identical, z': ((0.1, 0.2), SymmetryEps(0.1, 0.2)),
        'equal, float': ((0.1,), 0.1),
        }

    @parametrize('args,other', hash_.values(), ids=hash_)
    def test_hash_equality(self, args, other, subtests):
        """Check that hash(eps) == hash(other) when they're equal."""
        eps = SymmetryEps(*args)
        with subtests.test('identity'):
            assert eps is not other
        with subtests.test('hash equality'):
            assert hash(eps) == hash(other)
