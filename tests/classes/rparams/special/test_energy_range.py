"""Tests for module energy_range of viperleed.tleedmlib.rparams.special.

Created on 2023-10-27

@author: Michele Riva (@michele-riva)
"""

import itertools

import pytest
from pytest_cases import fixture, parametrize

from viperleed.tleedmlib.classes.rparams.special.energy_range import (
    EnergyRange
    )
from viperleed.tleedmlib.classes.rparams._defaults import NO_VALUE


class TestEnergyRange:
    """Collection of tests for the EnergyRange base class."""

    @fixture(name='make_range', scope='session')
    def factory_make_range(self):
        """Return an EnergyRange object."""
        def _make(name):
            return EnergyRange(*self.valid[name][0])
        return _make

    valid = {      # init value,     expected
        '1 2 0.1': ((1.0, 2.0, 0.1), (1.0, 2.0, 0.1)),
        '1 2 0.2': ((1.0, 2.0, 0.2), (1.0, 2.0, 0.2)),
        'start==stop': ((1.0, 1.0, 0.1), (1.0, 1.0, 0.1)),
        'small step': ((0.0, 1.0, 1e-10), (0.0, 1.0, 1e-10)),
        'swapped': ((3.0, 1.0, -0.1), (1.0, 3.0, 0.1)),
        'swapped no step': ((3.0, 1.0), (1.0, 3.0, NO_VALUE)),
        'no step': ((1.0, 2.0), (1.0, 2.0, NO_VALUE)),
        }
    invalid = {
        'wrong type': (('invalid_input',), TypeError),
        'start>stop': ((2.0, 1.0, 0.1), ValueError),
        'step zero': ((1.0, 3.0, 0.0), ValueError),
        'infinite': ((float('-inf'), 5.0, 1.2), ValueError),
        'nan': ((1.0, float('nan'), 0.4), ValueError),
        'inf step': ((1.0, 5.0, float('inf')), ValueError),
        'nan step': ((1.0, 5.0, float('nan')), ValueError),
        }

    @parametrize(name=valid)
    def test_equal(self, name, make_range):
        """Check correct result of the equality method."""
        assert make_range(name) == make_range(name)

    @parametrize(names=itertools.combinations(valid, 2))
    def test_different(self, names, make_range):
        """Check non-equality of different ranges."""
        assert make_range(names[0]) != make_range(names[1])

    non_comparable = ('string', 1.78)

    @parametrize(invalid=non_comparable)
    def test_different_not_comparable(self, invalid, make_range):
        """Check that an invalid type cannot be compared."""
        assert make_range('1 2 0.1') != invalid

    @parametrize('name,expected', valid.items(), ids=valid)
    def test_iter(self, name, expected, make_range):
        """Check correct iteration of an EnergyRange."""
        assert tuple(make_range(name)) == expected[1]

    @parametrize('name,expected', valid.items(), ids=valid)
    def test_properties(self, name, expected, make_range):
        """Check min and max."""
        energy_range = make_range(name)
        assert (energy_range.min, energy_range.max) == expected[1][:2]

    @parametrize('value,expected', valid.values(), ids=valid)
    def test_from_value(self, value, expected):
        """Check correct creation from a single sequence."""
        energy_range = EnergyRange.from_value(value)
        assert tuple(energy_range) == expected

    @parametrize('value,exc', invalid.values(), ids=invalid)
    def test_invalid_input(self, value, exc):
        """Check complaints when created from an invalid input."""
        with pytest.raises(exc):
            EnergyRange.from_value(value)
