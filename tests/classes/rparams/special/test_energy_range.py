"""Tests for module energy_range of viperleed.tleedmlib.rparams.special.

Created on 2023-10-27

@author: Michele Riva (@michele-riva)
"""

import itertools

import pytest
from pytest_cases import fixture, parametrize

from viperleed.tleedmlib.classes.rparams._defaults import NO_VALUE
from viperleed.tleedmlib.classes.rparams.special.energy_range import (
    EnergyRange
    )
from viperleed.tleedmlib.classes.rparams._defaults import NO_VALUE


class TestEnergyRange:
    """Collection of tests for the EnergyRange base class."""

    _class = EnergyRange

    @fixture(name='make_range', scope='class')
    def factory_make_range(self):
        """Return an EnergyRange object."""
        def _make(name):
            try:
                return self._class(*self.valid[name][0])
            except KeyError:
                pytest.skip(f'{name!r} is not valid')
                raise
        return _make

    valid = {      # init value,     expected
        '1 2 0.1': ((1.0, 2.0, 0.1), (1.0, 2.0, 0.1)),
        '1 2 0.2': ((1.0, 2.0, 0.2), (1.0, 2.0, 0.2)),
        'start==stop': ((1.0, 1.0, 0.1), (1.0, 1.0, 0.1)),
        'small step': ((1.0, 1.7, 1e-10), (1.0, 1.7, 1e-10)),
        'swapped': ((3.0, 1.0, -0.1), (1.0, 3.0, 0.1)),
        'swapped no step': ((3.0, -1.0), (-1.0, 3.0, NO_VALUE)),
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

    @parametrize('name,args', valid.items(), ids=valid)
    def test_equal(self, name, args, make_range):
        """Check correct result of the equality method."""
        *_, expected = args
        assert pytest.approx(make_range(name)) == expected
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

    @parametrize('name,args', valid.items(), ids=valid)
    def test_iter(self, name, args, make_range):
        """Check correct iteration of an EnergyRange."""
        *_, expected = args
        assert tuple(make_range(name)) == pytest.approx(expected)

    @parametrize('name,args', valid.items(), ids=valid)
    def test_properties(self, name, args, make_range):
        """Check min and max."""
        *_, expected = args
        energy_range = make_range(name)
        min_max = energy_range.min, energy_range.max
        assert min_max == pytest.approx(expected[:2])

    @parametrize('name,args', valid.items(), ids=valid)
    def test_from_value(self, name, args):
        """Check correct creation from a single sequence."""
        if name not in self.valid:
            pytest.skip(reason=f'{name!r} is not valid')
        value, expected = args
        energy_range = self._class.from_value(value)
        assert tuple(energy_range) == pytest.approx(expected)

    @parametrize('value,exc', invalid.values(), ids=invalid)
    def test_invalid_input(self, value, exc):
        """Check complaints when created from an invalid input."""
        with pytest.raises(exc):
            self._class.from_value(value)

    def test_defined(self, make_range):
        """Check correct value of defined @property."""
        assert make_range('swapped').defined
        assert not make_range('no step').defined
