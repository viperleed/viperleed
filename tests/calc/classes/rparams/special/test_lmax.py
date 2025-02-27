"""Tests for module l_max of viperleed.calc.classes.rparams.special."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-12-16'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture, parametrize

from viperleed.calc.classes.rparams.defaults import NO_VALUE
from viperleed.calc.classes.rparams.special.l_max import LMax

from .....helpers import not_raises


class TestLMaxValid:
    """Collection of tests for valid inputs to LMax objects."""

    valid = {
        'no value': tuple(),
        'single': (3,),
        'range': (2, 5),
        'swapped': (17, 3),
        'int-like float': (2.000000002, 8),
        }
    _min_max_bool = {
        'no value': (NO_VALUE, NO_VALUE, False),
        'single': (3, 3, True),
        'range': (2, 5, True),
        'swapped': (3, 17, True),
        }
    _str_repr = {
        'no value': ('', 'LMax()'),
        'single': ('3', 'LMax(3)'),
        'range': ('2-5', 'LMax(2, 5)'),
        }

    @fixture(name='make_lmax')
    def factory_make_lmax(self):
        """Return an LMax object from a key of self._valid."""
        def _make(key):
            return LMax.from_value(self.valid[key])
        return _make

    @parametrize(args=valid.values(), ids=valid)
    def test_creation(self, args):
        """Check successful creation of an LMax."""
        with not_raises(ValueError), not_raises(RuntimeError):
            LMax(*args)

    @parametrize(key=_min_max_bool)
    def test_bounds(self, key, make_lmax):
        """Check correct values of bounds of valid LMax objects."""
        lmax = make_lmax(key)
        *min_max, _ = self._min_max_bool[key]
        assert [lmax.min, lmax.max] == min_max

    @parametrize(key=_min_max_bool)
    def test_bool(self, key, make_lmax):
        """Check correct values of bounds of valid LMax objects."""
        lmax = make_lmax(key)
        *_, bool_ = self._min_max_bool[key]
        assert bool(lmax) is bool_

    @parametrize(key=_str_repr)
    def test_str(self, key, make_lmax):
        """Check correctness of string representation."""
        str_, *_ = self._str_repr[key]
        assert str(make_lmax(key)) == str_

    @parametrize(key=_str_repr)
    def test_repr(self, key, make_lmax):
        """Check correctness of string representation."""
        *_, repr_ = self._str_repr[key]
        assert repr(make_lmax(key)) == repr_

    _new_min = {  # new min, expected
        'increased': (3, (3, 5)),
        'at max': (5, (5, 5)),
        'beyond max': (12, (5, 12)),
        }

    @parametrize('new_min,expected', _new_min.values(), ids=_new_min)
    def test_set_min(self, new_min, expected, make_lmax):
        """Check correct assignment of a new minimum value."""
        lmax = make_lmax('range')
        lmax.min = new_min
        assert (lmax.min, lmax.max) == expected

    _new_max = {  # new max, expected
        'increased': (7, (2, 7)),
        'at min': (2, (2, 2)),
        'beyond min': (1, (1, 2)),
        }

    @parametrize('new_max,expected', _new_max.values(), ids=_new_max)
    def test_set_max(self, new_max, expected, make_lmax):
        """Check correct assignment of a new maximum value."""
        lmax = make_lmax('range')
        lmax.max = new_max
        assert (lmax.min, lmax.max) == expected


class TestLMaxInvalid:
    """Collection of tests for various invalid inputs for LMax objects."""

    invalid = {
        'non numeric': (('invalid',), TypeError),
        'not an integer': ((3.5, 17), TypeError),
        'negative': ((-5, 3), ValueError),
        'out of range': ((7, 100), ValueError),
        }

    @parametrize('args,exc', invalid.values(), ids=invalid)
    def test_invalid_creation(self, args, exc):
        """Check complaints at creation time with invalid arguments."""
        with pytest.raises(exc):
            LMax(*args)

    invalid_min = {  # new_min, exception
        'non numeric': ((1, 2), TypeError),
        'non int': (6.7, TypeError),
        'negative': (-5, ValueError),
        'out of range': (100, ValueError),
        }

    @parametrize('new_min,exc', invalid_min.values(), ids=invalid_min)
    def test_set_min_invalid(self, new_min, exc):
        """Check complaints when setting an invalid new minimum."""
        lmax = LMax(2)
        with pytest.raises(exc):
            lmax.min = new_min
