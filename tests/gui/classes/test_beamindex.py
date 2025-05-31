"""Tests for module beamindex of viperleed.gui.classes."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-03-03'
__license__ = 'GPLv3+'

from collections import namedtuple
import operator

import pytest
from pytest_cases import fixture, parametrize
from quicktions import Fraction

from viperleed.gui.classes.beamindex import BeamIndex


_IndexArgs = namedtuple('_IndexArgs', ('args', 'kwargs'))


@fixture(name='with_clean_cache', autouse=True)
def factory_clean_cache():
    """Yield BeamIndex with an empty cache."""
    old_cache = BeamIndex._cache.copy()
    BeamIndex.clear_cache()
    try:
        yield
    finally:
        BeamIndex._cache = old_cache


class TestBeamIndex:
    """Collection of tests for the BeamIndex class."""

    def test_cache(self, subtests):
        """Test caching functionality."""
        with subtests.test('no -1'):
            beam1 = BeamIndex(1, 1)
            beam2 = BeamIndex(1, 1)
            assert beam1 is beam2
        with subtests.test('with -1'):
            beam1 = BeamIndex(1, -1)
            beam2 = BeamIndex(1, -1)
            assert beam1 is not beam2

    def test_clear_cache(self):
        """Test cache clearing."""
        beam1 = BeamIndex(1, 1)
        beam1.clear_cache()
        beam2 = BeamIndex(1, 1)
        assert beam1 is not beam2

    _init = {  # args, kwargs, expect
        'one string': (
            _IndexArgs(('1/4,   -1/10',), {}),
            (Fraction(1, 4), Fraction(-1, 10))
            ),
        'sequence': (
            _IndexArgs(([1, Fraction(-3)],), {}),
            (Fraction(1), Fraction(-3))
            ),
        'two args': (
            _IndexArgs((2, 10), {}),
            (Fraction(2), Fraction(10))
            ),
        'numerators': (
            _IndexArgs((2, 3), {'from_numerators': True}),
            (Fraction(2), Fraction(3))
            ),
        'denominator': (
            _IndexArgs((-2/3, 6), {'denominator': 3}),
            (Fraction(-2, 3), Fraction(6))
            ),
        'denominator from numerators': (
            _IndexArgs((-2, 6), {'denominator': 3, 'from_numerators': True}),
            (Fraction(-2, 3), Fraction(2))
            ),
        }

    @parametrize('args,expect', _init.values(), ids=_init)
    def test_creation(self, args, expect):
        """Check creation of a valid BeamIndex."""
        beam = BeamIndex(*args.args, **args.kwargs)
        assert beam == expect

    _cache_problems = { # pairs of problematic beams
        ('2, -2', '-2, 2'): ((Fraction(2), Fraction(-2)),
                             (Fraction(-2), Fraction(2))),
        ('2, 2', '-2, -2'): ((Fraction(2), Fraction(2)),
                             (Fraction(-2), Fraction(-2))),
        ('3, -2', '-3, 0'): ((Fraction(3), Fraction(-2)),
                             (Fraction(-3), Fraction(0))),
        ('-3, -2', '3, 0'): ((Fraction(-3), Fraction(-2)),
                             (Fraction(3), Fraction(0))),
        ('3/2, 3/2', '-3/2, -3/2'): (
            (Fraction(3, 2), Fraction(3, 2)),
            (Fraction(-3, 2), Fraction(-3, 2)),
            ),
        ('3/2, -3/2', '-3/2, 3/2'): (
            (Fraction(3, 2), Fraction(-3, 2)),
            (Fraction(-3, 2), Fraction(3, 2)),
            ),
        }

    @parametrize('indices,expect', _cache_problems.items())
    def test_cache_problems(self, indices, expect):
        """Check correct handling of problematic hash values."""
        beams = tuple(BeamIndex(i) for i in indices)
        assert beams == expect
        
        # Invalidate cache, and do the opposite
        BeamIndex.clear_cache()
        beams = tuple(BeamIndex(i) for i in reversed(indices))
        assert beams == tuple(reversed(expect))

    _fmt = {  # str_arg, spec, expect
        'str': ('1|1', '', '1, 1'),
        'tuple format': ('3/2,5', '>10', '(    3/2, 5)'),
        'vbar auto fields': ('3/2,-5', 's', ' 3/2|-5  '),
        'vbar specified widths': ('3/2,-5', '(5,3)s', '    3/2  |   -5    '),
        'float': ('3/2,-5', '6.3f', '(     1.500,    -5.000)'),
        }

    @parametrize('arg,spec,expect', _fmt.values(), ids=_fmt)
    def test_format(self, arg, spec, expect):
        """Check expected outcome of __format__ method."""
        beam = BeamIndex(arg)
        assert format(beam, spec) == expect

    _fmt_width = {
        (1, -3): (2, 0),
        (Fraction(2, -5), Fraction(-40)): (3, 1),
        }

    @parametrize('args,expect', _fmt_width.items())
    def test_format_widths(self, args, expect):
        """Check correct identification of the char widths of num & den."""
        beam = BeamIndex(*args)
        assert beam.get_format_widths() == expect

    _num = {
        'integer': ((1, 2), (1, 2)),
        'fractional': ((Fraction(1, 5), Fraction(2, 9)), (1, 2)),
        }

    @parametrize('args,expect', _num.values(), ids=_num)
    def test_numerators(self, args, expect):
        """Check correctness of .numerators property."""
        beam = BeamIndex(*args)
        assert beam.numerators == expect

    _mul = {
        'int': (2, (3, 4)),
        'fraction': (Fraction(5, 4), (Fraction(15, 8), Fraction(5, 2))),
        }

    @parametrize('value,expect', _mul.values(), ids=_mul)
    def test_scaling(self, value, expect, subtests):
        """Check correct outcome of beam*value and value*beam."""
        with subtests.test('mul'):
            beam = BeamIndex('3/2,2')
            scaled = beam * value
            assert scaled == expect
        with subtests.test('rmul'):
            beam = BeamIndex('3/2,2')
            scaled = value * beam
            assert scaled == expect

    @parametrize('value,expect', _mul.values(), ids=_mul)
    def test_scaling_in_place(self, value, expect):
        """Check correct outcome of in-place scaling."""
        beam = BeamIndex('3/2|2')
        beam *= value
        assert beam == expect

    def test_str_repr(self):
        """Check __str__ and __repr__ methods."""
        beam = BeamIndex('1/2,3')
        assert str(beam) == '1/2, 3'
        assert repr(beam) == 'BeamIndex(1/2, 3)'


class TestBeamIndexRaises:
    """Collection of tests for invalid arguments to BeamIndex."""

    _init = {
        'too few': (_IndexArgs((1,), {}), TypeError),
        'too many': (_IndexArgs((1, 2, 3), {}), ValueError),
        'wrong separator': (_IndexArgs(('1;2',), {}), ValueError),
        'from numerator, not int': (
            _IndexArgs((1, None), {'from_numerators': True}),
            TypeError
            ),
        'wrong denominator': (
            _IndexArgs((1.5, -3.2), {'denominator': 3}),
            ValueError
            ),
        }

    @parametrize('args,exc', _init.values(), ids=_init)
    def test_creation(self, args, exc):
        """Check complaints when initializing with invalid arguments."""
        with pytest.raises(exc):
            BeamIndex(*args.args, **args.kwargs)

    _operation = {
        'add': (2,),
        'mul': 1.5,
        'imul': '2',
        }

    @parametrize('op,value', _operation.items(), ids=_operation)
    def test_operation_raises(self, op, value):
        """Check that BeamIndex cannot be extended."""
        beam = BeamIndex(1, 1)
        op_ = getattr(operator, op)
        with pytest.raises(TypeError):
            op_(beam, value)


