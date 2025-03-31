"""Tests for module symmetry_eps of viperleed.calc.classes.rparams.special."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-12-11'
__license__ = 'GPLv3+'

import operator

import pytest
from pytest_cases import fixture, parametrize

from viperleed.calc.classes.rparams.special.symmetry_eps import SymmetryEps


@fixture(name='value_and_z')
def factory_value_and_z():
    """Return the float and z values of a SymmetryEps."""
    def _make(eps):
        return float(eps), eps.z
    return _make


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
        'negative': ((-0.1,), ValueError),
        'negative z': ((0.1, -0.2), ValueError),
        }

    @parametrize('values,expected', valid.values(), ids=valid)
    def test_valid(self, values, expected, value_and_z):
        """Check correct interpretation of valid values."""
        eps = SymmetryEps(*values)
        assert value_and_z(eps) == expected

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

    def test_repr(self):
        """Check correct result of repr(eps)."""
        eps = SymmetryEps(0.1, 0.3)
        # pylint: disable=magic-value-comparison
        assert 'z=' in repr(eps)


class TestSymmetryEpsArithmetics:
    """Collection of tests for arithmetic operations on SymmetryEps objects."""

    valid = {  # self, other, operation, expected
        'add, no z': ((0.1,), SymmetryEps(0.1), operator.add,
                      SymmetryEps(0.2)),
        'add, z, float': ((0.1, 0.3), 0.1, operator.add,
                          SymmetryEps(0.2, 0.4)),
        'add, z': ((0.1, 0.3), SymmetryEps(0.15, 0.4), operator.add,
                   SymmetryEps(0.25, 0.7)),
        'sub, no z': ((0.2,), SymmetryEps(0.05), operator.sub,
                      SymmetryEps(0.15)),
        'sub, z, float': ((0.2, 0.1), SymmetryEps(0.05), operator.sub,
                          SymmetryEps(0.15, 0.05)),
        'sub, z': ((0.4, 0.3), SymmetryEps(0.1, 0.2), operator.sub,
                   SymmetryEps(0.3, 0.1)),
        'mul, no z': ((0.1,), 5.5, operator.mul, SymmetryEps(0.55)),
        'mul, z': ((0.1, 0.2), 3.2, operator.mul, SymmetryEps(0.32, 0.64)),
        'div, no z': ((0.5,), 5, operator.truediv, SymmetryEps(0.1)),
        'div, z': ((1.2, 0.3), 0.3, operator.truediv, SymmetryEps(4, 1)),
        'pow, no z': ((2,), 3, operator.pow, SymmetryEps(8)),
        'pow, z': ((2, 3), 3, operator.pow, SymmetryEps(8, 27)),
        'modulo, no z': ((4.2,), 2, operator.mod, SymmetryEps(0.2)),
        'modulo, z': ((4.2, 0.1), 2, operator.mod, SymmetryEps(0.2, 0.1)),
        }
    invalid = {  # self, other, operation, exc
        'type': ((0.2,), 'abc', operator.mul, TypeError),
        'mul': ((0.1,), SymmetryEps(0.1), operator.mul, TypeError),
        'div': ((0.1, 0.2), SymmetryEps(0.3), operator.floordiv, TypeError),
        'pow': ((0.3, 0.1), SymmetryEps(0.5), pow, TypeError),
        'neg result': ((0.1, 0.2), 0.5, operator.sub, ValueError),
        }

    @parametrize('args,other,operation,expected', valid.values(), ids=valid)
    # pylint: disable-next=too-many-arguments
    def test_valid(self, args, other, operation, expected, value_and_z):
        """Check that operation(eps, other) == expected."""
        eps = SymmetryEps(*args)
        result = operation(eps, other)
        assert value_and_z(result) == pytest.approx(value_and_z(expected))

    @parametrize('args,other,operation,exc', invalid.values(), ids=invalid)
    def test_invalid(self, args, other, operation, exc):
        """Check complaints when using invalid arithmetic operations."""
        eps = SymmetryEps(*args)
        with pytest.raises(exc):
            operation(eps, other)

    def test_reverse(self, value_and_z):
        """Check result of an arithmetic operation with swapped args."""
        eps = SymmetryEps(0.3, z=0.5)
        result = 3.2 + eps
        expected = SymmetryEps(3.5, 3.7)
        assert value_and_z(result) == pytest.approx(value_and_z(expected))

    valid_pos = {
        'with z': (2.7, 0.43),
        'without z': (3.5,),
        }

    @parametrize(args=valid_pos.values(), ids=valid_pos)
    def test_pos(self, args, value_and_z, subtests):
        """Check correct result of +eps."""
        eps = SymmetryEps(*args)
        result = +eps
        with subtests.test('result'):
            assert value_and_z(result) == pytest.approx(value_and_z(eps))
        with subtests.test('immutable'):
            assert result is not eps

    def test_neg(self):
        """Ensure -eps complains."""
        eps = SymmetryEps(2.7, 0.43)
        with pytest.raises(TypeError):
            _ = -eps

    three = {
        'round, no digits': ((2.7, 1.2), (None,), round, SymmetryEps(3, 1)),
        'round, digits': ((2.73, 1.26), (1,), round, SymmetryEps(2.7, 1.3)),
        'round no z': ((2.2,), (None,), round, SymmetryEps(2)),
        }

    @parametrize('args,others,operation,expected', three.values(), ids=three)
    # pylint: disable-next=too-many-arguments
    def test_valid_ternary(self, args, others, operation,
                           expected, value_and_z):
        """Check result of an operation with more three arguments."""
        eps = SymmetryEps(*args)
        result = operation(eps, *others)
        assert value_and_z(result) == pytest.approx(value_and_z(expected))

    def test_pow_modulo_raises(self):
        """Ensure that the 3-args form of pow() is unsupported."""
        eps = SymmetryEps(0.1)
        with pytest.raises(TypeError):
            pow(eps, 2, 3)
