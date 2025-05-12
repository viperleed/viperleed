"""Tests for module files/displacements/range."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-11'


import pytest

from viperleed_jax.files.displacements.range import RangeToken


@pytest.mark.parametrize(
    'input_str, exp_start, exp_stop',
    [
        ('-0.1 0.1', -0.1, 0.1),
        ('1.0 2.0', 1.0, 2.0),
        ('-3 4', -3.0, 4.0),
        ('4 -3', 4.0, -3.0),
        ('1e0 2E+0', 1.0, 2.0),
        ('  5\t6  ', 5.0, 6.0),
    ],
)
def test_init_from_string_without_step(input_str, exp_start, exp_stop):
    dr = RangeToken(input_str)
    assert dr.start == exp_start
    assert dr.stop == exp_stop
    assert dr.step is None
    assert dr.has_step is False


@pytest.mark.parametrize(
    'input_str, exp_start, exp_stop, exp_step',
    [
        ('-1.5 3.0 0.5', -1.5, 3.0, 0.5),
        ('2 4 1', 2.0, 4.0, 1.0),
        ('0.1 1e2 1e-1', 0.1, 100.0, 0.1),
        ('-2E1 3e1 0.5', -20.0, 30.0, 0.5),
        ('  7  8.0   0.25  ', 7.0, 8.0, 0.25),
    ],
)
def test_init_from_string_with_step(input_str, exp_start, exp_stop, exp_step):
    dr = RangeToken(input_str)
    assert dr.start == pytest.approx(exp_start)
    assert dr.stop == pytest.approx(exp_stop)
    assert dr.step == pytest.approx(exp_step)
    assert dr.has_step is True


def test_init_stripping_whitespace():
    dr = RangeToken('  0   10   2  ')
    assert (dr.start, dr.stop, dr.step) == (0.0, 10.0, 2.0)


def test_init_invalid_number_of_parts():
    with pytest.raises(ValueError) as excinfo:
        RangeToken('42')
    assert 'Expected format' in str(excinfo.value)

    with pytest.raises(ValueError):
        RangeToken('1 2 3 4')


def test_init_non_numeric():
    with pytest.raises(ValueError) as excinfo:
        RangeToken('a b c')
    assert 'Non-numeric value' in str(excinfo.value)


def test_negative_step():
    with pytest.raises(ValueError) as excinfo:
        RangeToken('-0.2 0.2 -0.05')
    assert 'Step must be positive' in str(excinfo.value)


def test_from_floats():
    dr = RangeToken.from_floats(0, 5)
    assert isinstance(dr, RangeToken)
    assert dr.start == 0.0
    assert dr.stop == 5.0
    assert dr.step is None
    assert dr.has_step is False

    dr2 = RangeToken.from_floats(1, 4, 0.5)
    assert dr2.step == 0.5
    assert dr2.has_step is True


def test_equality_and_epsilon():
    dr1 = RangeToken('0 1 0.1')
    dr2 = RangeToken.from_floats(0 + 1e-7, 1 - 1e-7, 0.1 + 1e-7)
    assert dr1 == dr2

    dr3 = RangeToken.from_floats(0, 1, 0.1001)
    assert not (dr1 == dr3)

    # Without step
    dr4 = RangeToken('2 3')
    dr5 = RangeToken.from_floats(2 + 1e-7, 3 - 1e-7)
    assert dr4 == dr5

    # Mismatch of has_step
    dr6 = RangeToken('2 3 1')
    assert dr6 != dr4

    # Comparing with non-instance
    assert not (dr1 == (0, 1, 0.1))


def test_repr():
    dr = RangeToken('1 2')
    assert repr(dr) == 'RangeToken(start=1.0, stop=2.0)'

    drs = RangeToken('1 2 0.25')
    assert repr(drs) == 'RangeToken(start=1.0, stop=2.0, step=0.25)'