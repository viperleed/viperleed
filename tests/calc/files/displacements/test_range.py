"""Tests for module files/displacements/range."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-11'


import pytest

from viperleed_jax.files.displacements.range import DisplacementsRange


def test_init_from_string_minimum():
    dr = DisplacementsRange('1 2')
    assert dr.start == 1.0
    assert dr.stop == 2.0
    assert dr.step is None
    assert dr.has_step is False


def test_init_from_string_with_step():
    dr = DisplacementsRange('-1.5 3.0 0.5')
    assert dr.start == -1.5
    assert dr.stop == 3.0
    assert dr.step == 0.5
    assert dr.has_step is True


def test_init_stripping_whitespace():
    dr = DisplacementsRange('  0   10   2  ')
    assert (dr.start, dr.stop, dr.step) == (0.0, 10.0, 2.0)


def test_init_invalid_number_of_parts():
    with pytest.raises(ValueError) as excinfo:
        DisplacementsRange('42')
    assert 'Expected format' in str(excinfo.value)

    with pytest.raises(ValueError):
        DisplacementsRange('1 2 3 4')


def test_init_non_numeric():
    with pytest.raises(ValueError) as excinfo:
        DisplacementsRange('a b c')
    assert 'Non-numeric value' in str(excinfo.value)


def test_from_floats():
    dr = DisplacementsRange.from_floats(0, 5)
    assert isinstance(dr, DisplacementsRange)
    assert dr.start == 0.0
    assert dr.stop == 5.0
    assert dr.step is None
    assert dr.has_step is False

    dr2 = DisplacementsRange.from_floats(1, 4, 0.5)
    assert dr2.step == 0.5
    assert dr2.has_step is True


def test_equality_and_epsilon():
    dr1 = DisplacementsRange('0 1 0.1')
    dr2 = DisplacementsRange.from_floats(0 + 1e-7, 1 - 1e-7, 0.1 + 1e-7)
    assert dr1 == dr2

    dr3 = DisplacementsRange.from_floats(0, 1, 0.1001)
    assert not (dr1 == dr3)

    # Without step
    dr4 = DisplacementsRange('2 3')
    dr5 = DisplacementsRange.from_floats(2 + 1e-7, 3 - 1e-7)
    assert dr4 == dr5

    # Mismatch of has_step
    dr6 = DisplacementsRange('2 3 1')
    assert dr6 != dr4

    # Comparing with non-instance
    assert not (dr1 == (0, 1, 0.1))


def test_repr():
    dr = DisplacementsRange('1 2')
    assert repr(dr) == 'DisplacementsRange(start=1.0, stop=2.0)'

    drs = DisplacementsRange('1 2 0.25')
    assert repr(drs) == 'DisplacementsRange(start=1.0, stop=2.0, step=0.25)'