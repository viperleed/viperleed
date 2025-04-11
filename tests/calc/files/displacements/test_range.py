"""Tests for module files/displacements/range."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-11'


import pytest

from viperleed_jax.files.displacements.range import DisplacementsRange


def test_equal_ranges_with_step():
    r1 = DisplacementsRange(0.0, 1.0, 0.1)
    r2 = DisplacementsRange(0.0, 1.0, 0.1)
    assert r1 == r2


def test_equal_ranges_with_step_within_eps():
    r1 = DisplacementsRange(0.0, 1.0, 0.1)
    r2 = DisplacementsRange(0.0 + 1e-7, 1.0 - 1e-7, 0.1 + 1e-7)
    assert r1 == r2


def test_ranges_with_step_and_without_step_not_equal():
    r1 = DisplacementsRange(0.0, 1.0, 0.1)
    r2 = DisplacementsRange(0.0, 1.0)
    assert r1 != r2


def test_equal_ranges_without_step():
    r1 = DisplacementsRange(0.0, 1.0)
    r2 = DisplacementsRange(0.0, 1.0)
    assert r1 == r2


def test_ranges_without_step_not_equal_different_stop():
    r1 = DisplacementsRange(0.0, 1.0)
    r2 = DisplacementsRange(0.0, 2.0)
    assert r1 != r2


def test_repr_with_step():
    r = DisplacementsRange(0.1, 0.9, 0.2)
    assert repr(r) == '(start=0.1, stop=0.9, step=0.2)'


def test_repr_without_step():
    r = DisplacementsRange(-1.0, 1.0)
    assert repr(r) == '(start=-1.0, stop=1.0)'
