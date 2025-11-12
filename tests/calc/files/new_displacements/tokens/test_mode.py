"""Tests for module viperleed.calc.files.new_displacements.tokens.mode."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-05-13'

import pytest

from viperleed.calc.classes.perturbation_mode import PerturbationMode
from viperleed.calc.files.new_displacements.tokens.mode import ModeToken


@pytest.mark.parametrize(
    'type_str, expected_enum',
    [
        ('geo', PerturbationMode.GEO),
        ('vib', PerturbationMode.VIB),
        ('occ', PerturbationMode.OCC),
        ('dom', PerturbationMode.DOM),
    ],
)
def test_valid_mode_strings(type_str, expected_enum):
    token = ModeToken(type_str)
    assert token.mode is expected_enum
    # str representation should mention the enum
    rep = str(token)
    assert 'ModeToken' in rep
    assert expected_enum.name in rep or expected_enum.value in rep


def test_equality_and_mode_mismatch():
    a = ModeToken('geo')
    b = ModeToken('geo')
    c = ModeToken('vib')
    assert a == b
    assert a != c
    # comparing to non-ModeToken always returns False
    assert a != 'geo'
    assert a is not None
