"""Tests for <type> token of DISPLACEMENTS file."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-05-13'

import pytest

from viperleed.calc.classes.perturbation_type import PerturbationType
from viperleed.calc.files.new_displacements.tokens.type import TypeToken


@pytest.mark.parametrize(
    'type_str, expected_enum',
    [
        ('geo', PerturbationType.GEO),
        ('vib', PerturbationType.VIB),
        ('occ', PerturbationType.OCC),
        ('dom', PerturbationType.DOM),
    ],
)
def test_valid_type_strings(type_str, expected_enum):
    token = TypeToken(type_str)
    assert token.type is expected_enum
    # str representation should mention the enum
    rep = str(token)
    assert 'TypeToken' in rep
    assert expected_enum.name in rep or expected_enum.value in rep


def test_equality_and_type_mismatch():
    a = TypeToken('geo')
    b = TypeToken('geo')
    c = TypeToken('vib')
    assert a == b
    assert a != c
    # comparing to non-TypeToken always returns False
    assert a != 'geo'
    assert a is not None
