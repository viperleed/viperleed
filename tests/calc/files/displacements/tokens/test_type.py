"""Tests for module perturbation_type."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'

import pytest

from viperleed_jax.files.displacements.tokens.type import TypeToken, TypeTokenParserError
from viperleed_jax.files.displacements.perturbation_type import PerturbationType


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
    tok = TypeToken(type_str)
    assert tok.type is expected_enum
    # repr should mention the enum
    rep = repr(tok)
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
    assert a != None
