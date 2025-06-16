"""Tests for module perturbation_type."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'

import pytest

from viperleed.calc.classes.perturbation_type import (
    PerturbationMode,
    PerturbationModeError,
)


@pytest.mark.parametrize(
    'input_str, expected_enum',
    [
        ('geo', PerturbationMode.GEO),
        ('vib', PerturbationMode.VIB),
        ('occ', PerturbationMode.OCC),
    ],
)
def test_valid_enum_from_string(input_str, expected_enum):
    """Test valid string-to-enum conversion."""
    assert PerturbationMode.from_string(input_str) is expected_enum


@pytest.mark.parametrize(
    'enum_member, expected_str',
    [
        (PerturbationMode.GEO, 'geo'),
        (PerturbationMode.VIB, 'vib'),
        (PerturbationMode.OCC, 'occ'),
    ],
)
def test_enum_str_values(enum_member, expected_str):
    """Test that enum values behave like strings."""
    assert str(enum_member) == expected_str
    assert enum_member == expected_str  # thanks to str inheritance


@pytest.mark.parametrize(
    'invalid_str',
    [
        'geoo',
        'VIB',
        'oc',
        '',
        '123',
    ],
)
def test_invalid_enum_from_string(invalid_str):
    """Test that invalid strings raise ValueError."""
    with pytest.raises(PerturbationModeError, match='Unknown perturbation'):
        PerturbationMode.from_string(invalid_str)
