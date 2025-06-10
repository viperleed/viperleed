"""Tests for <element> token."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-05-13'

# test_element_token.py
import pytest

from viperleed_jax.files.displacements.tokens.element import (
    ElementToken, ElementTokenParserError)


@pytest.mark.parametrize(
    'input_str, exp_z, exp_symbol',
    [
        ('H', 1, 'H'),
        ('h', 1, 'H'),
        (' He ', 2, 'He'),
    ],
)
def test_valid_elements(input_str, exp_z, exp_symbol):
    tok = ElementToken(input_str)
    assert tok.atomic_number == exp_z
    assert tok.symbol == exp_symbol


def test_invalid_element_raises():
    with pytest.raises(ElementTokenParserError) as exc:
        ElementToken('Xx')
    msg = str(exc.value)
    assert 'Could not parse chemical element' in msg
    assert 'Xx' in msg


def test_equality_and_type_mismatch():
    a = ElementToken('H')
    b = ElementToken('h')
    c = ElementToken('He')
    assert a == b  # same atomic number
    assert a != c
    assert a != 1  # different type
