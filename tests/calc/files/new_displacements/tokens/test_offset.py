"""Tests for module offset of viperleed.calc.files.new_displacements.tokens."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-12'
__license__ = 'GPLv3+'


import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.files.new_displacements.tokens.offset import (
    OffsetToken,
    OffsetTokenParserError,
)


@pytest.mark.parametrize(
    'input_str, exp_offset',
    [
        ('1.0', 1.0),
        ('  -2.5  ', -2.5),
        ('3e-1', 0.3),
        ('4E+2', 400.0),
        # 2D
        ('1.0 2.0', np.array([1.0, 2.0])),
        # 3D
        ('1.0 2.0 3.0', np.array([1.0, 2.0, 3.0])),
    ],
)
def test_init_valid(input_str, exp_offset):
    token = OffsetToken(input_str)
    assert isinstance(token, OffsetToken)
    assert token.offset == pytest.approx(exp_offset)
    # repr
    assert str(token) == f'OffsetToken(offset={token.offset})'


@parametrize('bad', ['', '1 2 3 4', 'a b c'])
def test_init_invalid_count(bad):
    with pytest.raises(OffsetTokenParserError):
        OffsetToken(bad)


@parametrize('bad', ['foo', '1.0x', '--3'])
def test_init_non_numeric(bad):
    with pytest.raises(OffsetTokenParserError) as exc:
        OffsetToken(bad)
    assert 'Non-numeric value' in str(exc.value)


def test_from_float():
    a = OffsetToken.from_floats(5.0)
    assert a.offset == pytest.approx(5.0)


def test_from_floats_and_equality():
    a = OffsetToken('5.0')
    b = OffsetToken('5.0000001')
    assert a == b
    c = OffsetToken.from_floats(5.1)
    assert a != c
