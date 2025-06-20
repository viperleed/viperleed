import pytest
import numpy as np

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


def test_init_invalid_count():
    for bad in ['', '1 2 3 4', 'a b c']:
        with pytest.raises(OffsetTokenParserError) as exc:
            OffsetToken(bad)


def test_init_non_numeric():
    for bad in ['foo', '1.0x', '--3']:
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
