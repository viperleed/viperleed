import pytest

from viperleed_jax.files.displacements.tokens.offset import (
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
    ],
)
def test_init_valid(input_str, exp_offset):
    tok = OffsetToken(input_str)
    assert isinstance(tok, OffsetToken)
    assert tok.offset == pytest.approx(exp_offset)
    # repr
    assert repr(tok) == f'OffsetToken(offset={tok.offset})'


def test_init_invalid_count():
    for bad in ['', '1 2', '1 2 3']:
        with pytest.raises(OffsetTokenParserError) as exc:
            OffsetToken(bad)
        assert 'Invalid offset format' in str(exc.value)


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
