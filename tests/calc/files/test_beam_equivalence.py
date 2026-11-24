"""Tests for module beam_equivalence of viperleed.calc.files."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize
from quicktions import Fraction

from viperleed.calc.files.beam_equivalence import parse_beam_label

_HALF = Fraction(1, 2)
_THIRD = Fraction(1, 3)


_VALID_LABELS = {  # label, (h, k)
    'spaces': ('(1 0)', (1, 0)),
    'comma': ('(1, 0)', (1, 0)),
    'comma no space': ('(1,0)', (1, 0)),
    'vbar': ('(1 | 0)', (1, 0)),
    'vbar no space': ('(1|0)', (1, 0)),
    'vbar padded': ('( 1| 0)', (1, 0)),
    'tab': ('(\t1\t-1 )', (1, -1)),
    'outer whitespace': ('  (2 3)  ', (2, 3)),
    'negative': ('(-1 -2)', (-1, -2)),
    'plus sign': ('(+1 0)', (1, 0)),
    'fraction': ('(1/2 -1)', (_HALF, -1)),
    'negative fraction': ('(-1/2, 1/3)', (-_HALF, _THIRD)),
    'decimal': ('(0.5 0)', (_HALF, 0)),
    'decimal no zero': ('(.5 1.)', (_HALF, 1)),
    'six decimals': ('(0.333333 -0.666667)', (_THIRD, -2*_THIRD)),
    'unsnapped decimal': ('(0.1234 0)', (Fraction('0.1234'), 0)),
    }

_INVALID_LABELS = {
    'no parentheses': '1 0',
    'no closing': '(1 0',
    'no opening': '1 0)',
    'square brackets': '[1 0]',
    'one index': '(1)',
    'three indices': '(1 0 0)',
    'empty': '()',
    'text': '(a 0)',
    'double separator': '(1,,0)',
    'comma and vbar': '(1, | 0)',
    'no separator': '(1-1)',
    'double sign': '(--1 0)',
    'decimal fraction': '(0.5/2 0)',
    'exponent': '(1e0 0)',
    'zero denominator': '(1/0 0)',
    'trailing text': '(1 0) x',
    'two labels': '(1 0)(0 1)',
    }


class TestParseBeamLabel:
    """Tests for parse_beam_label."""

    @parametrize('label,expect', _VALID_LABELS.values(),
                 ids=_VALID_LABELS)
    def test_valid(self, label, expect):
        """Check the indices of a valid label."""
        indices = parse_beam_label(label)
        assert indices == expect
        assert all(isinstance(index, Fraction) for index in indices)

    def test_notations_equal(self):
        """Check that one beam in different notations gives one key.

        Mutant 'no snapping' (parse_beam_label without
        _snap_to_fraction) fails this test.
        """
        labels = '(0.5 0)', '(1/2, 0)', '(0.500000 | 0)', '(0.5000001 0)'
        keys = {parse_beam_label(label) for label in labels}
        assert keys == {(_HALF, 0)}

    @parametrize(label=_INVALID_LABELS.values(), ids=_INVALID_LABELS)
    def test_invalid(self, label):
        """Check that an invalid label raises ValueError."""
        with pytest.raises(ValueError, match='Invalid beam label'):
            parse_beam_label(label)
