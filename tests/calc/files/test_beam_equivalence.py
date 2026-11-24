"""Tests for module beam_equivalence of viperleed.calc.files."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize
from quicktions import Fraction

from viperleed.calc.classes.beam_correspondence import BeamCorrespondence
from viperleed.calc.files.beam_equivalence import parse_beam_label
from viperleed.calc.files.beam_equivalence import read_beam_equivalence

_HALF = Fraction(1, 2)
_THIRD = Fraction(1, 3)


_VALID_LABELS = {  # pylint: disable=R6101
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


@fixture(name='read_text')
def fixture_read_text(tmp_path):
    """Return a function that reads BEAM_EQUIVALENCE contents."""
    def _read(contents):
        path = tmp_path / 'BEAM_EQUIVALENCE'
        path.write_text(contents, encoding='utf-8')
        return read_beam_equivalence(path)
    return _read


_VALID_FILES = {  # pylint: disable=R6101
    'one beam': ('(1 0) = (1 0)', {(1, 0): {(1, 0): 1}}),
    'no spaces': ('(1 0)=(1 0)+(0 1)', {(1, 0): {(1, 0): 1, (0, 1): 1}}),
    'many spaces': ('  (1 0)   =  (1 0)  +   (0 1)  ',
                    {(1, 0): {(1, 0): 1, (0, 1): 1}}),
    'weights': ('(1 0) = 3*(1 0) + 1 * (0 1)',
                {(1, 0): {(1, 0): 3, (0, 1): 1}}),
    'weight no spaces': ('(1 0) = 0.75*(1 0)+0.25*(0 1)',
                         {(1, 0): {(1, 0): 3, (0, 1): 1}}),
    'fraction weights': ('(1 0) = 2/3*(1 0) + 1/3*(0 1)',
                         {(1, 0): {(1, 0): 2, (0, 1): 1}}),
    'implicit weight 1': ('(1 0) = 0.5*(1 0) + (0 1)',
                          {(1, 0): {(1, 0): 1, (0, 1): 2}}),
    'equal weights': ('(1 0) = 2*(1 0) + 2*(0 1)',
                      {(1, 0): {(1, 0): 1, (0, 1): 1}}),
    'plus in label': ('(+1 0) = (+1 0) + (0 +1)',
                      {(1, 0): {(1, 0): 1, (0, 1): 1}}),
    'comments and blanks': (
        ('# header\n\n(1 0) = (1 0) ! one\n'
         '   \n% another\n(0 1) = (0 1) # two\n'),
                            {(1, 0): {(1, 0): 1}, (0, 1): {(0, 1): 1}}),
    'fractional beams': ('(1/2 0) = (0.5, 0) + (0 | 1/2)',
                         {(_HALF, 0): {(_HALF, 0): 1, (0, _HALF): 1}}),
    'shared beam': ('(1 0) = (1 0) + (0 1)\n(0 1) = (0 1)',
                    {(1, 0): {(1, 0): 1, (0, 1): 1}, (0, 1): {(0, 1): 1}}),
    'zero weight': ('(1 0) = (1 0) + 0*(0 1)', {(1, 0): {(1, 0): 1}}),
    'zero weight decimal': ('(1 0) = 0.0*(0 1) + (1 0)',
                            {(1, 0): {(1, 0): 1}}),
    }

# The rotational-domain example of issue #206, without and with
# weights (twice as much of the domain giving (1 0) and (-1 0))
_ISSUE_206 = (
    '! Rotational domains, see #206\n'
    '(1 0) = (1 0) + (-1 0) + (0 1) + (0 -1)\n'
    '(1 1) = (1 1) + (-1 -1) + (1 -1) + (-1 1)\n'
    )
_ISSUE_206_WEIGHTED = (
    '(1 0) = 2*(1 0) + 2*(-1 0) + (0 1) + (0 -1)\n'
    '(1 1) = (1 1) + (-1 -1) + (1 -1) + (-1 1)\n'
    )

_INVALID_FILES = {  # pylint: disable=R6101
    'no equal sign': ('(1 0) (1 0)', 1, 'Invalid line'),
    'two equal signs': ('(1 0) = (1 0) = (0 1)', 1, 'Invalid line'),
    'bad left label': ('(1) = (1 0)', 1, 'Invalid beam label'),
    'bad right label': ('(1 0) = (1 0 0)', 1, 'Invalid beam label'),
    'empty left': (' = (1 0)', 1, 'Invalid beam label'),
    'empty right': ('(1 0) = ', 1, 'Invalid right-hand side'),
    'trailing plus': ('(1 0) = (1 0) +', 1, 'Invalid right-hand side'),
    'double plus': ('(1 0) = (1 0) ++ (0 1)', 1, 'Invalid right-hand side'),
    'leading plus': ('(1 0) = + (1 0)', 1, 'Invalid right-hand side'),
    'minus term': ('(1 0) = (1 0) - (0 1)', 1, 'Invalid term'),
    'negative weight': ('(1 0) = -0.5*(1 0)', 1, 'Invalid term'),
    'minus after plus': ('(1 0) = (1 0) + -1*(0 1)', 1, 'Invalid term'),
    'weight after beam': ('(1 0) = (1 0)*0.5', 1, 'Invalid term'),
    'weight no star': ('(1 0) = 0.5(1 0)', 1, 'Invalid term'),
    'text weight': ('(1 0) = a*(1 0)', 1, 'Invalid term'),
    'exponent weight': ('(1 0) = 1e-1*(1 0)', 1, 'Invalid term'),
    'zero denominator': ('(1 0) = 1/0*(1 0)', 1, 'division by zero'),
    'missing star label': ('(1 0) = 0.5*', 1, 'Invalid term'),
    'all zero': ('(1 0) = (1 0)\n(0 1) = 0*(0 1) + 0.0*(1 0)',
                 2, 'no calculated beam with a non-zero weight'),
    'duplicate right': ('\n(1 0) = (1 0) + (0 1) + (1 0)',
                        2, 'more than once on the right-hand side'),
    'duplicate right snapped': ('(1 0) = (1/3 0) + (0.333333 0)',
                                1, 'more than once on the right-hand side'),
    'duplicate left': ('(1 0) = (1 0)\n# c\n(1 0) = (0 1)',
                       3, 'already assigned on line 1'),
    'duplicate left snapped': ('(1/3 0) = (1 0)\n(0.333333 | 0) = (0 1)',
                               2, 'already assigned on line 1'),
    'error after valid': ('(1 0) = (1 0)\n(0 1) = (0 1)\n(1 1) (1 1)',
                          3, 'Invalid line'),
    }


class TestReadBeamEquivalence:
    """Tests for read_beam_equivalence."""

    @parametrize('contents,expect', _VALID_FILES.values(), ids=_VALID_FILES)
    def test_valid(self, contents, expect, read_text):
        """Check the correspondence read from a valid file."""
        assert read_text(contents) == BeamCorrespondence(expect)

    def test_order(self, read_text):
        """Check that beam assignments keep the order of the file."""
        correspondence = read_text('(0 1) = (0 1)\n(1 0) = (1 0)')
        assert list(correspondence.assignments) == [(0, 1), (1, 0)]

    def test_zero_weight_dropped(self, read_text):
        """Check that a zero-weight beam is in no beam assignment."""
        correspondence = read_text('(1 0) = (1 0) + 0*(0 1)')
        assert correspondence.assignments == {(1, 0): {(1, 0): 1.0}}

    def test_issue_206(self, read_text):
        """Check the rotational-domain example of issue #206."""
        expect = BeamCorrespondence({
            (1, 0): {(1, 0): 1, (-1, 0): 1, (0, 1): 1, (0, -1): 1},
            (1, 1): {(1, 1): 1, (-1, -1): 1, (1, -1): 1, (-1, 1): 1},
            })
        correspondence = read_text(_ISSUE_206)
        assert correspondence == expect
        assert not correspondence.is_weighted

    def test_issue_206_weighted(self, read_text):
        """Check the weighted rotational-domain example of #206."""
        expect = BeamCorrespondence({
            (1, 0): {(1, 0): 1/3, (-1, 0): 1/3, (0, 1): 1/6, (0, -1): 1/6},
            (1, 1): {(1, 1): 1, (-1, -1): 1, (1, -1): 1, (-1, 1): 1},
            })
        correspondence = read_text(_ISSUE_206_WEIGHTED)
        assert correspondence == expect
        assert correspondence.is_weighted

    def test_path_str(self, tmp_path):
        """Check reading from a path given as a string."""
        path = tmp_path / 'BEAM_EQUIVALENCE'
        path.write_text('(1 0) = (1 0)', encoding='utf-8')
        expect = BeamCorrespondence({(1, 0): {(1, 0): 1}})
        assert read_beam_equivalence(str(path)) == expect

    @parametrize('contents,line_nr,cause', _INVALID_FILES.values(),
                 ids=_INVALID_FILES)
    def test_invalid(self, contents, line_nr, cause, read_text, tmp_path):
        """Check the error for an invalid line.

        Mutants 'duplicate left-hand side accepted' (no check for an
        experimental beam already in the file) and 'no snapping before
        the duplicate check' (keys compared before _snap_to_fraction)
        fail the 'duplicate left' and 'duplicate left snapped' cases.
        """
        path = tmp_path / 'BEAM_EQUIVALENCE'
        prefix = re.escape(f'{path}, line {line_nr}: ')
        with pytest.raises(ValueError, match=prefix) as exc_info:
            read_text(contents)
        message = str(exc_info.value)
        assert message.startswith(f'{path}, line {line_nr}: ')
        assert cause in message

    @parametrize(contents=('', '\n\n', '# only a comment\n  ! another\n'))
    def test_no_assignment(self, contents, read_text, tmp_path):
        """Check that a file without beam assignments is an error."""
        path = tmp_path / 'BEAM_EQUIVALENCE'
        with pytest.raises(ValueError, match='No beam assignment') as exc:
            read_text(contents)
        assert str(exc.value).startswith(str(path))
