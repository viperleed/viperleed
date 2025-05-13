"""Tests for module files/displacements/lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'

import pytest

from viperleed_jax.files.displacements.errors import (
    InvalidDisplacementsSyntaxError,
)
from viperleed_jax.files.displacements.tokens import (
    DirectionToken,
    ElementToken,
    OffsetToken,
    RangeToken,
    TargetToken,
    TokenParserError,
    TypeToken,
)

from viperleed_jax.files.displacements.tokens.direction import DirectionToken
from viperleed_jax.files.displacements.lines import (
    GeoDeltaLine,
    VibDeltaLine,
    OffsetsLine,
    separate_direction_from_targets
)

class TestGeoDeltaLine:
    @pytest.mark.parametrize(
        'line, exp_targets, exp_direction, exp_range',
        [
            pytest.param(
                'Fe z = -0.2 0.2',
                [TargetToken('Fe')],
                DirectionToken('z'),
                RangeToken.from_floats(-0.2, 0.2),
                id='single-target-no-step',
            ),
            pytest.param(
                'Fe xy[1 1] = -0.2 0.2',
                [TargetToken('Fe')],
                DirectionToken('xy[1 1]'),
                RangeToken.from_floats(-0.2, 0.2),
                id='single-target-no-step',
            ),
            pytest.param(
                'Fe 1 2, Cu_surf xy = -1.5 3.0 0.5',
                [TargetToken('Fe 1 2'), TargetToken('Cu_surf')],
                DirectionToken('xy'),
                RangeToken.from_floats(-1.5, 3.0, 0.5),
                id='multi-target-with-step',
            ),
            pytest.param(
                '  Fe_*    xyz=   0.    2E+1   0.25  ',
                [TargetToken('Fe_*')],
                DirectionToken('xyz'),
                RangeToken.from_floats(0.0, 20.0, 0.25),
                id='whitespace-scientific',
            ),
            pytest.param(
                'Fe L(1-3), Ni_sub L(2) z = -0.5 0.5',
                [TargetToken('Fe L(1-3)'), TargetToken('Ni_sub L(2)')],
                DirectionToken('z'),
                RangeToken.from_floats(-0.5, 0.5),
                id='layer-selector',
            ),
        ],
    )
    def test_geodelta_parametrized(
        self, line,exp_targets, exp_direction, exp_range):
        geo = GeoDeltaLine(line)
        # targets
        assert len(geo.targets) == len(exp_targets)
        for target, exp_target in zip(geo.targets, exp_targets):
            assert target == exp_target
        # direction
        assert isinstance(geo.direction, DirectionToken)
        assert geo.direction == exp_direction
        # range
        assert isinstance(geo.range, RangeToken)
        assert geo.range == exp_range

    @pytest.mark.parametrize(
        'bad_line, exp_error',
        [
            pytest.param('A =', InvalidDisplacementsSyntaxError, id='missing RHS'),
            pytest.param(
                'x = 0 1', InvalidDisplacementsSyntaxError, id='missing target'
            ),
            pytest.param(
                'A x 0 1', InvalidDisplacementsSyntaxError, id='not equal sign'
            ),
            pytest.param(
                'A y = foo bar',
                InvalidDisplacementsSyntaxError,
                id='non-numeric RHS',
            ),
        ],
    )
    def test_invalid_format_raises(self, bad_line, exp_error):
        with pytest.raises(exp_error):
            GeoDeltaLine(bad_line)

class TestVibDeltaLine:

    @pytest.mark.parametrize(
        "line, exp_targets, exp_range",
        [
            # single target, no step
            ("A = 0 1",
             [TargetToken("A")],
             RangeToken.from_floats(0.0, 1.0)
            ),
            # multiple targets comma separated, with step
            ("A1, B2 = -1.5 3.0 0.5",
             [TargetToken("A1"), TargetToken("B2")],
             RangeToken.from_floats(-1.5, 3.0, 0.5)
            ),
        ],
    )
    def test_valid(self, line, exp_targets, exp_range):
        vib = VibDeltaLine(line)
        # targets
        assert len(vib.targets) == len(exp_targets)
        for actual, expected in zip(vib.targets, exp_targets):
            assert actual == expected
        # no direction allowed
        assert not hasattr(vib, "direction")
        # range
        assert isinstance(vib.range, RangeToken)
        assert vib.range == exp_range


    @pytest.mark.parametrize(
        "bad_line, exp_error",
        [
            ("A x = 0 1", InvalidDisplacementsSyntaxError),
            ("A =", InvalidDisplacementsSyntaxError),
            ("A 0 1", InvalidDisplacementsSyntaxError),
            ("A = foo bar", InvalidDisplacementsSyntaxError),
            ("A = 0 1 = 2", InvalidDisplacementsSyntaxError),
        ],
    )
    def test_invalid(self, bad_line, exp_error):
        with pytest.raises(exp_error):
            VibDeltaLine(bad_line)

    def test_repr_shows_targets_and_range(self):
        line = "X,Y = 2 4 1"
        vib = VibDeltaLine(line)
        rep = repr(vib)
        assert "TargetToken" in rep
        assert "=" in rep
        assert "RangeToken" in rep


class TestOffsetLine:

    @pytest.mark.parametrize(
        'line, exp_type, exp_targets, exp_direction, exp_offset',
        [
            # geometric offset with direction
            pytest.param(
                'geo A x = 1.23',
                'geo',
                ['A'],
                'x',
                1.23,
                id='geo-with-direction',
            ),
            # vib offset, no direction
            pytest.param(
                'vib B = 2.0',
                'vib',
                ['B'],
                None,
                2.0,
                id='vib-no-direction',
            ),
            # # occ offset, multiple targets, no direction
            # pytest.param(
            #     'occ C1, D2 = -0.5',
            #     'occ',
            #     ['C1', 'D2'],
            #     None,
            #     -0.5,
            #     id='occ-multiple-targets',
            # ),
        ],
    )
    def test_offsets_valid(self, line, exp_type, exp_targets, exp_direction, exp_offset):
        off = OffsetsLine(line)
        # type
        assert isinstance(off.type, TypeToken)
        assert off.type == TypeToken(exp_type)
        # targets
        assert len(off.targets) == len(exp_targets)
        for tok, exp in zip(off.targets, exp_targets):
            assert isinstance(tok, TargetToken)
            assert tok.target_str == exp
        # direction
        if exp_direction is None:
            assert off.direction is None
        else:
            assert isinstance(off.direction, DirectionToken)
            assert off.direction == DirectionToken(exp_direction)
        # offset
        assert isinstance(off.offset, OffsetToken)
        assert off.offset == OffsetToken.from_floats(exp_offset)


    @pytest.mark.parametrize(
        'bad_line, exp_error',
        [
            # Missing direction for geo
            ('geo A = 1.0', InvalidDisplacementsSyntaxError),
            # Direction not allowed for non-geo
            ('vib A x = 2.0', InvalidDisplacementsSyntaxError),
            # Missing type and target
            ('A = 1.0', InvalidDisplacementsSyntaxError),
            # No equals sign
            ('geo A x 1.0', InvalidDisplacementsSyntaxError),
            # Non-numeric offset
            ('geo A x = foo', InvalidDisplacementsSyntaxError),
            # Multiple offset tokens
            ('geo A x = 1 2', InvalidDisplacementsSyntaxError),
            # Multiple '=' signs
            ('geo A x = 1.0 = 2.0', InvalidDisplacementsSyntaxError),
        ],
    )
    def test_offsets_invalid(self, bad_line, exp_error):
        with pytest.raises(exp_error):
            OffsetsLine(bad_line)



@pytest.mark.parametrize(
    'input_str, exp_targets, exp_direction',
    [
        # No direction at all
        ('A1 B2', 'A1 B2', ''),
        # Single-letter direction
        ('A x', 'A', 'x'),
        # Multi-letter direction
        ('Label xy', 'Label', 'xy'),
        # Bracketed vector direction
        ('T1 xy[1 0]', 'T1', 'xy[1 0]'),
        # 3D vector
        ('T2 xyz[1 2 -1]', 'T2', 'xyz[1 2 -1]'),
        # azi(...) form
        ('Foo azi(ab[1 2])', 'Foo', 'azi(ab[1 2])'),
        # r(...) form
        ('Bar r([1 0])', 'Bar', 'r([1 0])'),
        # Extra whitespace
        ('  A1,B2   x  ', 'A1,B2', 'x'),
        # Only direction, no targets
        ('xy', '', 'xy'),
        ('[1 0]', '', '[1 0]'),
        # Empty input
        ('', '', ''),
    ],
)
def test_separate_direction(input_str, exp_targets, exp_direction):
    targets, direction = separate_direction_from_targets(input_str)
    assert targets == exp_targets
    assert direction == exp_direction
