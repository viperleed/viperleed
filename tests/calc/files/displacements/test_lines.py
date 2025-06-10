"""Tests for module files/displacements/lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'

import numpy as np
import pytest

from viperleed_jax.files.displacements.errors import (
    InvalidDisplacementsSyntaxError,
)
from viperleed_jax.files.displacements.lines import (
    ConstraintLine,
    GeoDeltaLine,
    OccDeltaLine,
    OffsetsLine,
    VibDeltaLine,
    separate_direction_from_targets,
)
from viperleed_jax.files.displacements.tokens import (
    DirectionToken,
    ElementToken,
    LinearOperationToken,
    OffsetToken,
    RangeToken,
    TargetToken,
    TokenParserError,
    TypeToken,
)
from viperleed_jax.files.displacements.tokens.direction import DirectionToken


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


class TestOccDeltaLine:
    @pytest.mark.parametrize(
        'line, exp_targets, exp_elem_ranges',
        [
            # single target, single element-range
            (
                'A = H 0.0 1.0',
                [TargetToken('A')],
                [(ElementToken('H'), RangeToken.from_floats(0.0, 1.0))],
            ),
            # multiple targets
            (
                'A,B = Ni 0 1',
                [TargetToken('A'), TargetToken('B')],
                [(ElementToken('Ni'), RangeToken.from_floats(0.0, 1.0))],
            ),
            # two element-range pairs
            (
                'X = O 0.1 1, N 0 0.5',
                [TargetToken('X')],
                [
                    (ElementToken('O'), RangeToken.from_floats(0.1, 1.0)),
                    (ElementToken('N'), RangeToken.from_floats(0.0, 0.5)),
                ],
            ),
            # multiple targets and multiple pairs with step
            (
                'A_surf, B_def L(2-3) = Pt 0.1 0.2, Au 0.1 1 0.1',
                [TargetToken('A_surf'), TargetToken('B_def L(2-3)')],
                [
                    (ElementToken('Pt'), RangeToken.from_floats(0.1, 0.2)),
                    (ElementToken('Au'), RangeToken.from_floats(0.1, 1.0, 0.1)),
                ],
            ),
            # whitespace tolerance
            (
                '  D   =   Li  0.5 1.0  ,  Be 0.3 0.5 ',
                [TargetToken('D')],
                [
                    (ElementToken('Li'), RangeToken.from_floats(0.5, 1.0)),
                    (ElementToken('Be'), RangeToken.from_floats(0.3, 0.5)),
                ],
            ),
        ],
    )
    def test_occ_delta_valid(self, line, exp_targets, exp_elem_ranges):
        occ = OccDeltaLine(line)
        # check targets
        for target, exp_target in zip(occ.targets, exp_targets):
            assert target == exp_target
        # check element_ranges length and content
        assert hasattr(occ, 'element_ranges')
        assert len(occ.element_ranges) == len(exp_elem_ranges)
        for (act_e, act_r), (exp_e, exp_r) in zip(
            occ.element_ranges, exp_elem_ranges
        ):
            assert isinstance(act_e, ElementToken)
            assert act_e.atomic_number == exp_e.atomic_number
            assert isinstance(act_r, RangeToken)
            assert act_r == exp_r


    @pytest.mark.parametrize(
        'bad_line',
        [
            '',  # empty line
            'A H 1 2',  # missing '='
            'A x = H 1 2',  # direction on LHS
            'A =',  # missing RHS
            'A = H',  # missing range
            'A = H1 2',  # no space between element and range
            'A = H foo bar',  # non-numeric range
            'A = H 1 2 = N 0 1',  # multiple '='
            'A = Xx 0 1',  # invalid element
        ],
    )
    def test_occ_delta_invalid(self, bad_line):
        with pytest.raises(InvalidDisplacementsSyntaxError):
            OccDeltaLine(bad_line)


    def test_repr_contains_elements_and_ranges(self):
        line = 'A,B = Fe 0.0 1.0, Ni 0.0 0.2'
        occ = OccDeltaLine(line)
        rep = repr(occ)
        assert 'A' in rep and 'Fe' in rep and '0.0' in rep
        assert 'B' in rep and 'Ni' in rep and '0.2' in rep
        assert '=' in rep


class TestConstraintLine:
    @pytest.mark.parametrize(
        'line, exp_type, exp_targets, exp_link_target, exp_op',
        [
            pytest.param(
                'geo A_surf, A_def = linked',
                TypeToken('geo'),
                [TargetToken('A_surf'), TargetToken('A_def')],
                TargetToken('A_surf'),
                LinearOperationToken.from_array(np.eye(1)),
                id='geo-linked-shorthand',
            ),
            pytest.param(
                'geo A = B',
                TypeToken('geo'),
                [TargetToken('A')],
                TargetToken('B'),
                LinearOperationToken.from_array(np.eye(1)),
                id='geo-direct-link',
            ),
            pytest.param(
                'geo Fe 1 = Fe 2',
                TypeToken('geo'),
                [TargetToken('Fe 1')],
                TargetToken('Fe 2'),
                LinearOperationToken.from_array(np.eye(1)),
                id='geo-single-target-link',
            ),
            pytest.param(
                'geo A = [1 0 0] B',
                TypeToken('geo'),
                [TargetToken('A')],
                TargetToken('B'),
                LinearOperationToken('[1 0 0]'),
                id='geo-linear-operation',
            ),
            pytest.param(
                'geo Fe_surf, O_surf = [[1 0 0] [0 0 1] [0 1 0]] Fe_surf',
                TypeToken('geo'),
                [TargetToken('Fe_surf'), TargetToken('O_surf')],
                TargetToken('Fe_surf'),
                LinearOperationToken('[[1 0 0] [0 0 1] [0 1 0]]'),
                id='geo-linear-operation',
            ),
            pytest.param(
                'occ A, B = linked',
                TypeToken('occ'),
                [TargetToken('A'), TargetToken('B')],
                TargetToken('A'),
                LinearOperationToken.from_array(np.eye(1)),
                id='occ-linked-shorthand',
            ),
            pytest.param(
                'vib X, Y = 0.5 Z',
                TypeToken('vib'),
                [TargetToken('X'), TargetToken('Y')],
                TargetToken('Z'),
                LinearOperationToken.from_array(np.array(0.5)),
                id='vib-float-link',
            ),
        ],
    )
    def test_valid_constraint_lines(
        self, line, exp_type, exp_targets, exp_link_target, exp_op
    ):
        con = ConstraintLine(line)
        # Type
        assert con.type == exp_type
        # Targets
        if 'linked' in line:
            for target, exp_target in zip(con.targets[:-1], exp_targets):
                assert target == exp_target
        else:
            for target, exp_target in zip(con.targets, exp_targets):
                assert target == exp_target
        # Link target
        assert con.link_target == exp_link_target
        # Array
        assert con.linear_operation == exp_op

    @pytest.mark.parametrize(
        'bad_line',
        [
            'geo = linked',  # no targets
            'vib A =',  # empty rhs
            'geo A = [1 0 0',  # malformed op
            'vib = [1 0] B',  # missing lhs target
        ],
    )
    def test_invalid_constraint_lines(self, bad_line):
        with pytest.raises(InvalidDisplacementsSyntaxError):
            ConstraintLine(bad_line)


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
