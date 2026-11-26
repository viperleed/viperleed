"""Tests for module viperleed.calc.files.new_displacements.lines."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-04-10'
__license__ = 'GPLv3+'

import re

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.calc.files.new_displacements.errors import (
    DisplacementsSyntaxError,
)
from viperleed.calc.files.new_displacements.lines import (
    ConstraintLine,
    GeoDeltaLine,
    LoopMarkerLine,
    ModeHeaderLine,
    OccDeltaLine,
    OffsetsLine,
    VibDeltaLine,
    separate_direction_from_targets,
)
from viperleed.calc.files.new_displacements.tokens.direction import (
    CartesianDirectionToken,
    FractionalDirectionToken,
)
from viperleed.calc.files.new_displacements.tokens.element import ElementToken
from viperleed.calc.files.new_displacements.tokens.linear_operation import (
    LinearOperationToken,
)
from viperleed.calc.files.new_displacements.tokens.mode import ModeToken
from viperleed.calc.files.new_displacements.tokens.offset import OffsetToken
from viperleed.calc.files.new_displacements.tokens.range_ import RangeToken
from viperleed.calc.files.new_displacements.tokens.target import TargetToken
from viperleed.calc.files.new_displacements.tokens.total_occupation import (
    TotalOccupationToken,
)


class TestGeoDeltaLine:
    @parametrize(
        'line, exp_targets, exp_direction, exp_range',
        [
            pytest.param(
                'Fe z = -0.2 0.2',
                [TargetToken('Fe')],
                CartesianDirectionToken('z'),
                RangeToken.from_floats(-0.2, 0.2),
                id='single-target-no-step',
            ),
            pytest.param(
                'Fe xy[1 1] = -0.2 0.2',
                [TargetToken('Fe')],
                CartesianDirectionToken('xy[1 1]'),
                RangeToken.from_floats(-0.2, 0.2),
                id='single-target-no-step',
            ),
            pytest.param(
                'Fe 1 2, Cu_surf xy = -1.5 3.0 0.5',
                [TargetToken('Fe 1 2'), TargetToken('Cu_surf')],
                CartesianDirectionToken('xy'),
                RangeToken.from_floats(-1.5, 3.0, 0.5),
                id='multi-target-with-step',
            ),
            pytest.param(
                '  Fe_*    xyz=   0.    2E+1   0.25  ',
                [TargetToken('Fe_*')],
                CartesianDirectionToken('xyz'),
                RangeToken.from_floats(0.0, 20.0, 0.25),
                id='whitespace-scientific',
            ),
            pytest.param(
                'Fe L(1-3), Ni_sub L(2) z = -0.5 0.5',
                [TargetToken('Fe L(1-3)'), TargetToken('Ni_sub L(2)')],
                CartesianDirectionToken('z'),
                RangeToken.from_floats(-0.5, 0.5),
                id='layer-selector',
            ),
            pytest.param(
                'O 1 ab[1 2] = -0.1 0.1',
                [TargetToken('O 1')],
                FractionalDirectionToken('ab[1 2]'),
                RangeToken.from_floats(-0.1, 0.1),
                id='fractional-vector',
            ),
            pytest.param(
                'O 1 [0.5 1] = -0.1 0.1',
                [TargetToken('O 1')],
                FractionalDirectionToken('ab[1 2]'),
                RangeToken.from_floats(-0.1, 0.1),
                id='bare-fractional-vector',
            ),
            pytest.param(
                'O, Fe b = -0.1 0.1',
                [TargetToken('O'), TargetToken('Fe')],
                FractionalDirectionToken('ab[0 1]'),
                RangeToken.from_floats(-0.1, 0.1),
                id='fractional-b',
            ),
            pytest.param(
                'O ab = -0.1 0.1',
                [TargetToken('O')],
                FractionalDirectionToken('ab'),
                RangeToken.from_floats(-0.1, 0.1),
                id='fractional-plane',
            ),
        ],
    )
    def test_geodelta_parametrized(
        self, line, exp_targets, exp_direction, exp_range
    ):
        geo = GeoDeltaLine(line)
        # targets
        assert len(geo.targets) == len(exp_targets)
        for target, exp_target in zip(geo.targets, exp_targets):
            assert target == exp_target
        # direction
        assert isinstance(geo.direction, type(exp_direction))
        assert geo.direction == exp_direction
        # range
        assert isinstance(geo.range, RangeToken)
        assert geo.range == exp_range

    @parametrize(
        'bad_line, exp_error',
        [
            pytest.param('A =', DisplacementsSyntaxError, id='missing RHS'),
            pytest.param(
                'x = 0 1', DisplacementsSyntaxError, id='missing target'
            ),
            pytest.param(
                'A x 0 1', DisplacementsSyntaxError, id='not equal sign'
            ),
            pytest.param(
                'A y = foo bar',
                DisplacementsSyntaxError,
                id='non-numeric RHS',
            ),
            pytest.param('A c = 0 1', DisplacementsSyntaxError, id='c'),
            pytest.param('A abc = 0 1', DisplacementsSyntaxError, id='abc'),
            pytest.param('A az = 0 1', DisplacementsSyntaxError, id='az'),
            pytest.param(
                'A azi(ab[1 0]) = 0 1', DisplacementsSyntaxError, id='azi'
            ),
        ],
    )
    def test_invalid_format_raises(self, bad_line, exp_error):
        with pytest.raises(exp_error):
            GeoDeltaLine(bad_line)

    _invalid_cartesian = (
        'yx', 'xz', 'yz', 'zyx', 'xx[1 2]', 'xyx[1 0 5]', 'zz[1 1]',
        'zyx[1 0 0]', 'x[1]', 'z[-1]', 'x[0]', 'zz', 'xxy',
        'xy[0 0]', 'xy[1 0 0]',
    )

    @parametrize(direction=_invalid_cartesian)
    def test_invalid_cartesian_direction(self, direction):
        """Check that unordered or repeated x, y, z are rejected."""
        with pytest.raises(DisplacementsSyntaxError):
            GeoDeltaLine(f'O 1 {direction} = -0.1 0.1')


class TestVibDeltaLine:
    @parametrize(
        'line, exp_targets, exp_range',
        [
            # single target, no step
            ('A = 0 1', [TargetToken('A')], RangeToken.from_floats(0.0, 1.0)),
            # multiple targets comma separated, with step
            (
                'A1, B2 = -1.5 3.0 0.5',
                [TargetToken('A1'), TargetToken('B2')],
                RangeToken.from_floats(-1.5, 3.0, 0.5),
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
        assert not hasattr(vib, 'direction')
        # range
        assert isinstance(vib.range, RangeToken)
        assert vib.range == exp_range

    @parametrize(
        'bad_line, exp_error',
        [
            ('A x = 0 1', DisplacementsSyntaxError),
            ('A =', DisplacementsSyntaxError),
            ('A 0 1', DisplacementsSyntaxError),
            ('A = foo bar', DisplacementsSyntaxError),
            ('A = 0 1 = 2', DisplacementsSyntaxError),
        ],
    )
    def test_invalid(self, bad_line, exp_error):
        with pytest.raises(exp_error):
            VibDeltaLine(bad_line)

    def test_single_value_is_fixed_range(self):
        """A single VIB_DELTA value is read as a fixed range."""
        vib = VibDeltaLine('O 1 = 0.1')
        assert vib.range == VibDeltaLine('O 1 = 0.1 0.1').range
        assert vib.range == RangeToken.from_floats(0.1, 0.1)

    def test_single_value_not_numeric(self):
        """A non-numeric single VIB_DELTA value is rejected."""
        with pytest.raises(DisplacementsSyntaxError):
            VibDeltaLine('O 1 = foo')

    @parametrize(
        'value', ('nan', 'inf', '-inf', 'infinity', '1_0', '1e999', '0x1')
    )
    def test_single_value_not_plain_finite(self, value):
        """Non-finite and non-decimal single values are rejected."""
        with pytest.raises(DisplacementsSyntaxError):
            VibDeltaLine(f'O 1 = {value}')

    @parametrize('rhs', ('nan 1', '0 inf', '1_0 2', '0 1 nan'))
    def test_range_not_plain_finite(self, rhs):
        """Non-finite and non-decimal range values are rejected."""
        with pytest.raises(DisplacementsSyntaxError):
            VibDeltaLine(f'O 1 = {rhs}')

    def test_str_shows_targets_and_range(self):
        line = 'X,Y = 2 4 1'
        vib = VibDeltaLine(line)
        rep = str(vib)
        assert 'TargetToken' in rep
        assert '=' in rep
        assert 'RangeToken' in rep


@parametrize(
    'line_cls, line',
    (
        (GeoDeltaLine, 'O 0 z = -0.1 0.1'),
        (VibDeltaLine, 'O L(0) = -0.1 0.1'),
        (OccDeltaLine, 'O 0 = Fe 0 1'),
        (OffsetsLine, 'geo O L(0) z = 0.1'),
        (ConstraintLine, 'geo O 0 = O 1'),
    ),
)
def test_bad_target_number_in_message(line_cls, line):
    """Check that the error names the bad atom or layer number."""
    with pytest.raises(DisplacementsSyntaxError) as exc_info:
        line_cls(line)
    assert re.search(r'Invalid (atom|layer) number 0', exc_info.value.message)


def test_lowercase_layer_in_geo_line():
    """Check that a lowercase layer keeps the direction and range."""
    lower = GeoDeltaLine('O l(2) z = -0.1 0.1 0.05')
    upper = GeoDeltaLine('O L(2) z = -0.1 0.1 0.05')
    assert lower.direction == upper.direction
    assert lower.range == upper.range


@parametrize(
    'line_type, line, exp_offsets_line',
    (
        (GeoDeltaLine, 'O 1 z offset = 0.1', 'geo O 1 z = 0.1'),
        (GeoDeltaLine, 'O 1 L(2) z Offset = 0.1', 'geo O 1 L(2) z = 0.1'),
        (
            GeoDeltaLine,
            'O 1, N 2 xy offset = 0.1 0.2',
            'geo O 1, N 2 xy = 0.1 0.2',
        ),
        (GeoDeltaLine, 'O 1 z offset = clear', 'geo O 1 z = 0'),
        (VibDeltaLine, 'O 1 offset = 0.1', 'vib O 1 = 0.1'),
        (VibDeltaLine, 'O 1-2 L(1) offset = -0.05', 'vib O 1-2 L(1) = -0.05'),
    ),
)
def test_inline_offset_suggests_offsets_line(
    line_type, line, exp_offsets_line
):
    """Check that an in-line offset error quotes a valid line."""
    with pytest.raises(DisplacementsSyntaxError) as exc_info:
        line_type(line)
    message = str(exc_info.value)
    assert line in message
    assert f'"{exp_offsets_line}"' in message
    # The suggestion is itself a valid line in the OFFSETS block
    assert OffsetsLine(exp_offsets_line).raw_line == exp_offsets_line


_INLINE_OFFSET_CASES = {
    'geo_zero_no_direction': (GeoDeltaLine, 'O 1 offset = 0', 'Delete'),
    'geo_clear_no_direction': (GeoDeltaLine, 'O 1 offset = clear', 'Delete'),
    'geo_original_no_direction': (
        GeoDeltaLine, 'O 1 offset = original', 'Delete'
    ),
    'geo_original': (
        GeoDeltaLine, 'O 1 z offset = original', None, 'geo O 1 z = 0'
    ),
    'geo_range': (
        GeoDeltaLine,
        'O 1 z offset = -0.1 0.1 0.05',
        None,
        'geo O 1 z = -0.1',
    ),
    'geo_cont': (GeoDeltaLine, 'O 1 z offset = cont', 'Delete'),
    'geo_nonzero_no_direction': (
        GeoDeltaLine, 'O 1 offset = 0.1', 'needs a direction'
    ),
    'geo_multi_dof_one_value': (
        GeoDeltaLine, 'O 1 xy offset = 0.1', 'one value per component'
    ),
    'geo_suffix_zoffset': (
        GeoDeltaLine, 'O 1 zoffset = 0.1', None, 'geo O 1 z = 0.1'
    ),
    'geo_suffix_bracket': (
        GeoDeltaLine, 'O 1 [0 1]offset = 0.1', None, 'geo O 1 [0 1] = 0.1'
    ),
    'vib_range': (VibDeltaLine, 'O 1 offset = -0.05 0.05 0.01', 'VIB_DELTA'),
    'vib_cont': (VibDeltaLine, 'O 1 offset = cont', 'Delete'),
    'occ_element': (OccDeltaLine, 'O 1 offset = O 0.5', 'no element'),
    'occ_value': (OccDeltaLine, 'O 1 offset = 0.5', None, 'occ O 1 = 0.5'),
}


@parametrize('line_type, line, keyword, suggestion',
             [(*c, None)[:4] for c in _INLINE_OFFSET_CASES.values()],
             ids=list(_INLINE_OFFSET_CASES))
def test_inline_offset_advice(line_type, line, keyword, suggestion):
    """Check that quoted OFFSETS suggestions parse, else give advice."""
    with pytest.raises(DisplacementsSyntaxError) as exc_info:
        line_type(line)
    message = str(exc_info.value)
    assert line in message
    quoted = re.findall(r'Write "([^"]+)"', message)
    for offsets_line in quoted:
        assert OffsetsLine(offsets_line)
    if suggestion:
        assert quoted == [suggestion]
    else:
        assert not quoted
        assert keyword in message
    if keyword == 'VIB_DELTA':
        assert VibDeltaLine('O 1 = -0.05 0.05 0.01')


def test_geo_single_value_is_invalid():
    """Check that a GEO_DELTA line needs a range, not one value."""
    with pytest.raises(DisplacementsSyntaxError):
        GeoDeltaLine('O 1 z = 0.1')


class TestOffsetLine:
    @parametrize(
        'line, exp_targets, exp_direction, exp_offset',
        [
            # 1, 2 and 3D geometric offsets
            pytest.param(
                'geo A x = 1.23',
                ['A'],
                'x',
                1.23,
                id='geo-1D',
            ),
            pytest.param(
                'geo A xy= 1.23 4.56',
                ['A'],
                'xy',
                (1.23, 4.56),
                id='geo-2D',
            ),
            pytest.param(
                'geo A xyz = 1.23 4.56 7.89',
                ['A'],
                'xyz',
                (1.23, 4.56, 7.89),
                id='geo-3D',
            ),
            # vib offset, no direction
            pytest.param(
                'vib B = 2.0',
                ['B'],
                None,
                2.0,
                id='vib-no-direction',
            ),
            # occ offset, multiple targets, no direction
            pytest.param(
                'occ C1, D2 = -0.5',
                ['C1', 'D2'],
                None,
                -0.5,
                id='occ-multiple-targets',
            ),
        ],
    )
    def test_offsets_valid(self, line, exp_targets, exp_direction, exp_offset):
        off = OffsetsLine(line)
        # targets
        assert len(off.targets) == len(exp_targets)
        for tok, exp in zip(off.targets, exp_targets):
            assert isinstance(tok, TargetToken)
            assert tok.target_str == exp
        # direction
        if exp_direction is None:
            assert off.direction is None
        else:
            assert isinstance(off.direction, CartesianDirectionToken)
            assert off.direction == CartesianDirectionToken(exp_direction)
        # offset
        assert isinstance(off.offset, OffsetToken)
        assert off.offset == OffsetToken.from_floats(exp_offset)

    @parametrize(
        'bad_line, exp_error',
        [
            # Missing direction for geo
            ('geo A = 1.0', DisplacementsSyntaxError),
            # Direction not allowed for non-geo
            ('vib A x = 2.0', DisplacementsSyntaxError),
            # Missing type and target
            ('A = 1.0', DisplacementsSyntaxError),
            # No equals sign
            ('geo A x 1.0', DisplacementsSyntaxError),
            # Non-numeric offset
            ('geo A x = foo', DisplacementsSyntaxError),
            # Multiple offset tokens
            ('geo A x = 1 2', DisplacementsSyntaxError),
            # Multiple '=' signs
            ('geo A x = 1.0 = 2.0', DisplacementsSyntaxError),
        ],
    )
    def test_offsets_invalid(self, bad_line, exp_error):
        with pytest.raises(exp_error):
            OffsetsLine(bad_line)


class TestOccDeltaLine:
    @parametrize(
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
                    (
                        ElementToken('Au'),
                        RangeToken.from_floats(0.1, 1.0, 0.1),
                    ),
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
        assert len(occ.targets) == len(exp_targets)
        for target, exp_target in zip(occ.targets, exp_targets):
            assert target == exp_target
        # check element_ranges length and content
        assert hasattr(occ, 'element_ranges')
        assert len(occ.element_ranges) == len(exp_elem_ranges)
        for (act_e, act_r), (exp_e, exp_r) in zip(
            occ.element_ranges, exp_elem_ranges
        ):
            assert isinstance(act_e, ElementToken)
            assert act_e == exp_e
            assert isinstance(act_r, RangeToken)
            assert act_r == exp_r

    @parametrize(
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
        with pytest.raises(DisplacementsSyntaxError):
            OccDeltaLine(bad_line)

    def test_single_value_is_fixed_range(self):
        """A single OCC_DELTA value is read as a fixed range."""
        occ = OccDeltaLine('O 1 = O 0.8, N 0.2')
        explicit = OccDeltaLine('O 1 = O 0.8 0.8, N 0.2 0.2')
        assert occ.element_ranges == explicit.element_ranges
        assert occ.element_ranges[0][1] == RangeToken.from_floats(0.8, 0.8)

    @parametrize('value', ('1.5', '-0.1', '2', 'nan', 'inf', '1_0'))
    def test_single_value_invalid_occupation(self, value):
        """Occupations outside [0, 1] or non-finite are invalid."""
        with pytest.raises(DisplacementsSyntaxError):
            OccDeltaLine(f'O 1 = O {value}')

    def test_str_contains_elements_and_ranges(self):
        line = 'A,B = Fe 0.0 1.0, Ni 0.0 0.2'
        occ = OccDeltaLine(line)
        rep = str(occ)
        assert 'A' in rep and 'Fe' in rep and '0.0' in rep
        assert 'B' in rep and 'Ni' in rep and '0.2' in rep
        assert '=' in rep


class TestConstraintLine:
    @parametrize(
        'line, exp_type, exp_targets, exp_link_target, exp_op',
        [
            pytest.param(
                'geo A_surf, A_def = linked',
                ModeToken('geo'),
                [TargetToken('A_surf'), TargetToken('A_def')],
                TargetToken('A_surf'),
                LinearOperationToken.from_array(np.eye(1)),
                id='geo-linked-shorthand',
            ),
            pytest.param(
                'geo A = B',
                ModeToken('geo'),
                [TargetToken('A')],
                TargetToken('B'),
                LinearOperationToken.from_array(np.eye(1)),
                id='geo-direct-link',
            ),
            pytest.param(
                'geo Fe 1 = Fe 2',
                ModeToken('geo'),
                [TargetToken('Fe 1')],
                TargetToken('Fe 2'),
                LinearOperationToken.from_array(np.eye(1)),
                id='geo-single-target-link',
            ),
            pytest.param(
                'geo Fe_surf, O_surf = [[1 0 0] [0 0 1] [0 1 0]] Fe_surf',
                ModeToken('geo'),
                [TargetToken('Fe_surf'), TargetToken('O_surf')],
                TargetToken('Fe_surf'),
                LinearOperationToken('[[1 0 0] [0 0 1] [0 1 0]]'),
                id='geo-linear-operation',
            ),
            pytest.param(
                'occ A, B = linked',
                ModeToken('occ'),
                [TargetToken('A'), TargetToken('B')],
                TargetToken('A'),
                LinearOperationToken.from_array(np.eye(1)),
                id='occ-linked-shorthand',
            ),
            pytest.param(
                'occ A* = total 1.0',
                ModeToken('occ'),
                [TargetToken('A*')],
                None,
                TotalOccupationToken('1.0'),
                id='occ-total',
            ),
            pytest.param(
                'vib X, Y = 0.5 Z',
                ModeToken('vib'),
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
        assert con.mode_token == exp_type
        # Targets
        act_targets = con.targets
        if 'linked' in line:
            act_targets = con.targets[:-1]  # remove link target from list
        for target, exp_target in zip(act_targets, exp_targets):
            assert target == exp_target

        # Link target
        assert con.link_target == exp_link_target
        # Array
        assert con.linear_operation == exp_op

    @parametrize(
        'bad_line',
        [
            'geo = linked',  # no targets
            'vib A =',  # empty rhs
            'geo A = [1 0 0',  # malformed op
            'vib = [1 0] B',  # missing lhs target
        ],
    )
    def test_invalid_constraint_lines(self, bad_line):
        with pytest.raises(DisplacementsSyntaxError):
            ConstraintLine(bad_line)

    @parametrize(
        'line, exp_mode, exp_value, exp_element',
        (
            ('geo O 1 = 0.05', 'geo', 0.05, None),
            ('vib O 1 = 0.05', 'vib', 0.05, None),
            ('vib O 1 = -1', 'vib', -1.0, None),
            ('geo O 1 = -0.05', 'geo', -0.05, None),
            ('occ O 1 = O 0.8', 'occ', 0.8, ElementToken('O')),
            ('occ O 1, N 2 = Fe 0.25', 'occ', 0.25, ElementToken('Fe')),
        ),
    )
    def test_fixed_value(self, line, exp_mode, exp_value, exp_element):
        """Check the fixed value and element read from the rhs."""
        con = ConstraintLine(line)
        assert con.mode_token == ModeToken(exp_mode)
        assert con.fixed_value == pytest.approx(exp_value)
        assert con.fixed_element == exp_element
        assert con.link_target is None
        assert con.linear_operation is None

    @parametrize(
        'line',
        (
            'geo O 1 = 0.05 0.1',  # two values
            'occ O 1 = 0.8',  # element missing
            'occ O 1 = Xx 0.8',  # invalid element
            'occ O 1 = O 0.8 0.9',  # two values
        ),
    )
    def test_invalid_fixed_value(self, line):
        """Check that malformed fixed values are rejected."""
        with pytest.raises(DisplacementsSyntaxError):
            ConstraintLine(line)

    def test_link_has_no_fixed_value(self):
        """A link to another target is not a fixed value."""
        con = ConstraintLine('occ O 1 = O 2')
        assert con.fixed_value is None
        assert con.fixed_element is None
        assert con.link_target == TargetToken('O 2')

    def test_occ_integer_is_link(self):
        """An integer after the element is a link to an atom number."""
        con = ConstraintLine('occ O 1 = O 2')
        assert con.link_target == TargetToken('O 2')
        assert con.fixed_value is None
        assert con.fixed_element is None

    def test_occ_non_integer_freezes(self):
        """A non-integer after the element freezes the occupation."""
        con = ConstraintLine('occ O 1 = O 1.0')
        assert con.fixed_value == 1.0
        assert con.fixed_element == ElementToken('O')
        assert con.link_target is None

    @parametrize('line', ('occ O 1 = O 0', 'occ O 1 = N 0'))
    def test_occ_invalid_integer_asks_for_decimal(self, line):
        """An integer that is no valid target suggests a decimal."""
        with pytest.raises(DisplacementsSyntaxError, match=r'\.0'):
            ConstraintLine(line)

    @parametrize('mode', ('geo', 'vib', 'occ'))
    @parametrize(
        'value', ('nan', 'inf', '-inf', 'infinity', '1_0', '1e999', '0.5.1')
    )
    def test_fixed_value_not_plain_finite(self, mode, value):
        """Non-finite and non-decimal fixed values are rejected."""
        rhs = f'O {value}' if mode == 'occ' else value
        with pytest.raises(DisplacementsSyntaxError):
            ConstraintLine(f'{mode} O 1 = {rhs}')

    @parametrize('value', ('2.5', '-0.5', '1.0001'))
    def test_fixed_occupation_out_of_range(self, value):
        """Fixed occupations outside [0, 1] are rejected."""
        with pytest.raises(DisplacementsSyntaxError, match=r'\[0, 1\]'):
            ConstraintLine(f'occ O 1 = O {value}')

    @parametrize(
        'line, exp_str',
        (
            ('geo O 1 = 0', '= 0.0)'),
            ('vib O 1 = 0.05', '= 0.05'),
            ('occ O 1 = O 0.0', '= ElementToken'),
            ('geo O 1 = ind(3)', 'ind(3)'),
            ('geo O 1 = O 2', 'TargetToken'),
        ),
    )
    def test_str_of_value_constraints(self, line, exp_str):
        """The string form shows each kind of right-hand side."""
        text = str(ConstraintLine(line))
        assert exp_str in text
        assert 'None' not in text

    @parametrize('mode', ('geo', 'vib', 'occ'))
    @parametrize('rhs', ('ind(2)', 'IND(2)', 'Ind( 2 )'))
    def test_fixed_index(self, mode, rhs):
        """Check that ind(<n>) is read, in any case and spacing."""
        con = ConstraintLine(f'{mode} O 1 = {rhs}')
        assert con.fixed_index == 2
        assert con.link_target is None
        assert con.linear_operation is None
        assert con.fixed_value is None
        assert con.fixed_element is None
        assert 'ind(2)' in str(con)

    @parametrize(
        'rhs',
        ('ind()', 'ind(x)', 'ind(0)', 'ind(-1)', 'ind(1.5)', 'ind(2) 3'),
    )
    def test_invalid_fixed_index(self, rhs):
        """Check that the error for a bad index quotes the line."""
        line = f'geo O 1 = {rhs}'
        with pytest.raises(DisplacementsSyntaxError, match=re.escape(line)):
            ConstraintLine(line)

    @parametrize(
        'line', ('geo O 1 = O 2', 'occ O 1 = O 0.8', 'vib O 1 = 0.1')
    )
    def test_other_forms_have_no_fixed_index(self, line):
        """Check that only ind(<n>) sets a fixed index."""
        assert ConstraintLine(line).fixed_index is None

    @parametrize(
        'line',
        ('occ O 1 = O 0.8 Fe 0.2', 'geo O 1 = 0.1 0.2'),
    )
    def test_invalid_format_lists_all_forms(self, line):
        """The hint lists the link, value, ind(<n>) and total forms."""
        with pytest.raises(DisplacementsSyntaxError) as exc_info:
            ConstraintLine(line)
        message = exc_info.value.message
        assert line in message
        for form in ('<linear_operation>', '<value>', '<element> <value>',
                     'ind(<n>)', 'total <total_occupation>'):
            assert form in message

    def test_offset_in_constrain_line(self):
        with pytest.raises(DisplacementsSyntaxError,
                           match=r'do not support "offset".*OFFSETS block'):
            ConstraintLine('geo A = offset 0.1')


@parametrize(
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
        # Fractional directions
        ('T1 2 ab', 'T1 2', 'ab'),
        ('T1 L(1) ab[1 2]', 'T1 L(1)', 'ab[1 2]'),
        ('T1 [1 2]', 'T1', '[1 2]'),
        ('T1 [0.5 -1]', 'T1', '[0.5 -1]'),
        # Invalid frames are still split off (and rejected later)
        ('T1 bc', 'T1', 'bc'),
        ('T1 xb', 'T1', 'xb'),
        # 3D vector
        ('T2 xyz[1 2 -1]', 'T2', 'xyz[1 2 -1]'),
        # azi(...) form
        ('Foo azi(ab[1 2])', 'Foo', 'azi(ab[1 2])'),
        # r(...) form
        ('Bar r([1 0])', 'Bar', 'r([1 0])'),
        # Extra whitespace
        ('  A1,B2   x  ', 'A1,B2', 'x'),
        # Empty input
        ('', '', ''),
    ],
)
def test_separate_direction(input_str, exp_targets, exp_direction):
    targets, direction = separate_direction_from_targets(input_str)
    assert targets == exp_targets
    assert direction == exp_direction


@parametrize(header=('= CONSTRAINT', '=constraint', '= Constrain'))
def test_constraint_header_alias(header):
    """Check that CONSTRAINT is read as the CONSTRAIN section."""
    assert ModeHeaderLine(header).section == 'CONSTRAIN'


_DIRECTION_FIRST = (
    pytest.param(
        GeoDeltaLine,
        'Cu_def xy L(2-4) = -0.10 0.10 0.01',
        'Cu_def L(2-4) xy = -0.10 0.10 0.01',
        id='geo_layer',
    ),
    pytest.param(
        GeoDeltaLine,
        'O z 1 2 = -0.1 0.1',
        'O 1 2 z = -0.1 0.1',
        id='geo_atoms',
    ),
    pytest.param(
        OffsetsLine,
        'geo O ab[1 0] L(1) = 0.1',
        'geo O L(1) ab[1 0] = 0.1',
        id='offset_geo',
    ),
)


@parametrize(
    'line_cls, line, expected',
    _DIRECTION_FIRST,
)
def test_direction_before_targets_message(line_cls, line, expected):
    """Check that a direction before targets shows the fixed line."""
    with pytest.raises(DisplacementsSyntaxError) as exc_info:
        line_cls(line)
    assert f'Write "{expected}"' in exc_info.value.message
    assert f'"{line}"' in exc_info.value.message


@parametrize(
    'line',
    (
        'O 1 xy[1 0] = -0.1 0.1 0.05',
        'O L(1) z = -0.1 0.1 0.05',
        'O 1 2 ab[1 0] = -0.1 0.1 0.05',
        'Cu_def L(2-4) xy = -0.10 0.10 0.01',
    ),
)
def test_direction_after_targets_ok(line):
    """Check that lines with the direction last are still read."""
    GeoDeltaLine(line)


class TestLoopMarkerLine:
    """Tests for the LoopMarkerLine class."""

    _valid = (
        ('<loop>', 'start'),
        ('  <LOOP>  ', 'start'),
        ('< loop >', 'start'),
        ('</loop>', 'end'),
        ('<\\loop>', 'end'),
        ('< / loop >', 'end'),
        ('</Loop>', 'end'),
    )
    _invalid = (
        '<loop>junk',
        '<loop> GEO_DELTA',
        '</loop> foo',
        '<\\loop>x',
        'x <loop>',
        '<lo op>',
        '<loops>',
        'loop',
    )

    @parametrize('line,kind', _valid)
    def test_valid(self, line, kind):
        """Check that a plain or spaced marker is read."""
        marker = LoopMarkerLine(line)
        assert marker.kind == kind

    @parametrize(line=_invalid)
    def test_invalid(self, line):
        """Check that a marker with extra characters is rejected."""
        with pytest.raises(DisplacementsSyntaxError):
            LoopMarkerLine(line)
