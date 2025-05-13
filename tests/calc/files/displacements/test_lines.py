"""Tests for module files/displacements/lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'

import pytest

from viperleed_jax.files.displacements.tokens.direction import DirectionToken
from viperleed_jax.files.displacements.lines import GeoDeltaLine, separate_direction_from_targets
from viperleed_jax.files.displacements.tokens.range import RangeToken
from viperleed_jax.files.displacements.tokens.target import TargetingError, Targets
from viperleed_jax.files.displacements.errors import InvalidDisplacementsSyntaxError


@pytest.mark.parametrize(
    'line, exp_subtargets, exp_direction, exp_range',
    [
        pytest.param(
            'Fe z = -0.2 0.2',
            ['Fe'],
            DirectionToken('z'),
            RangeToken.from_floats(-0.2, 0.2),
            id='single-target-no-step',
        ),
        pytest.param(
            'Fe xy[1 1] = -0.2 0.2',
            ['Fe'],
            DirectionToken('xy[1 1]'),
            RangeToken.from_floats(-0.2, 0.2),
            id='single-target-no-step',
        ),
        pytest.param(
            'Fe 1 2, Cu_surf xy = -1.5 3.0 0.5',
            ['Fe 1 2', 'Cu_surf'],
            DirectionToken('xy'),
            RangeToken.from_floats(-1.5, 3.0, 0.5),
            id='multi-target-with-step',
        ),
        pytest.param(
            '  Fe_*    xyz=   0.    2E+1   0.25  ',
            ['Fe_*'],
            DirectionToken('xyz'),
            RangeToken.from_floats(0.0, 20.0, 0.25),
            id='whitespace-scientific',
        ),
        pytest.param(
            'Fe L(1-3), Ni_sub L(2) z = -0.5 0.5',
            ['Fe L(1-3)', 'Ni_sub L(2)'],
            DirectionToken('z'),
            RangeToken.from_floats(-0.5, 0.5),
            id='layer-selector',
        ),
    ],
)
def test_geodelta_parametrized(line, exp_subtargets, exp_direction, exp_range):
    geo = GeoDeltaLine(line)
    # targets
    actual_targets = [st.target_str for st in geo.targets.subtargets]
    assert actual_targets == exp_subtargets
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
def test_invalid_format_raises(bad_line, exp_error):
    with pytest.raises(exp_error):
        GeoDeltaLine(bad_line)


# def test_invalid_direction_raises(monkeypatch):
#     # simulate Direction throwing
#     from displacements.direction import Direction as RealDirection

#     def fake_init(self, s):
#         raise ValueError('bad dir')

#     monkeypatch.setattr('displacements.lines.Direction.__init__', fake_init)
#     with pytest.raises(ValueError) as excinfo:
#         GeoDeltaLine('A foo = 0 1')
#     assert 'Unable to parse direction' in str(excinfo.value)


# def test_invalid_range_raises(monkeypatch):
#     # simulate RangeToken throwing
#     from displacements.range import RangeToken as RealRange

#     def fake_init(self, s):
#         raise ValueError('bad range')

#     monkeypatch.setattr(
#         'displacements.lines.RangeToken.__init__', fake_init
#     )
#     with pytest.raises(ValueError) as excinfo:
#         GeoDeltaLine('A x = 0 1')
#     assert 'Bad range' in str(excinfo.value)



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
