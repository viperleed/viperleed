"""Tests for module files/displacements/lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-04-10'

import pytest

from viperleed_jax.files.displacements.direction import Direction
from viperleed_jax.files.displacements.lines import GeoDeltaLine
from viperleed_jax.files.displacements.range import DisplacementsRange
from viperleed_jax.files.displacements.targeting import TargetingError, Targets


@pytest.mark.parametrize(
    'line, exp_subtargets, exp_direction, exp_range',
    [
        pytest.param(
            'Fe z = -0.2 0.2',
            ['Fe'],
            Direction('z'),
            DisplacementsRange.from_floats(-0.2, 0.2),
            id='single-target-no-step',
        ),
        pytest.param(
            'Fe xy[1 1] = -0.2 0.2',
            ['Fe'],
            Direction('xy[1 1]'),
            DisplacementsRange.from_floats(-0.2, 0.2),
            id='single-target-no-step',
        ),
        pytest.param(
            'Fe 1 2, Cu_surf xy = -1.5 3.0 0.5',
            ['Fe 1 2', 'Cu_surf'],
            Direction('xy'),
            DisplacementsRange.from_floats(-1.5, 3.0, 0.5),
            id='multi-target-with-step',
        ),
        pytest.param(
            '  Fe_*    xyz=   0.    2E+1   0.25  ',
            ['Fe_*'],
            Direction('xyz'),
            DisplacementsRange.from_floats(0.0, 20.0, 0.25),
            id='whitespace-scientific',
        ),
        pytest.param(
            'Fe L(1-3), Ni_sub L(2) z = -0.5 0.5',
            ['Fe L(1-3)', 'Ni_sub L(2)'],
            Direction('z'),
            DisplacementsRange.from_floats(-0.5, 0.5),
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
    assert isinstance(geo.direction, Direction)
    assert geo.direction == exp_direction
    # range
    assert isinstance(geo.range, DisplacementsRange)
    assert geo.range == exp_range




# @pytest.mark.parametrize(
#     'bad_line',
#     [
#         'A =',  # missing RHS
#         'x = 0 1',  # missing target
#         'A = 0',  # only one number on RHS
#         'A x 0 1',  # missing '='
#         'A y = foo bar',  # non-numeric RHS
#     ],
# )
# def test_invalid_format_raises(bad_line):
#     with pytest.raises(ValueError):
#         GeoDeltaLine(bad_line)


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
#     # simulate DisplacementsRange throwing
#     from displacements.range import DisplacementsRange as RealRange

#     def fake_init(self, s):
#         raise ValueError('bad range')

#     monkeypatch.setattr(
#         'displacements.lines.DisplacementsRange.__init__', fake_init
#     )
#     with pytest.raises(ValueError) as excinfo:
#         GeoDeltaLine('A x = 0 1')
#     assert 'Bad range' in str(excinfo.value)