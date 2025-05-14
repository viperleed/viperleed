from pathlib import Path

import pytest
from pytest_cases import fixture

from viperleed_jax.files.displacements.lines import (
    ConstraintLine,
    GeoDeltaLine,
    LoopMarkerLine,
    OccDeltaLine,
    SearchHeaderLine,
    SectionHeaderLine,
    VibDeltaLine,
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
from viperleed_jax.files.displacements.reader import DisplacementsReader

_MOCK_DISPLACEMENTS_PATH = Path('tests/test_data/displacements/')
MOCK_DISPLACEMENTS_PATH = Path(
    'tests/test_data/displacements/DISPLACEMENTS_mixed'
)
_CU_111_DISPLACEMENTS_PATH = _MOCK_DISPLACEMENTS_PATH / 'Cu_111'
_CU_111_SIMPLE_PATH = _CU_111_DISPLACEMENTS_PATH / 'DISPLACEMENTS_simple'


@pytest.fixture()
def mock_displacements():
    path = MOCK_DISPLACEMENTS_PATH
    expected_lines = [
        (SearchHeaderLine, {'label': 'my_search_label'}),
        (SectionHeaderLine, {'section': 'GEO_DELTA'}),
        (
            GeoDeltaLine,
            {
                'targets': (TargetToken('O L(1-2)'),),
                'direction': DirectionToken('z'),
                'range': RangeToken('-0.05 0.05 0.005'),
            },
        ),
        (
            GeoDeltaLine,
            {
                'targets': (TargetToken('Ir L(1)'),),
                'direction': DirectionToken('z'),
                'range': RangeToken('0.03 -0.03 0.003'),
            },
        ),
        (SectionHeaderLine, {'section': 'VIB_DELTA'}),
        (
            VibDeltaLine,
            {
                'targets': (TargetToken('O 1'),),
                'range': RangeToken('-0.05 0.05 0.02'),
            },
        ),
        (
            VibDeltaLine,
            {
                'targets': (TargetToken('Ir_top'),),
                'range': RangeToken('-0.05 0.05 0.01'),
            },
        ),
        (
            VibDeltaLine,
            {
                'targets': (TargetToken('Si* 2'),),
                'range': RangeToken('0.05 -0.05'),
            },
        ),
        (SectionHeaderLine, {'section': 'OCC_DELTA'}),
        (
            OccDeltaLine,
            {
                'targets': (TargetToken('M_top'),),
                'element_ranges': (
                    (ElementToken('Fe'), RangeToken('0.4 0.6 0.05')),
                    (ElementToken('Ni'), RangeToken('0.6 0.4 0.05')),
                ),
            },
        ),
        (
            OccDeltaLine,
            {
                'targets': (TargetToken('O 1'),),
                'element_ranges': ((ElementToken('O'), RangeToken('0.8 1.0')),),
            },
        ),
        (SectionHeaderLine, {'section': 'CONSTRAIN'}),
        (
            ConstraintLine,
            {
                'type': TypeToken('vib'),
                'targets': (TargetToken('Ir_top'),),
                'linear_operation': LinearOperationToken.from_array([[1.0]]),
                'link_target': TargetToken('Ir_top'),
            },
        ),
        (
            ConstraintLine,
            {
                'type': TypeToken('geo'),
                'targets': (TargetToken('O L(1-2)'), TargetToken('Ir L(1)')),
                'linear_operation': LinearOperationToken.from_array([[1.0]]),
                'link_target': TargetToken('O L(1-2)'),
            },
        ),
    ]
    return path, expected_lines


# @fixture
# def mock_displacements_simple():
#     path = _CU_111_SIMPLE_PATH
#     expected_lines = [
#         SearchHeaderLine('test'),
#         SectionHeaderLine('GEO_DELTA'),
#         GeoDeltaLine('Cu_surf', None, 'z', -0.1, 0.1, 0.05,
#                      'Cu_surf z = -0.1 0.1 0.05'),
#         SectionHeaderLine('VIB_DELTA'),
#         VibDeltaLine('Cu_surf', None, -0.1, 0.1, 0.05,
#                      'Cu_surf = -0.1 0.1 0.05'),
#     ]
#     return path, expected_lines


def test_displacements_reader(mock_displacements):
    path, expected_lines = mock_displacements
    with DisplacementsReader(path) as reader:
        parsed_lines = list(reader)
    assert len(parsed_lines) == len(expected_lines)
    for parsed, expected in zip(parsed_lines, expected_lines):
        # check correct type
        assert isinstance(parsed, expected[0])
        # check correct attributes
        for attr, value in expected[1].items():
            assert getattr(parsed, attr) == value


# def test_displacements_reader_simple(mock_displacements_simple):
#     path, expected_lines = mock_displacements_simple
#     with DisplacementsReader(path) as reader:
#         parsed_lines = list(reader)
#     assert len(parsed_lines) == len(expected_lines)
#     for parsed, expected in zip(parsed_lines, expected_lines):
#         assert parsed == expected
