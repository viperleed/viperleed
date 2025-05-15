"""Conftest for files/displacements."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2025-05-15'

from pathlib import Path

import pytest
from pytest_cases import fixture

from viperleed_jax.files.displacements.lines import (
    ConstraintLine,
    GeoDeltaLine,
    OccDeltaLine,
    SearchHeaderLine,
    SectionHeaderLine,
    VibDeltaLine,
)
from viperleed_jax.files.displacements.tokens import (
    DirectionToken,
    ElementToken,
    LinearOperationToken,
    RangeToken,
    TargetToken,
    TypeToken,
)

_MOCK_DISPLACEMENTS_PATH = Path('tests/test_data/displacements/')
MOCK_DISPLACEMENTS_PATH = Path(
    'tests/test_data/displacements/DISPLACEMENTS_mixed'
)
_CU_111_DISPLACEMENTS_PATH = _MOCK_DISPLACEMENTS_PATH / 'Cu_111'
_CU_111_SIMPLE_PATH = _CU_111_DISPLACEMENTS_PATH / 'DISPLACEMENTS_simple'

@pytest.fixture
def mock_displacements_path_and_lines():
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


@fixture
def mock_displacements_cu_111_realistic_path_and_lines():
    path = _CU_111_SIMPLE_PATH
    expected_lines = [
        (SearchHeaderLine, {'label': 'test'}),
        (SectionHeaderLine, {'section': 'GEO_DELTA'}),
        (
            GeoDeltaLine,
            {
                'targets': (TargetToken('Cu_surf'),),
                'direction': DirectionToken('z'),
                'range': RangeToken('-0.1 0.1 0.05'),
            },
        ),
        (SectionHeaderLine, {'section': 'VIB_DELTA'}),
        (
            VibDeltaLine,
            {
                'targets': (TargetToken('Cu_surf'),),
                'range': RangeToken('-0.1 0.1 0.05'),
            },
        ),
    ]
    return path, expected_lines
