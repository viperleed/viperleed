"""Test configuration for calc/files/new_displacements."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-05-15'
__license__ = 'GPLv3+'

from operator import attrgetter

import pytest
from pytest_cases import fixture

from viperleed.calc.files.new_displacements.lines import (
    ConstraintLine,
    GeoDeltaLine,
    OccDeltaLine,
    SearchHeaderLine,
    SectionHeaderLine,
    VibDeltaLine,
)
from viperleed.calc.files.new_displacements.tokens import (
    CartesianDirectionToken,
    ElementToken,
    LinearOperationToken,
    ModeToken,
    RangeToken,
    TargetToken,
)

from ....helpers import TEST_DATA

_MOCK_DISPLACEMENTS_PATH = TEST_DATA / 'DISPLACEMENTS'
MOCK_DISPLACEMENTS_PATH = TEST_DATA / 'DISPLACEMENTS/DISPLACEMENTS_mixed'
_CU_111_DISPLACEMENTS_PATH = _MOCK_DISPLACEMENTS_PATH / 'Cu_111'
_CU_111_SIMPLE_PATH = _CU_111_DISPLACEMENTS_PATH / 'DISPLACEMENTS_simple'

DISPLACEMENTS_FILES = tuple(
    _MOCK_DISPLACEMENTS_PATH.glob('*/DISPLACEMENTS*')
)
_SUPPORTED_DISPLACEMENTS_FILES = tuple(
    p for p in DISPLACEMENTS_FILES if 'xfail' not in p.name
)
_UNSUPPORTED_DISPLACEMENTS_FILES = tuple(
    p for p in DISPLACEMENTS_FILES if 'xfail' in p.name
)

def displacements_file_id(file_path):
    """Return a case ID for a DISPLACEMENTS file."""
    rel_path = file_path.relative_to(_MOCK_DISPLACEMENTS_PATH).as_posix()
    return str(rel_path)


@fixture
@pytest.mark.parametrize(
    'file_path',
    _SUPPORTED_DISPLACEMENTS_FILES,
    ids=displacements_file_id,
)
def displacements_file_path(file_path):
    """Fixture to provide a path to a DISPLACEMENTS file."""
    return file_path

@fixture
@pytest.mark.parametrize(
    'file_path',
    _UNSUPPORTED_DISPLACEMENTS_FILES,
    ids=displacements_file_id,
)
def displacements_file_path_xfail(file_path):
    """Fixture to provide a path to a DISPLACEMENTS file."""
    return file_path


@fixture
def mock_displacements_path_and_lines():
    path = _MOCK_DISPLACEMENTS_PATH / 'DISPLACEMENTS_mixed'
    expected_lines = [
        (SearchHeaderLine, {'label': 'my_search_label'}),
        (SectionHeaderLine, {'section': 'GEO_DELTA'}),
        (
            GeoDeltaLine,
            {
                'targets': (TargetToken('O L(1-2)'),),
                'direction': CartesianDirectionToken('z'),
                'range': RangeToken('-0.05 0.05 0.005'),
            },
        ),
        (
            GeoDeltaLine,
            {
                'targets': (TargetToken('Ir L(1)'),),
                'direction': CartesianDirectionToken('z'),
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
                'targets': (TargetToken('Pt_sub'),),
                'range': RangeToken('-0.07 0.07 0.005'),
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
                'element_ranges': (
                    (ElementToken('O'), RangeToken('0.8 1.0')),
                ),
            },
        ),
        (SectionHeaderLine, {'section': 'CONSTRAIN'}),
        (
            ConstraintLine,
            {
                'mode_token': ModeToken('vib'),
                'targets': (TargetToken('Ir_top'),),
                'linear_operation': LinearOperationToken.from_array([[1.0]]),
                'link_target': TargetToken('Ir_top'),
            },
        ),
        (
            ConstraintLine,
            {
                'mode_token': ModeToken('geo'),
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
                'direction': CartesianDirectionToken('z'),
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
