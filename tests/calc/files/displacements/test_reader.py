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
        SearchHeaderLine('my_search_label'),
        SectionHeaderLine('GEO_DELTA'),
        GeoDeltaLine('O', 'L(1-2)', 'z', -0.05, 0.05, 0.005),
        GeoDeltaLine('Ir', 'L(1)', 'z', 0.03, -0.03, 0.003),
        SectionHeaderLine('VIB_DELTA'),
        VibDeltaLine('O', '1', -0.05, 0.05, 0.02),
        VibDeltaLine('Ir_top', None, -0.05, 0.05, 0.01),
        VibDeltaLine('Si*', '2', 0.05, -0.05, None),
        SectionHeaderLine('OCC_DELTA'),
        OccDeltaLine(
            'M_top', None, [('Fe', 0.4, 0.6, 0.05), ('Ni', 0.6, 0.4, 0.05)]
        ),
        OccDeltaLine('O', '1', [('O', 0.8, 1.0, None)]),
        SectionHeaderLine('CONSTRAIN'),
        ConstraintLine('vib', 'Ir_top', None, 'linked'),
        ConstraintLine('geo', 'O L(1-2), Ir L(1)', 'z', 'linked'),
    ]
    return path, expected_lines


@fixture
def mock_displacements_simple():
    path = _CU_111_SIMPLE_PATH
    expected_lines = [
        SearchHeaderLine('test'),
        SectionHeaderLine('GEO_DELTA'),
        GeoDeltaLine('Cu_surf', None, 'z', -0.1, 0.1, 0.05,
                     'Cu_surf z = -0.1 0.1 0.05'),
        SectionHeaderLine('VIB_DELTA'),
        VibDeltaLine('Cu_surf', None, -0.1, 0.1, 0.05,
                     'Cu_surf = -0.1 0.1 0.05'),
    ]
    return path, expected_lines


def test_displacements_reader(mock_displacements):
    path, expected_lines = mock_displacements
    with DisplacementsReader(path) as reader:
        parsed_lines = list(reader)
    assert len(parsed_lines) == len(expected_lines)
    for parsed, expected in zip(parsed_lines, expected_lines):
        assert parsed == expected


def test_displacements_reader_simple(mock_displacements_simple):
    path, expected_lines = mock_displacements_simple
    with DisplacementsReader(path) as reader:
        parsed_lines = list(reader)
    assert len(parsed_lines) == len(expected_lines)
    for parsed, expected in zip(parsed_lines, expected_lines):
        assert parsed == expected
