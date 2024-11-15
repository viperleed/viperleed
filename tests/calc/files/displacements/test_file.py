from pathlib import Path

import pytest

from viperleed_jax.files.displacements.file import (
    DisplacementsFile,
    SearchBlock,
)
from viperleed_jax.files.displacements.lines import (
    ConstraintLine,
    GeoDeltaLine,
    LoopMarkerLine,
    OccDeltaLine,
    VibDeltaLine,
)
from viperleed_jax.files.displacements.reader import DisplacementFileSections

MOCK_DISPLACEMENTS_PATH = Path(
    'tests/test_data/displacements/DISPLACEMENTS_mixed'
)


@pytest.fixture()
def mock_displacements_file():
    path = MOCK_DISPLACEMENTS_PATH
    expected_blocks = [SearchBlock('my_search_label')]

    # Expected block structure
    expected_blocks[0].add_line(
        DisplacementFileSections.GEO_DELTA,
        GeoDeltaLine('O', 'L(1-2)', 'z', -0.05, 0.05, 0.005),
    )
    expected_blocks[0].add_line(
        DisplacementFileSections.GEO_DELTA,
        GeoDeltaLine('Ir', 'L(1)', 'z', 0.03, -0.03, 0.003),
    )
    expected_blocks[0].add_line(
        DisplacementFileSections.VIB_DELTA,
        VibDeltaLine('O', '1', -0.05, 0.05, 0.02),
    )
    expected_blocks[0].add_line(
        DisplacementFileSections.VIB_DELTA,
        VibDeltaLine('Ir_top', None, -0.05, 0.05, 0.01),
    )
    expected_blocks[0].add_line(
        DisplacementFileSections.VIB_DELTA,
        VibDeltaLine('Si', None, 0.1, None, None),
    )
    expected_blocks[0].add_line(
        DisplacementFileSections.OCC_DELTA,
        OccDeltaLine(
            'M_top', None, [('Fe', 0.4, 0.6, 0.05), ('Ni', 0.6, 0.4, 0.05)]
        ),
    )
    expected_blocks[0].add_line(
        DisplacementFileSections.OCC_DELTA,
        OccDeltaLine('O', '1', [('O', 0.8, None, None)]),
    )
    expected_blocks[0].add_line(
        DisplacementFileSections.CONSTRAIN,
        ConstraintLine('vib', ['Ir_top'], 'linked'),
    )
    expected_blocks[0].add_line(
        DisplacementFileSections.CONSTRAIN,
        ConstraintLine('geo', ['O L(1-2)', 'Ir L(1)'], 'linked'),
    )

    return path, expected_blocks


def test_displacements_file_read(mock_displacements_file):
    """Test the DisplacementsFile read method using a multi-line file input."""
    path, expected_blocks = mock_displacements_file

    # Read the file and parse
    displacements_file = DisplacementsFile()
    displacements_file.read(path)

    # Check that the expected search blocks are present
    assert len(displacements_file.blocks) == len(expected_blocks)

    for parsed_block, expected_block in zip(
        displacements_file.blocks, expected_blocks
    ):
        assert parsed_block.label == expected_block.label
        assert parsed_block.sections == expected_block.sections
