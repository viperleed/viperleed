"""Tests for the DisplacementsFile class and its components."""

__authors__ = ('Alexander M. Imre (@amimre)',)

import pytest

from viperleed_jax.files.displacements.errors import (
    InvalidSearchBlocksError,
    InvalidSearchLoopError,
    OffsetsNotAtBeginningError,
)
from viperleed_jax.files.displacements.file import (
    DisplacementsFile,
    OffsetsBlock,
    SearchBlock,
)
from viperleed_jax.files.displacements.lines import (
    GeoDeltaLine,
    LoopMarkerLine,
    OffsetsLine,
)
from viperleed_jax.files.displacements.reader import (
    DisplacementFileSections,
    LoopMarker,
)


# Mock backend
class DummyBackend:
    def __init__(self, name='Dummy'):
        self.name = name

    def can_handle_search_block(self, offsets_block, block):
        return True, None


class TestDispalacementsFileSyntax:
    def test_searchblock_add_line_and_access(self):
        sb = SearchBlock('test')
        line = GeoDeltaLine('A_surf z = -0.1 0.1')
        sb.add_line(DisplacementFileSections.GEO_DELTA, line)

        assert sb.label == 'test'
        assert sb.geo_delta == [line]
        assert sb.vib_delta == []
        assert sb.occ_delta == []
        assert sb.explicit_constraints == []
        assert 'test' in repr(sb)


    def test_searchblock_add_line_invalid_section(self):
        sb = SearchBlock('bad')
        with pytest.raises(ValueError):
            sb.add_line('INVALID_SECTION', GeoDeltaLine('A_surf z = -0.1 0.1'))


    def test_offsets_block_add_line(self):
        ob = OffsetsBlock()
        line = OffsetsLine('vib A_surf = 0.1')
        ob.add_line(line)
        assert ob.lines == [line]


    def test_offset_block_detection(self):
        df = DisplacementsFile()
        offset_block = OffsetsBlock()
        df.blocks = [offset_block]
        assert df.offsets_block() is offset_block

        df.blocks = [SearchBlock('A')]
        assert df.offsets_block() is None


    def test_displacementsfile_loop_marker_detection(self):
        df = DisplacementsFile()
        df.blocks = [LoopMarkerLine(LoopMarker.LOOP_START)]
        assert df.unclosed_loop() is True

        df.blocks.append(LoopMarkerLine(LoopMarker.LOOP_END))
        assert df.unclosed_loop() is False


    def test_displacementsfile_check_valid_blocks_and_loops(self):
        df = DisplacementsFile()
        with pytest.raises(InvalidSearchBlocksError):
            df.check_valid()

        df.blocks = [LoopMarkerLine(LoopMarker.LOOP_START)]
        with pytest.raises(InvalidSearchLoopError):
            df.check_valid()

        df.blocks = [OffsetsBlock(), SearchBlock('main')]
        assert df.check_valid() is True

        df.blocks = [SearchBlock('main'), OffsetsBlock()]
        with pytest.raises(OffsetsNotAtBeginningError):
            df.check_valid()


    def test_displacementsfile_first_block_excludes_loopmarkers(self):
        df = DisplacementsFile()
        block = SearchBlock('first')
        df.blocks = [
            LoopMarkerLine(LoopMarker.LOOP_START),
            block,
            LoopMarkerLine(LoopMarker.LOOP_END),
        ]
        assert df.first_block() is block


    def test_displacementsfile_parse_replaces_blocks(self):
        df = DisplacementsFile()
        block = SearchBlock('main')
        df.blocks = [block]
        df.has_been_read = True
        df.has_been_parsed = False
        df.check_valid = lambda: True  # patch to bypass
        df.parse(DummyBackend())
        assert df.has_been_parsed is True


def test_parse_already_called():
    df = DisplacementsFile()
    df.has_been_parsed = True
    with pytest.raises(ValueError):
        df.parse(DummyBackend())


def test_read_from_file(mock_displacements_cu_111_realistic_path_and_lines):
    """Test reading a file and checking its validity."""
    path, _ = (
        mock_displacements_cu_111_realistic_path_and_lines
    )

    df = DisplacementsFile()
    assert df.has_been_read is False
    df.read(path)
    assert df.has_been_read is True
    assert df.check_valid()

    with pytest.raises(ValueError):
        # check that read() cannot be called again
        df.read(path)
