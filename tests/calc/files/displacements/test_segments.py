"""Test for segments."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-06-18"
__license__ = "GPLv3+"


from viperleed.calc.files.new_displacements.lines import (
    ConstraintLine,
    GeoDeltaLine,
    LoopMarkerLine,
    OccDeltaLine,
    OffsetsLine,
    SearchHeaderLine,
    SectionHeaderLine,
    VibDeltaLine,
)
from viperleed.calc.files.new_displacements.segments import (
    LoopBlock,
    SearchBlock,
    )

def test_search_block():

    header = SearchHeaderLine("== SEARCH test")
    lines = iter(
        [
            SectionHeaderLine("= GEO_DELTA"),
            GeoDeltaLine("Fe 1 z = 0.0 0.1"),
            SectionHeaderLine("= CONSTRAIN"),
            ConstraintLine("geo Fe 1 = Fe 2"),
            SectionHeaderLine("= VIB_DELTA"),
            VibDeltaLine("Fe 1 = 0.0 0.1"),
        ]
    )
    search_block = SearchBlock(header)
    assert search_block.label == "test"

    returned_lines = search_block.read_lines(lines)
    returned_lines = list(returned_lines)
    # make sure all lines were processed
    assert len(returned_lines) == 0

def test_layered_loop_block():
    header = LoopMarkerLine('<loop>')
    lines = iter(
        [
            LoopMarkerLine("<loop>"),
            SearchHeaderLine("== SEARCH test"),
            SectionHeaderLine("= GEO_DELTA"),
            GeoDeltaLine("Fe 1 z = 0.0 0.1"),
            LoopMarkerLine("</loop>"),
            LoopMarkerLine("</loop>"),
        ]
    )
    loop_block = LoopBlock(header)

    returned_lines = loop_block.read_lines(lines)
    returned_lines = list(returned_lines)
    # make sure all lines were processed
    assert len(returned_lines) == 0
