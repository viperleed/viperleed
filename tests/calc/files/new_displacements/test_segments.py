"""Test for module viperleed.calc.files.new_displacements.segments."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-06-18'
__license__ = 'GPLv3+'

import pytest

from viperleed.calc.files.new_displacements.errors import (
    DisplacementsSyntaxError,
)
from viperleed.calc.files.new_displacements.lines import (
    ConstraintLine,
    GeoDeltaLine,
    LoopMarkerLine,
    ModeHeaderLine,
    SearchHeaderLine,
    VibDeltaLine,
)
from viperleed.calc.files.new_displacements.segments import (
    LoopBlock,
    SearchBlock,
)


def test_search_block():
    header = SearchHeaderLine('== SEARCH test')
    lines = iter(
        [
            ModeHeaderLine('= GEO_DELTA'),
            GeoDeltaLine('Fe 1 z = 0.0 0.1'),
            ModeHeaderLine('= CONSTRAIN'),
            ConstraintLine('geo Fe 1 = Fe 2'),
            ModeHeaderLine('= VIB_DELTA'),
            VibDeltaLine('Fe 1 = 0.0 0.1'),
        ]
    )
    search_block = SearchBlock(header)
    assert search_block.label == 'test'

    returned_lines = search_block.read_lines(lines)
    returned_lines = list(returned_lines)
    # make sure all lines were processed
    assert len(returned_lines) == 0


def test_layered_loop_block():
    header = LoopMarkerLine('<loop>')
    lines = iter(
        [
            LoopMarkerLine('<loop>'),
            SearchHeaderLine('== SEARCH z'),
            ModeHeaderLine('= GEO_DELTA'),
            GeoDeltaLine('Fe 1 z = 0.0 0.1'),
            LoopMarkerLine('</loop>'),
            SearchHeaderLine('== SEARCH y'),
            ModeHeaderLine('= GEO_DELTA'),
            GeoDeltaLine('Fe 1 xy[0 1] = 0.0 0.1'),
            LoopMarkerLine('</loop>'),
        ]
    )
    loop_block = LoopBlock(header)

    returned_lines = loop_block.read_lines(lines)
    returned_lines = list(returned_lines)
    # make sure all lines were processed
    assert len(returned_lines) == 0


def test_empty_loop_fails():
    header = LoopMarkerLine('<loop>')
    lines = iter(
        [
            LoopMarkerLine('</loop>'),
        ]
    )
    loop_block = LoopBlock(header)

    with pytest.raises(
        DisplacementsSyntaxError,
        match='Loop block cannot be empty.',
    ):
        loop_block.read_lines(lines)


def test_loop_belongs_to_me_with_non_marker_last_line():
    # a non-marker last line must not raise AttributeError on '.kind'
    loop_block = LoopBlock(LoopMarkerLine('<loop>'))
    loop_block._lines.append(SearchHeaderLine('== SEARCH x'))
    assert loop_block._belongs_to_me(LoopMarkerLine('</loop>'))


def test_loop_with_only_loop_child_complains():
    # test that a loop whose only child is another loop raises an error
    header = LoopMarkerLine('<loop>')
    lines = iter(
        [
            LoopMarkerLine('<loop>'),
            SearchHeaderLine('== SEARCH test'),
            ModeHeaderLine('= GEO_DELTA'),
            GeoDeltaLine('Fe 1 z = 0.0 0.1'),
            LoopMarkerLine('</loop>'),
            LoopMarkerLine('</loop>'),
        ]
    )
    loop_block = LoopBlock(header)

    with pytest.raises(
        DisplacementsSyntaxError,
        match='Loop block cannot contain only another loop block',
    ):
        loop_block.read_lines(lines)


def test_loop_reset_iterator():
    """Check that reset_iterator restarts a finished loop."""
    header = LoopMarkerLine('<loop>')
    lines = iter(
        [
            SearchHeaderLine('== SEARCH first'),
            ModeHeaderLine('= GEO_DELTA'),
            GeoDeltaLine('Fe 1 z = 0.0 0.1'),
            SearchHeaderLine('== SEARCH second'),
            ModeHeaderLine('= VIB_DELTA'),
            VibDeltaLine('Fe 1 = 0.0 0.1'),
            LoopMarkerLine('</loop>'),
        ]
    )
    loop_block = LoopBlock(header)
    loop_block.read_lines(lines)
    first, second = loop_block.children
    # The first pass improves the R factor, so the loop repeats
    # once, then stops when the R factor stays the same.
    expected = (first, second, first, second)
    returned = tuple(loop_block.next(0.5) for _ in expected)
    assert returned == expected
    with pytest.raises(StopIteration):
        loop_block.next(0.5)

    # After a reset, the loop behaves as if freshly read: it forgets
    # the R factor of the earlier passes and repeats once more.
    loop_block.reset_iterator()
    returned = tuple(loop_block.next(0.5) for _ in expected)
    assert returned == expected
    with pytest.raises(StopIteration):
        loop_block.next(0.5)
