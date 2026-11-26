"""Tests for module viperleed.calc.files.new_displacements.file."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__license__ = 'GPLv3+'


import pytest

from viperleed.calc.files.new_displacements.errors import (
    AlreadyReadError,
)
from viperleed.calc.files.new_displacements.file import DisplacementsFile
from viperleed.calc.files.new_displacements.segments import LoopBlock
from viperleed.calc.files.new_displacements.segments import SearchBlock


def test_read_from_file(displacements_file_path, subtests):
    """Test reading a file and checking its validity."""

    df = DisplacementsFile()
    assert df._has_been_read is False
    df.read(displacements_file_path)
    assert df._has_been_read is True

    with subtests.test('read() called again'), pytest.raises(AlreadyReadError):
        # check that read() cannot be called again
        df.read(displacements_file_path)


def test_read_nested_loops(tmp_path):
    """Check that <loop> markers read from a file give nested LoopBlocks."""
    path = tmp_path / 'DISPLACEMENTS_nested_loops'
    path.write_text(
        '<loop>\n'
        '== SEARCH a\n'
        '= GEO_DELTA\n'
        'Fe 1 z = 0.0 0.1\n'
        '<loop>\n'
        '== SEARCH b\n'
        '= GEO_DELTA\n'
        'Fe 1 z = 0.0 0.1\n'
        '</loop>\n'
        '== SEARCH c\n'
        '= GEO_DELTA\n'
        'Fe 1 z = 0.0 0.1\n'
        '</loop>\n'
        '== SEARCH d\n'
        '= GEO_DELTA\n'
        'Fe 1 z = 0.0 0.1\n'
    )
    df = DisplacementsFile()
    df.read(path)

    outer, last = df.children
    assert isinstance(outer, LoopBlock)
    assert isinstance(last, SearchBlock)
    assert last.label == 'd'
    first, inner, third = outer.children
    assert isinstance(first, SearchBlock)
    assert isinstance(inner, LoopBlock)
    assert isinstance(third, SearchBlock)
    assert [c.label for c in (first, third)] == ['a', 'c']
    assert [type(c) for c in inner.children] == [SearchBlock]
    assert inner.children[0].label == 'b'


@pytest.mark.xfail(
    reason='Fractional directions not yet supported', strict=False
)
def test_read_from_unsupported_file(displacements_file_path_xfail):
    """Test reading a file and checking its validity."""

    # try and read the unsupported file, which should raise an error
    df = DisplacementsFile()
    df.read(displacements_file_path_xfail)

def test_next(displacements_file_path):
    """Test the next() method of DisplacementsFile."""

    df = DisplacementsFile()
    df.read(displacements_file_path)

    # try calling next() until we reach StopIteration
    with pytest.raises(StopIteration):
        while True:
            line = df.next(current_rfac=0.1)
            assert line is not None


def _exhaust(displacements_file):
    """Return all search blocks that next() yields until it stops."""
    blocks = []
    with pytest.raises(StopIteration):
        while True:
            blocks.append(displacements_file.next(current_rfac=0.1))
    return blocks


def test_reset_iterator(displacements_file_path):
    """Check that reset_iterator restarts a finished file."""
    df = DisplacementsFile()
    df.read(displacements_file_path)
    first_pass = _exhaust(df)
    df.reset_iterator()
    assert _exhaust(df) == first_pass
