"""Tests for module viperleed.calc.files.new_displacements.file."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__license__ = 'GPLv3+'


import pytest
from pytest_cases import parametrize

from viperleed.calc.files.new_displacements.errors import (
    AlreadyReadError,
)
from viperleed.calc.files.new_displacements.file import DisplacementsFile
from viperleed.calc.files.new_displacements.segments import LoopBlock
from viperleed.calc.files.new_displacements.segments import SearchBlock

from ....helpers import TEST_DATA

NESTED_LOOPS_PATH = (
    TEST_DATA / 'DISPLACEMENTS/Cu_111/DISPLACEMENTS_nested_loops'
)


def test_read_from_file(displacements_file_path, subtests):
    """Test reading a file and checking its validity."""

    df = DisplacementsFile()
    assert df._has_been_read is False
    df.read(displacements_file_path)
    assert df._has_been_read is True

    with subtests.test('read() called again'), pytest.raises(AlreadyReadError):
        # check that read() cannot be called again
        df.read(displacements_file_path)


def test_read_nested_loops():
    """Check that <loop> markers in a file give nested LoopBlocks."""
    df = DisplacementsFile()
    df.read(NESTED_LOOPS_PATH)

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


_NESTED_LOOPS_SEQUENCES = {
    # Each loop runs twice: the second pass does not improve R.
    'constant R': ((0.5,)*9, 'abbcabbcd'),
    # The inner loop repeats as long as R improves. When the outer
    # loop repeats, the inner one starts afresh and runs twice.
    'improving R': ((0.5, 0.4, 0.3) + (0.2,)*7, 'abbbcabbcd'),
    # R improves during the outer loop's second pass: a third pass.
    'outer improves': ((0.5, 0.5, 0.4, 0.4, 0.4) + (0.3,)*8, 'abbcabbcabbcd'),
}


@parametrize(
    'rfacs,expected',
    _NESTED_LOOPS_SEQUENCES.values(),
    ids=_NESTED_LOOPS_SEQUENCES,
)
def test_next_nested_loops_sequence(rfacs, expected):
    """Check the search blocks next() returns for nested loops."""
    df = DisplacementsFile()
    df.read(NESTED_LOOPS_PATH)
    labels = ''.join(df.next(current_rfac=rfac).label for rfac in rfacs)
    assert labels == expected
    with pytest.raises(StopIteration):
        df.next(current_rfac=rfacs[-1])


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
