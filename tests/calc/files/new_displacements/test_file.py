"""Tests for module viperleed.calc.files.new_displacements.file."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__license__ = 'GPLv3+'


import pytest
from pytest_cases import parametrize

from viperleed.calc.files.new_displacements.errors import (
    AlreadyReadError,
    OffsetsNotAtBeginningError,
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


_OFFSETS = '== OFFSETS\n    geo Cu_surf z = 0.05\n'
_SEARCH = '== SEARCH a\n    = GEO_DELTA\n    Cu_surf z = -0.1 0.1 0.05\n'


@parametrize(
    'contents, msg',
    (
        (_OFFSETS + _OFFSETS + _SEARCH,
         'There can only be one OFFSETS block at the beginning'),
        (_SEARCH + _OFFSETS,
         'The OFFSETS block must be at the beginning of the file'),
    ),
)
def test_read_offsets_misplaced(contents, msg, tmp_path):
    """Check the error for a duplicated or late OFFSETS block."""
    path = tmp_path / 'DISPLACEMENTS'
    path.write_text(contents)
    with pytest.raises(OffsetsNotAtBeginningError, match=msg):
        DisplacementsFile().read(path)


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


_NESTED_LOOPS_SEQUENCES = (
    # Each loop runs twice: the second pass does not improve R.
    pytest.param((0.5,)*9, 'abbcabbcd', id='constant R'),
    # The inner loop repeats as long as R improves. When the outer
    # loop repeats, the inner one starts afresh and runs twice.
    pytest.param((0.5, 0.4, 0.3) + (0.2,)*7, 'abbbcabbcd', id='improving R'),
    # R improves during the outer loop's second pass: a third pass.
    pytest.param((0.5, 0.5, 0.4, 0.4, 0.4) + (0.3,)*8,
                 'abbcabbcabbcd',
                 id='outer improves'),
)


@parametrize('rfacs,expected', _NESTED_LOOPS_SEQUENCES)
def test_next_nested_loops_sequence(rfacs, expected):
    """Check the search blocks next() returns for nested loops."""
    disp_file = DisplacementsFile()
    disp_file.read(NESTED_LOOPS_PATH)
    labels = ''.join(disp_file.next(current_rfac=rfac).label
                     for rfac in rfacs)
    assert labels == expected
    with pytest.raises(StopIteration):
        disp_file.next(current_rfac=rfacs[-1])


@pytest.mark.xfail(
    reason='Fractional directions not yet supported', strict=False
)
def test_read_from_unsupported_file(displacements_file_path_xfail):
    """Test reading a file and checking its validity."""

    # try and read the unsupported file, which should raise an error
    df = DisplacementsFile()
    df.read(displacements_file_path_xfail)


def _exhaust(displacements_file):
    """Return all search blocks that next() yields until it stops."""
    blocks = []
    while True:
        try:
            blocks.append(displacements_file.next(current_rfac=0.1))
        except StopIteration:
            return blocks


def test_next(displacements_file_path):
    """Test the next() method of DisplacementsFile."""
    df = DisplacementsFile()
    df.read(displacements_file_path)
    assert all(block is not None for block in _exhaust(df))


def test_reset_iterator(displacements_file_path):
    """Check that reset_iterator restarts a finished file."""
    df = DisplacementsFile()
    df.read(displacements_file_path)
    first_pass = _exhaust(df)
    df.reset_iterator()
    assert _exhaust(df) == first_pass
