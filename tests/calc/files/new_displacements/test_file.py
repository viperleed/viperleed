"""Tests for module viperleed.calc.files.new_displacements.file."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)

import pytest

from viperleed.calc.files.new_displacements.errors import (
    AlreadyReadError,
)
from viperleed.calc.files.new_displacements.file import DisplacementsFile


def test_read_from_file(displacements_file_path, subtests):
    """Test reading a file and checking its validity."""

    df = DisplacementsFile()
    assert df._has_been_read is False
    df.read(displacements_file_path)
    assert df._has_been_read is True

    with subtests.test('read() called again'), pytest.raises(AlreadyReadError):
        # check that read() cannot be called again
        df.read(displacements_file_path)


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
