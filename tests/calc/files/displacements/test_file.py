"""Tests for the DisplacementsFile class and its components."""

__authors__ = ('Alexander M. Imre (@amimre)',)

import pytest

from viperleed.calc.files.new_displacements.errors import (
    InvalidSearchBlocksError,
    InvalidSearchLoopError,
    OffsetsNotAtBeginningError,
)
from viperleed.calc.files.new_displacements.file import DisplacementsFile

def test_read_from_file(displacements_file_path, subtests):
    """Test reading a file and checking its validity."""

    df = DisplacementsFile()
    assert df.has_been_read is False
    df.read(displacements_file_path)
    assert df.has_been_read is True

    with subtests.test("read() called again"):
        with pytest.raises(ValueError):
            # check that read() cannot be called again
            df.read(displacements_file_path)
