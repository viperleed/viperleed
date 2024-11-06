"""Tests for module viperleed.calc.files.iorfactor."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2024 ViPErLEED developers"
__created__ = "2024-11-06"
__license__ = "GPLv3+"

import pytest
from io import StringIO
from viperleed.calc.files.input_reader import (
    InputStreamReader,
    InputFileReader,
    ShouldSkipLineError,
)


class MockInputStreamReader(InputStreamReader):
    """
    Mock implementation of InputStreamReader for testing purposes.

    This class provides a concrete implementation of the abstract
    `_read_one_line` method from InputReader to facilitate testing of the
    InputStreamReader functionality.
    """

    def __init__(self, source, noisy=True):
        super().__init__(source, noisy)

    def _read_one_line(self, line):
        if "skip" in line:
            raise ShouldSkipLineError("Line contains 'skip'")
        return line.strip()


class MockInputFileReader(InputFileReader):
    """
    Mock implementation of InputFileReader for testing purposes.

    This class provides a concrete implementation of the abstract
    `_read_one_line` method from InputReader to facilitate testing of the
    InputFileReader functionality.
    """

    def __init__(self, filename, noisy=True):
        super().__init__(filename, noisy)

    def _read_one_line(self, line):
        if "skip" in line:
            raise ShouldSkipLineError("Line contains 'skip'")
        return line.strip()


class TestInputStreamReader:
    """Tests for the InputStreamReader class."""

    def test_valid_input(self):
        source = StringIO("line1\nline2\nline3\n")
        reader = MockInputStreamReader(source)
        lines = [line for line in reader]
        assert lines == ["line1", "line2", "line3"]

    def test_invalid_input(self):
        with pytest.raises(TypeError):
            MockInputStreamReader("not a stream")

    def test_skip_line(self):
        source = StringIO("line1\nskip this line\nline3\n")
        reader = MockInputStreamReader(source)
        lines = [line for line in reader]
        assert lines == ["line1", "line3"]


class TestInputFileReader:
    """Tests for the InputFileReader class."""

    def test_file_not_found(self, tmp_path):
        non_existent_file = tmp_path / "non_existent.txt"
        with pytest.raises(FileNotFoundError):
            MockInputFileReader(non_existent_file)

    def test_context_manager(self, tmp_path):
        # Create a temporary test file
        test_file = tmp_path / "test.txt"
        test_file.write_text("line1\nline2\nline3\n")

        with MockInputFileReader(test_file) as reader:
            lines = [line for line in reader]
        assert lines == ["line1", "line2", "line3"]

    def test_skip_line(self, tmp_path):
        # Create a temporary test file
        test_file = tmp_path / "test_skip.txt"
        test_file.write_text("line1\nskip this line\nline3\n")

        with MockInputFileReader(test_file) as reader:
            lines = [line for line in reader]
        assert lines == ["line1", "line3"]
