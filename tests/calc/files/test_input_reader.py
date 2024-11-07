"""Tests for module viperleed.calc.files.iorfactor."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-11-06'
__license__ = 'GPLv3+'

from io import StringIO

import pytest
from pytest_cases import parametrize

from viperleed.calc.files.input_reader import InputStreamReader
from viperleed.calc.files.input_reader import InputFileReader
from viperleed.calc.files.input_reader import ShouldSkipLineError

from ...helpers import not_raises


def raise_on_skip(self, line):  # pylint: disable=unused-argument
    """Raise if line contains 'skip'."""
    # pylint: disable-next=magic-value-comparison
    if 'skip' in line:
        raise ShouldSkipLineError('Line contains \'skip\'')
    return line.strip()


# pylint: disable-next=too-few-public-methods  # Only one abstract
class MockInputFileReader(InputFileReader):
    """A concrete InputFileReader for testing."""

    _read_one_line = raise_on_skip


# pylint: disable-next=too-few-public-methods  # Only one abstract
class MockInputStreamReader(InputStreamReader):
    """A concrete InputStreamReader for testing."""

    _read_one_line = raise_on_skip


class TestInputStreamReader:
    """Tests for the InputStreamReader class."""

    _valid = {
        'no skip': ('line1\nline2\nline3\n', ('line1', 'line2', 'line3')),
        'skip one': ('line1\nskip this line\nline3\n', ('line1', 'line3')),
        }

    @parametrize('lines,expect', _valid.values(), ids=_valid)
    def test_iteration(self, lines, expect):
        """Check that iteration over a stream yields the expected lines."""
        reader = MockInputStreamReader(StringIO(lines))
        assert tuple(reader) == expect

    def test_invalid_input(self):
        """Check complaints when a non-stream object is given."""
        with pytest.raises(TypeError):
            MockInputStreamReader('not a stream')


class TestInputFileReader:
    """Tests for the InputFileReader class."""

    def test_file_not_found(self):
        """Check complaints when accessing a non-existing file."""
        non_existent_file = '_this_file__DOES_non_exist.funny_ext'
        reader = MockInputFileReader(non_existent_file)
        with pytest.raises(FileNotFoundError), reader:
            pass

    _files = {
        'no skip': ('line1\nline2\nline3\n', ('line1', 'line2', 'line3')),
        'skip one': ('line1\nskip this line\nline3\n', ('line1', 'line3')),
        }

    @parametrize('contents,expect', _files.values(), ids=_files)
    def test_context_manager(self, contents, expect, tmp_path):
        """Check correct reading from a file."""
        test_file = tmp_path / 'test.txt'
        test_file.write_text(contents)

        with MockInputFileReader(test_file) as reader:
            lines = tuple(reader)
        assert lines == expect

    def test_exit_before_enter(self, tmp_path):
        """Check that exiting before entering a context does not complain."""
        file = tmp_path / 'test.txt'
        file.touch()
        reader = MockInputFileReader(file)
        with not_raises(AttributeError):
            reader.__exit__(None, None, None)

    def test_current_line_updated(self, tmp_path):
        """Check that the line numbers a tracked."""
        file = tmp_path / 'test.txt'
        file.write_text('line1\nline2\nline3\n')
        with MockInputFileReader(file) as reader:
            # pylint: disable=protected-access            # OK in tests
            assert not reader._current_line
            next(reader)
            assert reader._current_line
        with MockInputFileReader(file) as reader:
            # pylint: disable-next=protected-access       # OK in tests
            assert not reader._current_line  # Reset after __enter__
