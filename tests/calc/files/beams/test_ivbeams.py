"""Tests for readIVBEAMS of viperleed.calc.files.beams."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-06-22'
__license__ = 'GPLv3+'

import logging
import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.files.beams import readIVBEAMS

_MODULE = 'viperleed.calc.files.beams'


class MockBeam:
    """A replacement for the Beam class."""

    def __init__(self, *indices):
        """Initialize beam from indices."""
        if len(indices) == 1:
            indices = indices[0]
        self.indices = indices

    def __eq__(self, other):
        """Return whether this beam is the same as another one."""
        try:
            return self.indices == other.indices
        except AttributeError:
            return NotImplemented


class TestReadIVBEAMS:
    """Tests for the readIVBEAMS function."""

    @fixture(name='make_ivbeams')
    def factory_make_ivbeams(self, tmp_path):
        """Create an IVBEAMS file in a temporary directory."""
        def _write(contents, fname='IVBEAMS'):
            ivbeams_file = tmp_path/fname
            ivbeams_file.write_text(contents, encoding='utf-8')
            return ivbeams_file
        return _write

    @fixture(name='mock_beam_cls', autouse=True)
    def fixture_mock_beam_cls(self, mocker):
        """Replace the Beam class with a mock."""
        mocker.patch(f'{_MODULE}.Beam', MockBeam)

    _empty = {
        'truly empty': '',
        'comments only': '''
!Commented header
#(1 | 0) commented beam
%(-5/2 2) another commented beam
''',
        'header only': 'header line',
        }

    @parametrize(contents=_empty.values(), ids=_empty)
    def test_empty_file(self, contents, make_ivbeams, caplog):
        """Check that no beams are read from an empty file."""
        caplog.set_level(logging.INFO)  # Except debug
        result = readIVBEAMS(make_ivbeams(contents))
        assert not result
        assert not caplog.text

    def test_file_not_found(self):
        """Check complaints if IVBEAMS is not found."""
        with pytest.raises(FileNotFoundError):
            readIVBEAMS('nonexistent_file')

    def test_first_line_has_data(self, make_ivbeams, caplog):
        """Check skipping of data in the header line."""
        result = readIVBEAMS(make_ivbeams('1 2\n3 4\n'))
        expect = [MockBeam(3, 4)]
        assert result == expect
        assert re.match(r'.*first line.* will not be read', caplog.text)

    _invalid = {
        'float': 'not_a_float 3',
        'fraction': '3/5 4/a',
        }

    @parametrize(invalid_line=_invalid.values(), ids=_invalid)
    def test_invalid_line(self, invalid_line, make_ivbeams, caplog):
        """Check exceptions when a line containing a non-float is found."""
        file = make_ivbeams(f'HEADER\n{invalid_line}')
        with pytest.raises(ValueError):
            readIVBEAMS(file)
        expect_log = f'Error reading IVBEAMS line: {invalid_line}'
        assert expect_log in caplog.text

    @parametrize(one_item=('only_one_item', '12345'))
    def test_one_element_line(self, one_item, make_ivbeams, caplog):
        """Check warnings when one line contains too little data."""
        file = make_ivbeams(
            'HEADER\n'
            f'{one_item}\n'  # Skipped
            '3 4'
            )
        result = readIVBEAMS(file)
        expect = [MockBeam(3, 4)]
        assert result == expect
        assert re.match('.*line.*only one element.*skipped', caplog.text)

    def test_open_fails(self, mocker):
        """Check complaints if IVBEAMS is not found."""
        mocker.patch('builtins.open', side_effect=Exception)
        with pytest.raises(Exception):
            readIVBEAMS()

    _header_line = (
        'HEADER LINE',
        'one/two HEADER LINE WITH A SLASH',
        )

    @parametrize(header=_header_line)
    def test_valid_data(self, header, make_ivbeams, caplog):
        """Check the correct reading of data from IVBEAMS."""
        caplog.set_level(0)  # All messages
        ivbeams_path = make_ivbeams(
            f'{header}\n'
            '1 2\n'
            '(3) | 4 more data here\n'
            '1/2 3/4\n'
            '1.0 2.000\n'  # Duplicate — should be skipped
            '   \n'        # Empty - should be skipped
            )
        result = readIVBEAMS(ivbeams_path)
        expect = [MockBeam(1.0, 2.0),
                  MockBeam(3.0, 4.0),
                  MockBeam(0.5, 0.75)]
        assert result == expect
        expect_log = 'IVBEAMS file was read successfully'
        assert expect_log in caplog.text
