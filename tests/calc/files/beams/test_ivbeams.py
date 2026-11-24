"""Tests for readIVBEAMS of viperleed.calc.files.beams."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2025-06-22'
__license__ = 'GPLv3+'

import itertools
import logging
import re
from types import SimpleNamespace

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize
from quicktions import Fraction

from viperleed.calc.classes.beam import Beam
from viperleed.calc.files.beams import readIVBEAMS
from viperleed.calc.files.beams import writeIVBEAMS

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


class TestReadIVBEAMSRoundTrip:
    """Tests for reading back IVBEAMS files written by writeIVBEAMS."""

    _fractions = (
        Fraction(0),
        Fraction(1, 2),
        Fraction(1, 3),
        Fraction(2, 3),
        Fraction(1, 4),
        Fraction(1, 6),
        Fraction(1, 7),
        Fraction(-1, 3),
        Fraction(3, 2),
        )

    @fixture(name='written_beams')
    def fixture_written_beams(self, tmp_path, mocker):
        """Return an IVBEAMS file from writeIVBEAMS and its beams."""
        mocker.patch(f'{_MODULE}.leedbase.getLEEDdict', return_value={})
        hkfracs = itertools.product(self._fractions, repeat=2)
        rpars = SimpleNamespace(
            AVERAGE_BEAMS=False,
            expbeams=[Beam(hkfrac) for hkfrac in hkfracs],
            )
        ivbeams_file = tmp_path/'IVBEAMS'
        written = writeIVBEAMS(None, rpars, filename=ivbeams_file)
        assert len(written) == len(self._fractions)**2
        return ivbeams_file, written

    def test_same_beams(self, written_beams):
        """Check that beams read back with their fractional indices."""
        ivbeams_file, written = written_beams
        read_back = readIVBEAMS(ivbeams_file)
        assert [b.hkfrac for b in read_back] == [b.hkfrac for b in written]
        assert [b.label for b in read_back] == [b.label for b in written]

    def test_same_float_indices(self, written_beams):
        """Check that float indices match the decimals in the file."""
        ivbeams_file, written = written_beams
        read_back = readIVBEAMS(ivbeams_file)
        assert len(read_back) == len(written)
        for beam, written_beam in zip(read_back, written):
            from_file = tuple(float(f'{hk:.6f}') for hk in written_beam.hk)
            assert beam.hk == pytest.approx(from_file, abs=1e-6)

    def test_duplicate_decimal_and_fraction(self, tmp_path):
        """Check that a decimal and a fraction of one beam merge."""
        ivbeams_file = tmp_path/'IVBEAMS'
        ivbeams_file.write_text('HEADER\n'
                                '  0.333333  0.000000\n'
                                '1/3 0\n',
                                encoding='utf-8')
        read_back = readIVBEAMS(ivbeams_file)
        assert len(read_back) == 1
        assert read_back[0].hkfrac == (Fraction(1, 3), Fraction(0))

    @parametrize(indices=(('1/127', '0'), ('1/128', '0'), ('0.005', '0')))
    def test_large_denominator(self, indices, tmp_path):
        """Check indices far from fractions with small denominators."""
        ivbeams_file = tmp_path/'IVBEAMS'
        ivbeams_file.write_text(f'HEADER\n{indices[0]} {indices[1]}\n',
                                encoding='utf-8')
        read_back = readIVBEAMS(ivbeams_file)
        # Beam from the float values of the indices, as before
        expect = Beam(tuple(float(Fraction(index)) for index in indices))
        assert len(read_back) == 1
        beam = read_back[0]
        assert beam.hk == pytest.approx(expect.hk, abs=1e-6)
        assert beam.hkfrac == expect.hkfrac
        assert beam.label == expect.label

    def test_distinct_large_denominators(self, tmp_path):
        """Check that beams with denominators above 99 do not merge."""
        ivbeams_file = tmp_path/'IVBEAMS'
        ivbeams_file.write_text('HEADER\n'
                                '1/127 0\n'
                                '1/128 0\n',
                                encoding='utf-8')
        read_back = readIVBEAMS(ivbeams_file)
        assert len(read_back) == 2
