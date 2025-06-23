"""Tests for readBEAMLIST of viperleed.calc.files.beams."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-06-22'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.files.beams import readBEAMLIST


class TestReadBEAMLIST:
    """Tests for the readBEAMLIST function."""

    @staticmethod
    def read(file_path=None):
        """Call readBEAMLIST at `file_path`."""
        if file_path is None:
            return readBEAMLIST()
        return readBEAMLIST(file_path)

    def test_file_not_found(self, caplog):
        """Check complaints when the BEAMLIST file is not found."""
        with pytest.raises(FileNotFoundError):
            self.read()
        expect_log = 'Error opening BEAMLIST file.'
        assert expect_log in caplog.text

    def test_os_error(self, mocker, caplog):
        """Check complaints when the BEAMLIST file is not found."""
        mocker.patch('builtins.open', side_effect=OSError)
        with pytest.raises(OSError):
            self.read()
        expect_log = 'Error opening BEAMLIST file.'
        assert expect_log in caplog.text

    @parametrize(fname=('BEAMLIST', 'other_beamlist'))
    def test_success(self, fname, tmp_path, mocker, caplog):
        """Check the successful reading of the BEAMLIST file."""
        caplog.set_level(0)  # All messages
        beamlist_path = tmp_path / fname
        content = ['line 1\n', 'line 2\n', 'line 3\n']
        beamlist_path.write_text(''.join(content))
        assert self.read(beamlist_path) == content

        expect_log = 'BEAMLIST file was read successfully'
        assert expect_log in caplog.text
