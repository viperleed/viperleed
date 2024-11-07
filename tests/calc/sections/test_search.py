"""Tests for viperleed.calc section search."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from unittest.mock import MagicMock
from unittest.mock import patch

import pytest
from pytest_cases import parametrize

from viperleed.calc.sections.search import InconsistentV0ImagError
from viperleed.calc.sections.search import MaxIntensitiesError
from viperleed.calc.sections.search import NotEnoughSlotsError
from viperleed.calc.sections.search import ProcessKilledError
from viperleed.calc.sections.search import SearchError
from viperleed.calc.sections.search import SigbusError
from viperleed.calc.sections.search import _check_search_log

from ...helpers import not_raises


def patch_read(log_contents):
    """Make Path.read_text return `log_contents`."""
    return patch('pathlib.Path.read_text',
                 MagicMock(return_value=log_contents))


class TestSearchAg100:
    """Check the successful outcome of a structure optimization for Ag(100)."""

    def test_successful_run(self, search_files_ag100):
        """Check that structure search exits without errors."""
        assert not search_files_ag100.failed
        assert search_files_ag100.records is not None
        assert search_files_ag100.records.get_last_section_state('search')

    @parametrize('expected_file', ('search.steu',))
    def test_search_input_exist(self, search_files_ag100, expected_file):
        """Make sure input files were generated."""
        assert search_files_ag100.expected_file_exists(expected_file)

    @parametrize('expected_file', ('SD.TL', 'control.chem'))
    def test_search_raw_files_exist(self, search_files_ag100, expected_file):
        """Make sure the expected output files are present."""
        assert search_files_ag100.expected_file_exists(expected_file)

    @parametrize('expected_file', ('Search-report.pdf', 'Search-progress.pdf'))
    def test_search_pdf_files_exist(self, search_files_ag100, expected_file):
        """Make sure the expected PDF report files are present."""
        assert search_files_ag100.expected_file_exists(expected_file)


class TestCheckSearchLog:
    """Tests for the _check_search_log function."""

    @parametrize(path=(None, ''))
    def test_no_log(self, path):
        """Check no complaints if there is no log file."""
        with not_raises(SearchError):
            _check_search_log(path)

    def test_cannot_open_file(self, caplog):
        """Check logging when failing to open a log file."""
        msg = 'Could not read search log'
        with patch_read('Some log data') as mock_read:
            mock_read.side_effect = OSError
            _check_search_log('does_not_exist')
            assert msg in caplog.text

    _log_faulty = {
        'max intensity': (
            MaxIntensitiesError,
            'MAX. INTENS. IN THEOR. BEAM <kdsv1jb<v IS SUSPECT bcjab\n'
            '****** STOP PROGRAM ****** <sò48+jvblk< vh',
            ),
        'sigbus': (
            SigbusError,
            'received signal SIGBUS <sn6ò undefined portion of a memory object'
            ),
        'v0i': (
            InconsistentV0ImagError,
            'Average optical potential value in rf.info is incorrect:',
            ),
        'killed Intel': (
            ProcessKilledError,
            'APPLICATION TERMINATED WITH THE EXIT STRING: Killed (signal 9)',
            ),
        'killed GNU': (ProcessKilledError,
                       '=   KILLED BY SIGNAL: 9 (Killed)'),
        'slots': (NotEnoughSlotsError,
                  'There are not enough slots available in the system'),
        }

    @parametrize('exc,log_contents', _log_faulty.values(), ids=_log_faulty)
    def test_raises(self, exc, log_contents):
        """Check that `exc` is raised for given log-file contents."""
        with pytest.raises(exc), patch_read(log_contents):
            _check_search_log('some_log_file')

    _no_errors = {
        'truly no error': 'No errors',
        'max intensity, truncated_1': 'MAX. INTENS. IN THEOR. BEAM',
        'max intensity, truncated_2': 'BEAM <kdsv1jb<v IS SUSPECT',
        'max intensity, truncated_3': 'b\n****** STOP PROGRAM',
        'sigbus, truncated_1': 'received signal SIGBUS',
        'sigbus, truncated_2': 'SIGBUS <sn6ò undefined portion',
        'v0i, truncated': 'Average optical potential incorrect:',
        'killed Intel, truncated': 'EXIT STRING: Killed (signal 9)',
        'killed GNU, truncated': 'BY SIGNAL: 9',
        'slots, truncated': 'not enough slots',
        }

    @parametrize(log_contents=_no_errors.values(), ids=_no_errors)
    def test_no_errors(self, log_contents):
        """Check successful processing of an error-free log."""
        with not_raises(SearchError), patch_read(log_contents):
            _check_search_log('some_log_file')
