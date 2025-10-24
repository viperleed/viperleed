"""Tests for module deltas of viperleed.calc.sections.

This module defines tests for running actual delta-amplitude
calculations on simple systems.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

import re

from viperleed.calc.constants import LOG_PREFIX


class TestDeltasAg100:
    """Test the successful outcome of a delta-amplitudes calculation."""

    def test_successful_run(self, delta_files_ag100):
        """Check that delta-amplitude calculation exits without errors."""
        assert not delta_files_ag100.failed
        assert delta_files_ag100.records is not None
        assert delta_files_ag100.records.get_last_state_for_section('delta')

    def test_delta_input_written(self, delta_files_ag100):
        """Check that an input file was correctly written."""
        assert delta_files_ag100.expected_file_exists('delta-input')

    def test_deltas_zip_created(self, delta_files_ag100):
        """Check that an archive with delta-amplitude files was created."""
        assert delta_files_ag100.expected_file_exists('Deltas/Deltas_001.zip')


class TestDeltasDomains:
    """Test execution of a multi-domain delta-amplitudes calculation."""

    def test_compile_folders_handled(self, delta_domains):
        """Ensure that compile folders/logs were processed correctly."""
        log_file = next(delta_domains.work_path.glob(f'{LOG_PREFIX}*.log'))
        # This regex covers
        #  (i) deletion of compilation folder
        # (ii) moving of compile_log file
        assert not re.search(r'WARNING.*Delta_Compile',
                             log_file.read_text())

    def test_successful_run(self, delta_domains):
        """Check that delta-amplitude calculation exits without errors."""
        assert not delta_domains.failed
        assert delta_domains.records is not None
        assert delta_domains.records.get_last_state_for_section('delta')
