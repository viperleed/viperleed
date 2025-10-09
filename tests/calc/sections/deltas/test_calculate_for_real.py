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

    def test_successful_run(self, delta_domains):
        """Check that delta-amplitude calculation exits without errors."""
        assert not delta_domains.failed
        assert delta_domains.records is not None
        assert delta_domains.records.get_last_state_for_section('delta')
