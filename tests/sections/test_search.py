"""Tests for section search.

Created on 2023-07-28

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

import pytest


class TestSearchAg100:
    """Check the successful outcome of a structure optimization for Ag(100)."""

    def test_successful_run(self, search_files_ag100):
        """Check that structure search exits without errors."""
        assert not search_files_ag100.failed

    @pytest.mark.parametrize('expected_file', ('search.steu',))
    def test_search_input_exist(self, search_files_ag100, expected_file):
        """Make sure input files were generated."""
        assert search_files_ag100.expected_file_exists(expected_file)

    @pytest.mark.parametrize('expected_file', ('SD.TL', 'control.chem'))
    def test_search_raw_files_exist(self, search_files_ag100, expected_file):
        """Make sure the expected output files are present."""
        assert search_files_ag100.expected_file_exists(expected_file)

    @pytest.mark.parametrize(
        'expected_file',
        ('Search-report.pdf', 'Search-progress.pdf')
        )
    def test_search_pdf_files_exist(self, search_files_ag100, expected_file):
        """Make sure the expected PDF report files are present."""
        assert search_files_ag100.expected_file_exists(expected_file)
