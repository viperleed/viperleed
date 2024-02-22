"""Tests for section refcalc.

Created on 2023-07-28

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

import pytest


class TestRefCalc:
    """Test the successful outcome of a full-dynamical calculation."""

    def test_successful_run(self, refcalc_files):
        """Check that refcalc exits without errors."""
        assert not refcalc_files.failed

    @pytest.mark.parametrize('expected_file', ('THEOBEAMS.csv',))
    def test_refcalc_files_present(self, refcalc_files, expected_file):
        """Ensure the expected output files were generated."""
        assert refcalc_files.expected_file_exists(expected_file)
