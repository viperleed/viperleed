"""Tests for a full calc chain with TensErLEED 1.6.1 and the newest.

Runs INITIALIZATION, REFCALC, DELTAS, and SEARCH one after the
other, for TensErLEED v1.6.1, which differs from the newer versions
in the format of several files, and for the most recent version. This
covers the branches for older versions, as well as those for newer
ones, in the modules that write the inputs of TensErLEED. No stored
TensErLEED data is used.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-28'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from .conftest import CHAIN_R_FACTORS

R_FACTOR_TOLERANCE = 5e-4


class TestChain:
    """Check the outcome of INITIALIZATION to SEARCH."""

    @parametrize(section=('initialization', 'refcalc', 'delta', 'search'))
    def test_section_ran(self, chain_files_ag100, section):
        """Check that each section ran, and that nothing failed."""
        assert chain_files_ag100.section_ran(section)

    @parametrize(
        expected_file=(
            'refcalc-FIN',       # Refcalc input
            'THEOBEAMS.csv',     # Refcalc output
            'delta-input',       # Delta input
            'Deltas/Deltas_001.zip',
            'search.steu',       # Search input
            'SD.TL',             # Search output
            'control.chem',
            'Search-report.pdf',
            'Search-progress.pdf',
            )
        )
    def test_file_exists(self, chain_files_ag100, expected_file):
        """Ensure the expected input and output files were generated."""
        assert chain_files_ag100.expected_file_exists(expected_file)

    def test_does_not_write_out_suffixed(self, chain_files_ag100):
        """Check that none of the files generated has an _OUT suffix."""
        assert not chain_files_ag100.out_suffixed_files()

    @parametrize(section=('refcalc', 'superpos'))
    def test_final_r_factor(self, chain_files_ag100, chain_tl_version,
                            section):
        """Check the R factor of a section against the stored value."""
        r_factor = chain_files_ag100.final_r_factor(section)
        assert r_factor is not None
        assert r_factor == pytest.approx(
            CHAIN_R_FACTORS[str(chain_tl_version)][section],
            abs=R_FACTOR_TOLERANCE,
            )
