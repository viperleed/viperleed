"""Tests for a full calc chain with TensErLEED 1.6.1.

Runs INITIALIZATION, REFCALC, DELTAS, and SEARCH one after the
other, for TensErLEED v1.6.1, which differs from the newer versions
in the format of several files. This covers the branches for older
versions in the modules that write the inputs of TensErLEED. No
stored TensErLEED data is used.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-28'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

# R factors of the chain, as in the calc log (4 decimals). They were
# identical in 17 runs, serial and with pytest-xdist. Mistakes in the
# inputs written for TensErLEED v1.6.1 that the other tests miss moved
# them by 0.0016 to 0.85.
R_FACTORS = {'refcalc': 0.1811, 'superpos': 0.1810}
R_FACTOR_TOLERANCE = 5e-4


class TestChainTL161:
    """Check the outcome of INITIALIZATION to SEARCH with TL 1.6.1."""

    @parametrize(section=('initialization', 'refcalc', 'delta', 'search'))
    def test_section_ran(self, chain_files_ag100_v161, section):
        """Check that each section ran, and that nothing failed."""
        assert chain_files_ag100_v161.section_ran(section)

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
    def test_file_exists(self, chain_files_ag100_v161, expected_file):
        """Ensure the expected input and output files were generated."""
        assert chain_files_ag100_v161.expected_file_exists(expected_file)

    def test_does_not_write_out_suffixed(self, chain_files_ag100_v161):
        """Check that none of the files generated has an _OUT suffix."""
        assert not chain_files_ag100_v161.out_suffixed_files()

    @parametrize(section=R_FACTORS)
    def test_final_r_factor(self, chain_files_ag100_v161, section):
        """Check the R factor of a section against the stored value."""
        r_factor = chain_files_ag100_v161.final_r_factor(section)
        assert r_factor is not None
        assert r_factor == pytest.approx(
            R_FACTORS[section],
            abs=R_FACTOR_TOLERANCE,
            )
