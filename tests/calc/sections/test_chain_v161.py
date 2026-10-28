"""Tests for a full calc chain with TensErLEED 1.6.1.

Runs INITIALIZATION, REFCALC, DELTAS, and SEARCH one after the
other, for TensErLEED v1.6.1 (see D-013), which differs from the
newer versions in the format of several files. This covers the
branches for older versions in the modules that write the inputs
of TensErLEED. No stored TensErLEED data is used.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-28'
__license__ = 'GPLv3+'

from pytest_cases import parametrize


class TestChainTL161:
    """Check the outcome of INITIALIZATION to SEARCH with TL 1.6.1."""

    @parametrize(section=('initialization', 'refcalc', 'delta', 'search'))
    def test_section_ran(self, chain_files_ag100_v161, section):
        """Check that each section ran, and that nothing failed."""
        assert not chain_files_ag100_v161.failed
        records = chain_files_ag100_v161.records
        assert records is not None
        assert records.get_last_state_for_section(section)

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
        out_suffixed = chain_files_ag100_v161.work_path.rglob('*_OUT*')
        # Skip R_OUT files
        # pylint: disable-next=magic-value-comparison
        out_suffixed = (f for f in out_suffixed if 'R=' not in f.name)
        assert not any(out_suffixed)
