"""Tests for the parsing of the raw output of eeasisss in psgen.

The eeasisss executable is replaced by a mock that 'produces' the raw
output files (PS.r.*) stored in the test data. The phase shifts that
viperleed.calc.psgen reads from them are compared to stored references.
No TensErLEED tree is needed. The stored outputs were produced by
eeasisss for the stored eeasisss inputs, but with a narrower energy
range (see ENERGIES) to keep them small.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from pytest_cases import fixture
from pytest_cases import parametrize_with_cases

from viperleed.calc import psgen

from ..poscar_slabs import CasePOSCARSlabs as POSCARSlabs
from .stored_eeasisss import RAW_OUTPUT_GLOB
from .stored_eeasisss import REFERENCE
from .stored_eeasisss import STORED_PATH
from .stored_eeasisss import produce_stored_outputs
from .stored_eeasisss import run_psgen

ENERGIES = 50, 100, 10  # Of these stored outputs, not of stored_eeasisss

# One case with a single site, and one with mixed-element sites.
CASES = (
    POSCARSlabs.case_poscar_ag100,
    POSCARSlabs.case_poscar_ptrh_100_1x3_o,
    )


@fixture(name='parsed_phaseshifts')
def fixture_parsed_phaseshifts(stub_tree, tmp_path, mocker):
    """Return a function that writes PHASESHIFTS from stored outputs."""
    def _parse(slab, rpars, stored_output):
        mock_run = mocker.patch.object(
            psgen.subprocess,
            'run',
            side_effect=produce_stored_outputs(stored_output),
            )
        work, *_ = run_psgen(slab, rpars, stub_tree, tmp_path,
                             energies=ENERGIES)
        mock_run.assert_called_once()
        assert not any(work.glob(RAW_OUTPUT_GLOB))  # Read and removed
        return work / REFERENCE
    return _parse


@parametrize_with_cases('args', cases=CASES)
def test_phaseshifts_from_stored_output(args, parsed_phaseshifts,
                                        check_reference, current_cases):
    """Check the PHASESHIFTS read from stored eeasisss outputs."""
    slab, rpars, *_ = args
    stored_output = STORED_PATH / current_cases['args'].id
    written = parsed_phaseshifts(slab, rpars, stored_output)
    check_reference(written, stored_output / REFERENCE)
