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

import shutil
import subprocess

from pytest_cases import fixture
from pytest_cases import parametrize_with_cases

from viperleed.calc import psgen
from viperleed.calc.files.phaseshifts import writePHASESHIFTS
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.time_utils import DateTimeFormat

from ...helpers import PSGEN_SEED
from ...helpers import private_random
from ..poscar_slabs import CasePOSCARSlabs as POSCARSlabs
from .conftest import SYSTEM_NAME

OUTPUT_FOLDER = 'stored_inputs/eeasisss_output'
RAW_OUTPUT_GLOB = 'PS.r.*'
REFERENCE = 'PHASESHIFTS'
ENERGIES = 50, 100, 10

# One case with a single site, and one with mixed-element sites.
CASES = (
    POSCARSlabs.case_poscar_ag100,
    POSCARSlabs.case_poscar_ptrh_100_1x3_o,
    )


@fixture(name='mock_eeasisss')
def fixture_mock_eeasisss(mocker):
    """Return a function making eeasisss 'produce' stored outputs."""
    def _mock(stored_output):
        def _run(args, **_):
            for raw_output in stored_output.glob(RAW_OUTPUT_GLOB):
                shutil.copy2(raw_output, raw_output.name)
            return subprocess.CompletedProcess(args, 0, 'stdout', '')
        return mocker.patch.object(psgen.subprocess, 'run',
                                   side_effect=_run)
    return _mock


@fixture(name='parsed_phaseshifts')
def fixture_parsed_phaseshifts(stub_tree, tmp_path, mock_eeasisss):
    """Return a function that writes PHASESHIFTS from stored outputs."""
    def _parse(slab, rpars, stored_output):
        rpars.paths.tensorleed = stub_tree
        rpars.systemName = SYSTEM_NAME
        rpars.timestamp = DateTimeFormat.FILE_SUFFIX.now()
        rpars.initTheoEnergies()
        rpars.THEO_ENERGIES = rpars.THEO_ENERGIES.from_value(ENERGIES)
        work = tmp_path / 'work'
        work.mkdir()
        mock_run = mock_eeasisss(stored_output)
        with execute_in_dir(work), private_random(psgen, PSGEN_SEED):
            firstline, phaseshifts = psgen.runPhaseshiftGen_old(slab, rpars)
        mock_run.assert_called_once()
        assert not any(work.glob(RAW_OUTPUT_GLOB))  # Read and removed
        written = work / REFERENCE
        writePHASESHIFTS(firstline, phaseshifts, written)
        return written
    return _parse


@parametrize_with_cases('args', cases=CASES)
def test_phaseshifts_from_stored_output(args, parsed_phaseshifts,
                                        check_reference, data_path,
                                        current_cases):
    """Check the PHASESHIFTS read from stored eeasisss outputs."""
    slab, rpars, *_ = args
    stored_output = (data_path / OUTPUT_FOLDER
                     / current_cases['args'].id)
    written = parsed_phaseshifts(slab, rpars, stored_output)
    check_reference(written, stored_output / REFERENCE)
