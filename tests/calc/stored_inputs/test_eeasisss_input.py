"""Tests for the eeasisss-input file written by viperleed.calc.psgen.

The input of eeasisss is compared to stored references. No TensErLEED
tree is needed: the executable and the data files are replaced by
stubs, and eeasisss itself is never run.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize_with_cases

from viperleed.calc import psgen
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.periodic_table import PERIODIC_TABLE
from viperleed.calc.lib.time_utils import DateTimeFormat

from ...helpers import CustomTestException
from ...helpers import PSGEN_SEED
from ...helpers import private_random
from ..poscar_slabs import CasePOSCARSlabs as POSCARSlabs
from .conftest import SYSTEM_NAME
from .conftest import pin_clock

INPUT_NAME = 'eeasisss-input'

# Cases of POSCARSlabs that carry parameters (presets or mixed sites)
CASES = (
    POSCARSlabs.case_poscar_ag100,
    POSCARSlabs.case_poscar_diamond,
    POSCARSlabs.case_poscar_fe3o4_001_cod,
    POSCARSlabs.case_poscar_lsmo_001_rt2,
    POSCARSlabs.case_poscar_ptrh_100_1x3_o,
    POSCARSlabs.case_poscar_ru0001_rt3_te,
    )


@fixture(name='stub_tree')
def fixture_stub_tree(tmp_path):
    """Return a fake TensErLEED directory with stub eeasisss/data."""
    tree = tmp_path / 'TensErLEED'
    tree.mkdir()
    (tree / 'eeasisss').touch()
    (tree / 'seSernelius').touch()
    for element in PERIODIC_TABLE:
        element_dir = tree / 'atom_density_files' / element
        element_dir.mkdir(parents=True)
        (element_dir / f'chgden{element}').touch()
    return tree


@fixture(name='written_input')
def fixture_written_input(stub_tree, tmp_path, mocker, monkeypatch):
    """Return a function that writes the eeasisss input for a slab."""
    pin_clock(monkeypatch)

    def _write(slab, rpars):
        rpars.paths.tensorleed = stub_tree
        rpars.systemName = SYSTEM_NAME
        rpars.timestamp = DateTimeFormat.FILE_SUFFIX.now()
        rpars.initTheoEnergies()
        work = tmp_path / 'work'
        work.mkdir()
        # Stop right after the input is written and 'run'
        mocker.patch.object(psgen.subprocess, 'run',
                            side_effect=CustomTestException)
        with execute_in_dir(work), private_random(psgen, PSGEN_SEED), \
                pytest.raises(CustomTestException):
            psgen.runPhaseshiftGen_old(slab, rpars)
        return work / INPUT_NAME
    return _write


@parametrize_with_cases('args', cases=CASES)
def test_eeasisss_input(args, written_input, check_reference, data_path,
                        current_cases):
    """Check the eeasisss-input file written for a slab."""
    slab, rpars, *_ = args
    case_id = current_cases['args'].id
    reference = data_path / 'stored_inputs/eeasisss' / case_id / INPUT_NAME
    check_reference(written_input(slab, rpars), reference)
