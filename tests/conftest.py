"""Module conftest of viperleed.tests.

Created on 2023-02-28

@author: Michele Riva
@author: Alexander M. Imre

Contains some useful general definitions that can be used when creating
or running tests.
"""


# Think about a decorator for injecting fixtures.
# Some ideas at
# https://github.com/pytest-dev/pytest/issues/2424
# https://github.com/pytest-dev/pytest/issues/6322
# https://github.com/nteract/testbook/issues/4

from pathlib import Path
import os
import sys
import shutil
from zipfile import ZipFile

import pytest
import numpy as np

vpr_path = str(Path(__file__).parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

from viperleed.tleedm import run_tleedm
from viperleed.tleedmlib import symmetry
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.classes.slab import Slab
from viperleed.tleedmlib.files import parameters, poscar
from viperleed.tleedmlib.files.vibrocc import readVIBROCC
from viperleed.tleedmlib.files.displacements import readDISPLACEMENTS, readDISPLACEMENTS_block

from .helpers import TEST_DATA, POSCAR_PATH


_EXAMPLE_POSCAR_EXPECTATIONS = [("POSCAR_Ag(100)", 6, 'p4m', 0),
                                ("POSCAR_STO(110)-4x1", 136, 'pm', 0),
                                ("POSCAR_TiO2", 540, 'pmm', -1),
                                ("POSCAR_diamond", 96, 'pm', 89),
                                ("POSCAR_36C_p6m", 36, 'p6m', 0),
                                ("POSCAR_36C_cm", 36,'cm', 0),]
                               #("POSCAR_Fe3O4_SCV", 83, 'cmm', 50)]            #TODO: Phaseshift generation fails. Why? @Fkraushofer (worked in fkpCurie:Florian_OldLocalTests/Fe3O4-001-SCV/history/t000.r013_211220-133452)


_EXAMPLE_POSCARs = [file.name for file in POSCAR_PATH.glob('POSCAR*')]

TENSORLEED_PATH = Path(vpr_path) / "viperleed" / "tensorleed"


@pytest.fixture(scope='session')
def poscars_path():
    """Return the Path to the directory containing POSCAR files."""
    return POSCAR_PATH


@pytest.fixture(scope='session', name='data_path')
def fixture_data_path():
    """Return the Path to the top-level folder containing test data."""
    return TEST_DATA


@pytest.fixture(scope='session', name='tensorleed_path')
def fixture_tensorleed_path():
    """Return the Path to the top-level tree with tensor-LEED source code."""
    return TENSORLEED_PATH

@pytest.fixture()
def ag100_parameters_example():
    # read Ag(100) POSCAR and PARAMETERS files
    slab = poscar.read(TEST_DATA / 'Ag(100)' / 'initialization' / 'POSCAR')
    rpars = parameters.readPARAMETERS(TEST_DATA / 'Ag(100)' / 'initialization' / 'PARAMETERS')
    # interpret PARAMETERS file
    interpreter = parameters.ParameterInterpreter(rpars)
    interpreter.interpret(slab)
    symmetry.findSymmetry(slab, rpars)
    symmetry.findBulkSymmetry(slab, rpars)
    return (rpars, slab)


@pytest.fixture(scope="function", params=_EXAMPLE_POSCARs)
def example_poscars(request):
    file_path = POSCAR_PATH / request.param
    slab = poscar.read(file_path)
    return slab

@pytest.fixture(scope="function", params=_EXAMPLE_POSCAR_EXPECTATIONS)
def slab_and_expectations(request):
    filename, expected_n_atoms, expected_pg, offset_at = request.param
    file_path = POSCAR_PATH / filename
    pos_slab = poscar.read(file_path)
    return (pos_slab, expected_n_atoms, expected_pg, offset_at)

@pytest.fixture(scope="function")
def slab_pg_rp(slab_and_expectations):
    slab, *_ = slab_and_expectations
    rp = Rparams()
    slab.fullUpdate(rp)
    pg = symmetry.findSymmetry(slab, rp, output=False)
    symmetry.enforceSymmetry(slab, rp)
    return slab, pg, rp


@pytest.fixture(scope='function')
def manual_slab_3_atoms():
    slab = Slab()
    slab.ucell = np.diag([3., 4., 5.])
    positions = (np.array([-0.25, 0, 0]),
                 np.array([0.00, 0, 0]),
                 np.array([0.25, 0, 0]))
    slab.atlist = [Atom('C', pos, i+1, slab)
                   for i, pos in enumerate(positions)]
    param = Rparams()
    slab.fullUpdate(param)
    return slab


@pytest.fixture()
def manual_slab_1_atom_trigonal():
    slab = Slab()
    slab.ucell = np.array([[ 1, 0, 0],
                           [-2.3, 3, 0],
                           [ 1, 2, 3]],dtype=float).T
    slab.atlist = [Atom('C', np.array([0.2, 0.7, 0.1]), 1, slab),]  # "random" position
    param = Rparams()
    slab.fullUpdate(param)
    return slab


@pytest.fixture()
def ag100_slab_param(poscars_path):
    slab = poscar.read(poscars_path /"POSCAR_Ag(100)")
    param = Rparams()
    param.N_BULK_LAYERS = 1
    slab.fullUpdate(param)
    return slab, param


@pytest.fixture()
def ag100_slab_with_displacements_and_offsets(ag100_slab_param, data_path):
    slab, param = ag100_slab_param
    vibrocc_path = data_path / "Ag(100)" / "mergeDisp" / "VIBROCC"
    displacements_path = data_path / "Ag(100)" / "mergeDisp" / "DISPLACEMENTS_mixed"
    readVIBROCC(param, slab, str(vibrocc_path))
    readDISPLACEMENTS(param, str(displacements_path))
    readDISPLACEMENTS_block(param, slab, param.disp_blocks[param.search_index])
    return slab, param
