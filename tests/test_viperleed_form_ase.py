"""Run tests for functionality in viperleed_from_ase.

Created on 2023-02-23

@author: Alex M. Imre
@author: Michele Riva

Define fixtures and test cases appropriate for the functionality
available in the viperleed_from_ase module of viperleed.
"""

from io import StringIO
from pathlib import Path
import sys

import ase.build
import numpy as np
import pytest

VPR_PATH = str(Path(__file__).resolve().parents[2])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Unfortunately no way to do this the correct way till we have
# an installable version of viperleed. The reason is the VPR_PATH
# bit above.
from viperleed import viperleed_from_ase as vpr_ase
from viperleed.tleedmlib.base import angle
from viperleed.tleedmlib.classes.slab import Slab
from viperleed.tleedmlib.files.beams import readOUTBEAMS
from viperleed.tleedmlib.files.poscar import readPOSCAR
# pylint: enable=wrong-import-position


INPUTS_ORIGIN = Path(__file__).parent / "fixtures"
INPUTS_ASE = INPUTS_ORIGIN / "from_ase"




# @Michele: put transformations to test here
_TRANSFORMATIONS_FOR_REFCALC = (
    # (orthogonal_transformation, isotropic_scaling)
    (None, None),
    )


@pytest.fixture(name="ase_Ni_100_1x1_cell")
def fixture_ase_nickel_cell():
    element = 'Ni'
    cell_1x1 = ase.build.fcc110(element, size=(1,1,6), vacuum=3)
    return cell_1x1

def test_ase_cell_correct(ase_Ni_100_1x1_cell):
    cell = ase_Ni_100_1x1_cell
    # make sure there are 6 atoms
    assert len(cell.positions) == 6

# @Michele: replace with Slab.from_ase method
def test_Ni_slab_from_ase(ase_Ni_100_1x1_cell):
    ase_cell = ase_Ni_100_1x1_cell
    slab = Slab(ase_cell)
    assert len(ase_cell.positions) == len(slab.atlist)


# @Michele: replace with new function - currently not working
def test_rot_mat_c(ase_Ni_100_1x1_cell):
    slab = Slab(ase_Ni_100_1x1_cell)
    a_before, b_before = slab.ucell[:,0], slab.ucell[:,1]
    theta = 30 # degrees
    rot_mat = vpr_ase.rot_mat_c(theta)
    slab.apply_matrix_transformation(rot_mat)
    a_after, b_after = slab.ucell[:,0], slab.ucell[:,1]
    get_angle = lambda x,y: np.arccos(np.dot(x,y)/np.linalg.norm(x)/np.linalg.norm(y))
    assert np.isclose(np.degrees(get_angle(a_before, a_after)), theta)
    assert np.isclose(np.degrees(get_angle(b_before, b_after)), theta)


@pytest.fixture(name="init_Ni_from_ase")
def fixture_run_from_ase_initialization(ase_Ni_100_1x1_cell, tmp_path_factory):
    ase_cell = ase_Ni_100_1x1_cell
    exec_path = tmp_path_factory.mktemp(basename='from_ase_Ni_100_init', numbered=True)
    inputs_path = INPUTS_ASE / "initialization"
    results = vpr_ase.run_from_ase(
        exec_path=exec_path,
        ase_object=ase_cell,
        inputs_path=inputs_path,
        cut_cell_c_fraction=0.0
    )
    return results, exec_path

@pytest.fixture(name="refcalc_Ni_from_ase", params=_TRANSFORMATIONS_FOR_REFCALC)
def fixture_run_from_ase_refcalc(ase_Ni_100_1x1_cell, tmp_path_factory, request):
    uc_transformation_matrix, uc_scaling = request.param
    ase_cell = ase_Ni_100_1x1_cell
    exec_path = tmp_path_factory.mktemp(basename='from_ase_Ni_100_init', numbered=True)
    inputs_path = INPUTS_ASE / "refcalc"
    results = vpr_ase.run_from_ase(
        exec_path=exec_path,
        ase_object=ase_cell,
        inputs_path=inputs_path,
        cut_cell_c_fraction=0.0
    )
    return results, exec_path

@pytest.fixture(name="refcalc_Ni_from_ase_beamlist")
def fixture_refcalc_thoeobeams(refcalc_Ni_from_ase):
    (theobeams, *_), _ = refcalc_Ni_from_ase
    theobeams_list = readOUTBEAMS(StringIO(theobeams))
    return theobeams_list

def test_returns_v0i(init_Ni_from_ase):
    (*_, v0i), _ = init_Ni_from_ase
    assert isinstance(v0i, float)

def test_init_writes_POSCAR(init_Ni_from_ase):
    _, exec_path = init_Ni_from_ase
    assert (exec_path / "POSCAR").is_file()

def test_init_writes_sensible_POSCAR(init_Ni_from_ase, ase_Ni_100_1x1_cell):
    _, exec_path = init_Ni_from_ase
    poscar_path = exec_path / "POSCAR"
    slab = readPOSCAR(poscar_path)
    assert len(slab.atlist) == len(ase_Ni_100_1x1_cell.positions)

def test_init_generates_VIBROCC(init_Ni_from_ase):
    _, exec_path = init_Ni_from_ase
    assert (exec_path / "work" / "VIBROCC").is_file()

@pytest.mark.parametrize('expected_file', (('BEAMLIST'), ('VIBROCC'), ('IVBEAMS')))
def test_run_from_ase_refcalc_files(refcalc_Ni_from_ase, expected_file):
    _, exec_path = refcalc_Ni_from_ase
    work_path = exec_path / "work"
    assert (work_path / expected_file).is_file()

def test_refcalc_output_not_empty(refcalc_Ni_from_ase):
    (theobeams, *_), _ = refcalc_Ni_from_ase
    assert theobeams != ""

def test_refcalc_output_correct_len(refcalc_Ni_from_ase_beamlist):
    theobeams_list = refcalc_Ni_from_ase_beamlist
    assert len(theobeams_list) == 4 # as defined in IVBEAMS

def test_refcalc_output_correct_energies(refcalc_Ni_from_ase_beamlist):
    theobeams_list = refcalc_Ni_from_ase_beamlist
    assert np.allclose(theobeams_list[0].energies, np.linspace(50, 70, 11))

def test_refcalc_output_non_nan(refcalc_Ni_from_ase_beamlist):
    theobeams_list = refcalc_Ni_from_ase_beamlist
    assert not all(np.isnan(theobeams_list[0].energies))


