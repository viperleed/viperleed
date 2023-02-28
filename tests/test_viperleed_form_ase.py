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


_ASE_ATOMS = (
    "ase_ni_100_1x1_cell",
    "ase_ni_100_1x1_cell",
    )


@pytest.fixture(name="ase_Ni_100_1x1_cell")
def fixture_ase_nickel_cell():
    element = 'Ni'
    cell_1x1 = ase.build.fcc110(element, size=(1,1,6), vacuum=3)
    return cell_1x1


@pytest.mark.parametrize("ase_atoms, n_atoms", (("ase_Ni_100_1x1_cell", 6),))
def test_ase_n_atoms(ase_atoms, n_atoms):
    """Make sure `ase_atoms` has `n_atoms` atoms."""
    assert len(ase_atoms.positions) == n_atoms


# TODO: will need to replace with Slab.from_ase method
def slab_from_ase(ase_atoms):
    """Return a Slab from an ase.Atoms object."""
    return Slab(ase_atoms)


@pytest.mark.parametrize("ase_atoms", _ASE_ATOMS)
def test_n_atoms_from_ase(ase_atoms):
    slab = slab_from_ase(ase_atoms)
    assert len(ase_atoms.positions) == len(slab.atlist)


THETA = 14.7  # degrees

@pytest.mark.parametrize("ase_atoms", _ASE_ATOMS)
def test_rot_mat_c(ase_atoms):
    slab = slab_from_ase(ase_atoms)
    rot_mat = vpr_ase.rot_mat_c(THETA)
    ucell_before = slab.ucell.T.copy()
    slab.apply_matrix_transformation(rot_mat)
    ucell_after = slab.ucell.T

    # a and b unit vectors
    for _before, _after in zip(ucell_before[:2], ucell_after[:2]):
        assert np.isclose(_before[2], _after[2])  # Same z
        assert np.isclose(np.degrees(angle(_before, _after)), THETA)
    # c unit vector
    assert np.allclose(ucell_before[2], ucell_after[2])


@pytest.fixture(name="init_Ni_from_ase")
def fixture_run_from_ase_initialization(ase_Ni_100_1x1_cell, tmp_path_factory):
    ase_atoms = ase_Ni_100_1x1_cell
    exec_path = tmp_path_factory.mktemp(basename='from_ase_Ni_100_init',
                                        numbered=True)
    inputs_path = INPUTS_ASE / "initialization"
    results = vpr_ase.run_from_ase(
        exec_path=exec_path,
        ase_object=ase_atoms,
        inputs_path=inputs_path,
        cut_cell_c_fraction=0.0
    )
    return results, exec_path

@pytest.fixture(name="refcalc_Ni_from_ase", params=_TRANSFORMATIONS_FOR_REFCALC)
def fixture_run_from_ase_refcalc(ase_Ni_100_1x1_cell, tmp_path_factory, request):
    uc_transformation_matrix, uc_scaling = request.param
    ase_atoms = ase_Ni_100_1x1_cell
    exec_path = tmp_path_factory.mktemp(basename='from_ase_Ni_100_init',
                                        numbered=True)
    inputs_path = INPUTS_ASE / "refcalc"
    results = vpr_ase.run_from_ase(
        exec_path=exec_path,
        ase_object=ase_atoms,
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


