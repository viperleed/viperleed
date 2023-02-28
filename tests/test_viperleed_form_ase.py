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


@pytest.fixture(name="ase_ni_100_1x1_cell")
def fixture_ase_nickel_cell():
    """Return an ase.Atoms Ni(100)-1x1 with 6 layers."""
    element = 'Ni'
    cell_1x1 = ase.build.fcc110(element, size=(1,1,6), vacuum=3)
    return cell_1x1


@pytest.mark.parametrize("ase_atoms, n_atoms", (("ase_ni_100_1x1_cell", 6),))
def test_ase_n_atoms(ase_atoms, n_atoms):
    """Make sure `ase_atoms` has `n_atoms` atoms."""
    assert len(ase_atoms.positions) == n_atoms


# TODO: will need to replace with Slab.from_ase method
def slab_from_ase(ase_atoms):
    """Return a Slab from an ase.Atoms object."""
    return Slab(ase_atoms)


@pytest.mark.parametrize("ase_atoms", _ASE_ATOMS)
def test_n_atoms_from_ase(ase_atoms):
    """Make sure the number of atoms in Slab match those in ase.Atoms."""
    slab = slab_from_ase(ase_atoms)
    assert len(ase_atoms.positions) == len(slab.atlist)


def test_rotation_matrices():
    """Test that rotation matrices are as expected."""
    raise NotImplementedError


THETA = 14.7  # degrees

@pytest.mark.parametrize("ase_atoms", _ASE_ATOMS)
def test_rot_mat_c(ase_atoms):
    """Assert correct rotation around z axis."""
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


# TODO: find a nice way to generate multiple fixtures dynamically.
# See also notes in helpers.py. Difficulties:
# - make temp path dependent on the name of the first argument in
#   the signature
# - pick the right folder in "initialization"
@pytest.fixture(name="run_from_ase_initialization")
def fixture_run_from_ase_initialization(ase_ni_100_1x1_cell, tmp_path_factory):
    """Return the results of an initialization run.

    Parameters
    ----------
    ase_ni_100_1x1_cell : pytest.fixture
        The ase.Atoms object from which to run.
    tmp_path_factory : pytest.fixture
        The pytest default name for the temporary directory factory.

    Returns
    -------
    results : tuple
        The return value of the call to run_from_ase.
    exec_path : Path
        The pat to the temporary directory created where the
        calculation was run.
    ase_atoms : ase.Atoms
        The Atoms object fed to run_from_ase.
    """
    ase_atoms = ase_ni_100_1x1_cell
    exec_path = tmp_path_factory.mktemp(basename='from_ase_Ni_100_init',
                                        numbered=True)
    inputs_path = INPUTS_ASE / "initialization"
    results = vpr_ase.run_from_ase(
        exec_path=exec_path,
        ase_object=ase_atoms,
        inputs_path=inputs_path,
        cut_cell_c_fraction=0.0
        )
    return results, exec_path, ase_atoms


# TODO: find a nice way to generate multiple fixtures dynamically.
@pytest.fixture(name="run_from_ase_refcalc",
                params=_TRANSFORMATIONS_FOR_REFCALC)
def fixture_run_from_ase_refcalc(ase_ni_100_1x1_cell,
                                 tmp_path_factory, request):
    """Return the results of a reference calculation run.

    Parameters
    ----------
    ase_ni_100_1x1_cell : pytest.fixture
        The ase.Atoms object from which to run.
    tmp_path_factory : pytest.fixture
        The pytest default name for the temporary directory factory.
    request : pytest.fixture
        The pytest default `equest` fixture for accessing
        fixture-level parameters. Used to access the transform
        to be applied to the slab.

    Returns
    -------
    results : tuple
        The return value of the call to run_from_ase.
    exec_path : Path
        The pat to the temporary directory created where the
        calculation was run.
    ase_atoms : ase.Atoms
        The Atoms object fed to run_from_ase.
    """
    ase_atoms = ase_ni_100_1x1_cell
    exec_path = tmp_path_factory.mktemp(basename='from_ase_Ni_100_init',
                                        numbered=True)
    inputs_path = INPUTS_ASE / "refcalc"
    results = vpr_ase.run_from_ase(
        exec_path=exec_path,
        ase_object=ase_atoms,
        inputs_path=inputs_path,
        cut_cell_c_fraction=0.0
        )
    return results, exec_path, ase_atoms


@pytest.fixture(name="refcalc_thoeobeams")
def fixture_refcalc_thoeobeams(run_from_ase_refcalc):
    """Return a list of full-dynamically calculated beams."""
    (theobeams, *_), *_ = run_from_ase_refcalc
    return readOUTBEAMS(StringIO(theobeams))


def test_init_output_empty(run_from_ase_initialization):
    """Ensure that run_from_ase initialization returns an empty output."""
    (*all_except_v0i, _), *_ = run_from_ase_initialization
    assert all(not s for s in all_except_v0i)


def test_returns_v0i(run_from_ase_initialization):
    """Make sure the last return of run_from_ase is V0i."""
    (*_, v0i), *_ = run_from_ase_initialization
    assert isinstance(v0i, float)


@pytest.mark.parametrize('file', ('POSCAR', 'VIBROCC'))
def test_init_writes_file(run_from_ase_initialization, file):
    """Ensure that run_from_ase writes `file` during initialization."""
    _, exec_path, _ = run_from_ase_initialization
    assert (exec_path / file).is_file()


@pytest.mark.parametrize('file', ('POSCAR', 'VIBROCC'))
def test_init_writes_workfile(run_from_ase_initialization, file):
    """Ensure that run_from_ase writes work/`file` during initialization."""
    _, exec_path, _ = run_from_ase_initialization
    assert (exec_path / "work"/ file).is_file()


# TODO: perhaps it would be even better to store somewhere
# checksums for files to be generated, and check them over here
def test_init_writes_sensible_poscar(run_from_ase_initialization):
    """Ensure that written POSCAR has the right number of atoms."""
    _, exec_path, ase_atoms = run_from_ase_initialization
    slab = readPOSCAR(exec_path / "POSCAR")
    assert len(slab.atlist) == len(ase_atoms.positions)


@pytest.mark.parametrize('file', ('BEAMLIST', 'VIBROCC', 'IVBEAMS'))
def test_refcalc_writes_file(run_from_ase_refcalc, file):
    """Ensure that run_from_ase writes work `file` during refcalc."""
    _, exec_path, _ = run_from_ase_refcalc
    assert (exec_path / "work" / file).is_file()


def test_refcalc_output_not_empty(run_from_ase_refcalc):
    """Ensure that run_from_ase ref-calc returns a non-empty output."""
    (theobeams, *_), *_ = run_from_ase_refcalc
    assert theobeams


def test_refcalc_output_correct_len(refcalc_thoeobeams):
    """Ensure ref-calc produced the correct number of beams."""
    assert len(refcalc_thoeobeams) == 4  # As per IVBEAMS


def test_refcalc_output_correct_energies(refcalc_thoeobeams):
    """Ensure ref-calc produced the correct energies."""
    assert all(refcalc_thoeobeams[0].energies == np.linspace(50, 70, 11))


def test_refcalc_output_non_nan(refcalc_thoeobeams):
    """Ensure ref-calc energies are numbers."""
    assert not any(np.isnan(refcalc_thoeobeams[0].energies))
