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
    )


@pytest.fixture(name="ase_ni_100_1x1_cell")
def fixture_ase_nickel_cell():
    """Return an ase.Atoms Ni(100)-1x1 with 6 layers."""
    element = 'Ni'
    cell_1x1 = ase.build.fcc110(element, size=(1,1,6), vacuum=3)
    return cell_1x1


# TODO: find a better way, perhaps a decorator that does
# the request + extracting fixture value.
@pytest.mark.parametrize("fixture, n_atoms", (("ase_ni_100_1x1_cell", 6),))
def test_ase_n_atoms(fixture, n_atoms, request):
    """Make sure `ase_atoms` has `n_atoms` atoms."""
    ase_atoms = request.getfixturevalue(fixture)
    assert len(ase_atoms.positions) == n_atoms


# TODO: will need to replace with Slab.from_ase method
def slab_from_ase(ase_atoms):
    """Return a Slab from an ase.Atoms object."""
    return Slab(ase_atoms)


@pytest.mark.parametrize("fixture", _ASE_ATOMS)
def test_n_atoms_from_ase(fixture, request):
    """Make sure the number of atoms in Slab match those in ase.Atoms."""
    ase_atoms = request.getfixturevalue(fixture)
    slab = slab_from_ase(ase_atoms)
    assert len(ase_atoms.positions) == len(slab.atlist)


def test_rotation_matrices():
    """Test that rotation matrices are as expected."""
    raise NotImplementedError


THETA = 14.7  # degrees

@pytest.mark.parametrize("fixture", _ASE_ATOMS)
def test_rot_mat_c(fixture, request):
    """Assert correct rotation around z axis."""
    ase_atoms = request.getfixturevalue(fixture)
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


class TestRaises:
    """Container for test that check exceptions."""

    @staticmethod
    def test_non_existing_exec_path():
        """Test exception for non existing execution path."""
        missing_path = INPUTS_ASE / "__th_is__do_es_not_ex_is_t__"
        with pytest.raises(ValueError) as exc:
            vpr_ase.run_from_ase(missing_path, None)
        assert exc.match("exec_path")

    @staticmethod
    def test_non_existing_parameters():
        """Test exception for non-existing PARAMETERS file."""
        with pytest.raises(RuntimeError) as exc:
            vpr_ase.run_from_ase(INPUTS_ASE, None)
        assert exc.match("PARAMETERS")

    @staticmethod
    def test_out_of_plane_ab(ase_ni_100_1x1_cell):
        """Test exception raised for a slab with non-xy a/b vectors."""
        rot_mat = vpr_ase.rot_mat_a(20)
        ase_atoms = ase_ni_100_1x1_cell
        exec_path = INPUTS_ASE / "initialization"  # Will not run
        with pytest.raises(RuntimeError) as exc:
            vpr_ase.run_from_ase(
                exec_path,
                ase_atoms,
                uc_transformation_matrix=rot_mat
                )
        assert exc.match("z component")

    # In principle we are also raising a RuntimeError in case
    # run_tleedm raises any exception. In practice, this should
    # not happen as the code currently is, since run_tleedm
    # swallows all exceptions....


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
    # The "initialization" folder contains only a PARAMETERS file,
    # but the required IVBEAMS or EXPBEAMS are not present, so the
    # next run_from_ase call will actually FAIL. This is fine though
    inputs_path = INPUTS_ASE / "initialization"
    results = vpr_ase.run_from_ase(
        exec_path=exec_path,
        ase_object=ase_atoms,
        inputs_path=inputs_path,
        cut_cell_c_fraction=0.0
        )
    return results, exec_path, ase_atoms


class TestFailingInitialization:
    """Tests for an "initialization" run that will fail.

    The failure comes from the fact the the "initialization" input
    directory is missing the required IVBEAMS file.
    """

    @pytest.fixture(autouse=True)                                              # TODO: would be nice to have scope="class", but tmp_path... needs scope="function"
    def run_init(self, run_from_ase_initialization):
        """Run the initialization once for the whole class."""
        (self.init_results,
         self.exec_path,
         self.ase_atoms) = run_from_ase_initialization

    def test_output_empty(self):
        """Ensure that run_from_ase initialization returns an empty output."""
        *all_except_v0i, _ = self.init_results
        assert all(not s for s in all_except_v0i)

    def test_returns_v0i(self):
        """Make sure the last return of run_from_ase is V0i."""
        *_, v0i = self.init_results
        assert isinstance(v0i, float)

    @pytest.mark.parametrize('file', ('POSCAR', 'work/VIBROCC'))
    def test_writes_file(self, file):
        """Ensure that run_from_ase writes `file` during initialization."""
        assert (self.exec_path / file).is_file()

    # TODO: perhaps it would be even better to store somewhere
    # checksums for files to be generated, and check them over here
    def test_writes_sensible_poscar(self):
        """Ensure that written POSCAR has the right number of atoms."""
        slab = readPOSCAR(self.exec_path / "POSCAR")
        assert len(slab.atlist) == len(self.ase_atoms.positions)


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


class TestSuccessfulRefcalc:
    """Tests for a "reference calculation" run with successful outcome."""

    @pytest.fixture(autouse=True, name="run_refcalc")                           # TODO: would be nice to have scope="class", but tmp_path... needs scope="function"
    def fixture_run_refcalc(self, run_from_ase_refcalc):
        """Run the ref-calc once for the whole class."""
        (self.refcalc_results,
         self.exec_path,
         self.ase_atoms) = run_from_ase_refcalc

    @pytest.fixture(autouse=True)                                               # TODO: would be nice to have scope="class", but tmp_path... needs scope="function"
    def read_theobeams_from_results(self, run_refcalc):
        """Store a list of full-dynamically calculated beams."""
        theobeams_content, *_ = self.refcalc_results
        self.theobeams = readOUTBEAMS(StringIO(theobeams_content))

    @pytest.mark.parametrize('file', ('BEAMLIST', 'VIBROCC', 'IVBEAMS'))
    def test_writes_file(self, file):
        """Ensure that run_from_ase writes work `file` during refcalc."""
        assert (self.exec_path / "work" / file).is_file()

    def test_output_not_empty(self):
        """Ensure that run_from_ase ref-calc returns a non-empty output."""
        theobeams_content, *_ = self.refcalc_results
        assert theobeams_content

    def test_correct_nr_beams(self):
        """Ensure ref-calc produced the correct number of beams."""
        assert len(self.theobeams) == 4  # As per IVBEAMS

    def test_correct_energies(self):
        """Ensure ref-calc produced the correct energies."""
        assert all(self.theobeams[0].energies == np.linspace(50, 70, 11))

    def test_no_nan_energies(self):
        """Ensure ref-calc energies are numbers."""
        assert not any(np.isnan(self.theobeams[0].energies))
