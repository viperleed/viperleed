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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  IMPORTANT NOTICE: all the fixtures below are class-scoped. This means the  #
#  calculations will only run once per class. This also means that if new     #
#  tests are added that modify the objects, each of the test sets working     #
#  with one modified object should be collected into its own class.           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


INPUTS_ORIGIN = Path(__file__).parent / "fixtures"
INPUTS_ASE = INPUTS_ORIGIN / "from_ase"


def _make_refcalc_ok_transforms():
    """Yield slab transforms (or sequences thereof) and names (for ids)."""
    # A single transform-nothing, on its own and as a sequence:
    no_transform = vpr_ase.SlabTransform(cut_cell_c_fraction=0.)
    yield no_transform, "no_transform"
    yield (no_transform,), "no_transform, sequence"

    # A simple cut
    cut_only = vpr_ase.SlabTransform(cut_cell_c_fraction=0.2)
    yield cut_only, "cut only"

    # A 90-degrees in-plane rotation, no cutting
    rot_90_z = vpr_ase.SlabTransform(
        orthogonal_matrix=vpr_ase.rot_mat_z(90),
        cut_cell_c_fraction=0.
        )
    yield rot_90_z, "90deg z rotation"


def _make_refcalc_fail_transforms():
    """Yield slab transforms and names that cause the refcalc to fail."""
    # A single transform-nothing.
    # Fails because too few bulk layers.
    cut_default = vpr_ase.SlabTransform()
    yield cut_default, "cut default"

    # A 90-degrees in-plane rotation and cut default.
    # Fails because of too few bulk layers.
    rot_90_z = vpr_ase.SlabTransform(
        orthogonal_matrix=vpr_ase.rot_mat_z(90),
        )
    yield rot_90_z, "90deg z rotation, cut default"


_ASE_ATOMS = (
    "ase_ni_100_1x1_cell",
    )


@pytest.fixture(name="ase_ni_100_1x1_cell", scope="class")
def fixture_ase_nickel_cell():
    """Return an ase.Atoms Ni(100)-1x1 with 6 layers."""
    element = 'Ni'
    return ase.build.fcc110(element, size=(1,1,6), vacuum=3)


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


class TestRotationMatrices:
    """Test correctness of some simple rotation matrices."""

    _angles = (0, 12.3, 34, 95, 129.7, 256, 316)

    @staticmethod
    def test_rot_z_90():
        """Test 90deg rotation around z."""
        assert np.allclose(vpr_ase.rot_mat_z(90),
                           ((0, -1, 0), (1, 0, 0), (0, 0, 1)))

    @staticmethod
    def test_rot_x_90():
        """Test 90deg rotation around x."""
        assert np.allclose(vpr_ase.rot_mat_x(90),
                           ((1, 0, 0), (0, 0, -1), (0, 1, 0)))

    @staticmethod
    @pytest.mark.parametrize("theta", _angles)
    def test_rot_axis_x(theta):
        """Test that rotation around [1,0,0] is the same as around x."""
        assert np.allclose(vpr_ase.rot_mat_axis([1, 0, 0], theta),
                           vpr_ase.rot_mat_x(theta))

    @staticmethod
    @pytest.mark.parametrize("theta", _angles)
    def test_rot_axis_z(theta):
        """Test that rotation around [1,0,0] is the same as around x."""
        assert np.allclose(vpr_ase.rot_mat_axis([0, 0, 3], theta),
                           vpr_ase.rot_mat_z(theta))

    @staticmethod
    @pytest.mark.parametrize("theta", _angles)
    def test_apply_twice_rot_x(theta):
        """Ensure that rotating twice around x is the same as 2*theta."""
        _rot_theta = vpr_ase.rot_mat_x(theta)
        _rot_2theta = vpr_ase.rot_mat_x(2*theta)
        assert np.allclose(_rot_theta.dot(_rot_theta), _rot_2theta)

    @staticmethod
    @pytest.mark.parametrize("theta", _angles)
    def test_apply_twice_rot_axis(theta):
        """Ensure that rotating twice around x is the same as 2*theta."""
        _rot_theta = vpr_ase.rot_mat_axis((-4, 3, 12), theta)
        _rot_2theta = vpr_ase.rot_mat_axis((-4, 3, 12), 2*theta)
        assert np.allclose(_rot_theta.dot(_rot_theta), _rot_2theta)

    @staticmethod
    @pytest.mark.parametrize("theta", _angles)
    def test_orthogonal_rot_x(theta):
        """Ensure R(theta).T == inv(R(theta)) == R(-theta)."""
        _rot_theta = vpr_ase.rot_mat_x(theta)
        _rot_minus_theta = vpr_ase.rot_mat_x(-theta)
        assert np.allclose(_rot_theta.dot(_rot_theta.T), np.identity(3))
        assert np.allclose(_rot_theta.T, _rot_minus_theta)

    @staticmethod
    @pytest.mark.parametrize("theta", _angles)
    def test_orthogonal_rot_axis(theta):
        """Ensure R(theta).T == inv(R(theta)) == R(-theta)."""
        _rand_axis = np.random.rand(3)
        _rot_theta = vpr_ase.rot_mat_axis(_rand_axis, theta)
        _rot_minus_theta = vpr_ase.rot_mat_axis(_rand_axis, -theta)
        assert np.allclose(_rot_theta.dot(_rot_theta.T), np.identity(3))
        assert np.allclose(_rot_theta.T, _rot_minus_theta)


THETA = 14.7  # degrees

@pytest.mark.parametrize("fixture", _ASE_ATOMS)
def test_rot_mat_z(fixture, request):
    """Assert correct rotation around z axis."""
    ase_atoms = request.getfixturevalue(fixture)
    slab = slab_from_ase(ase_atoms)
    rot_mat = vpr_ase.rot_mat_z(THETA)
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
        with pytest.raises(FileNotFoundError) as exc:
            vpr_ase.run_from_ase(missing_path, None)
        assert exc.match("exec_path")

    @staticmethod
    def test_non_existing_parameters():
        """Test exception for non-existing PARAMETERS file."""
        with pytest.raises(FileNotFoundError) as exc:
            vpr_ase.run_from_ase(INPUTS_ASE, None)
        assert exc.match("PARAMETERS")

    @staticmethod
    def test_out_of_plane_ab(ase_ni_100_1x1_cell):
        """Test exception raised for a slab with non-xy a/b vectors."""
        transform = vpr_ase.SlabTransform(
            orthogonal_matrix=vpr_ase.rot_mat_x(20)
            )
        ase_atoms = ase_ni_100_1x1_cell
        exec_path = INPUTS_ASE / "initialization"  # Will not run
        with pytest.raises(ValueError) as exc:
            vpr_ase.run_from_ase(
                exec_path,
                ase_atoms,
                slab_transforms=transform
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
@pytest.fixture(name="run_from_ase_initialization", scope="class")
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
    no_cut = vpr_ase.SlabTransform(cut_cell_c_fraction=0.)
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
        slab_transforms=no_cut,
        )
    return results, exec_path, ase_atoms


class TestFailingInitialization:
    """Tests for an "initialization" run that will fail.

    The failure comes from the fact the the "initialization" input
    directory is missing the required IVBEAMS file.
    """

    @pytest.fixture(autouse=True)
    def run_init(self, run_from_ase_initialization):
        """Run the initialization once for the whole class."""
        # pylint: disable=attribute-defined-outside-init
        # Could in principle add an __init__ for these, then have a
        # custom test collection, e.g., as suggested in
        # https://github.com/pytest-dev/pytest/issues/7033
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


def make_refcalc_fixture(name, slab_transforms_and_ids, **kwargs):
    """Return a pytest.fixture for a refcalc with name and kwargs.

    Parameters
    ----------
    name : str
        The name to give to the fixture returned.
    slab_transforms_and_ids : Iterable
        Should yield pairs of the form (transforms, id) to be
        used for tests.
    kwargs : dict
        Keyword arguments to pass to the fixture decorator.

    Returns
    -------
    fixture : pytest.fixture
        The fixture that can be used for running a refcalc.
        The fixture is class-scoped by default. This can be
        changed by passing an appropriate keyword argument.
    """
    params_and_ids_dict = dict(zip(("params", "ids"),
                                   zip(*slab_transforms_and_ids)))
    params_and_ids_dict.update(kwargs)

    @pytest.fixture(name=name, scope="class", **params_and_ids_dict)
    def _fixture(ase_ni_100_1x1_cell, tmp_path_factory, request):
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
            slab_transforms=request.param,
            )
        return results, exec_path, ase_atoms
    return _fixture


fixture_run_from_ase_refcalc = make_refcalc_fixture(
    "run_from_ase_refcalc",
    _make_refcalc_ok_transforms()
    )


class TestSuccessfulRefcalc:
    """Tests for a "reference calculation" run with successful outcome."""

    @pytest.fixture(autouse=True, name="run_refcalc")
    def fixture_run_refcalc(self, run_from_ase_refcalc):
        """Run the ref-calc once for the whole class."""
        # pylint: disable=attribute-defined-outside-init
        # Could in principle add an __init__ for these, then have a
        # custom test collection, e.g., as suggested in
        # https://github.com/pytest-dev/pytest/issues/7033
        (self.refcalc_results,
         self.exec_path,
         self.ase_atoms) = run_from_ase_refcalc

    @pytest.fixture(autouse=True)
    def read_theobeams_from_results(self, run_refcalc):
        """Store a list of full-dynamically calculated beams."""
        # pylint: disable=attribute-defined-outside-init
        # See note in fixture_run_refcalc
        _ = run_refcalc  # Otherwise unused-argument
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


fixture_run_from_ase_refcalc_fails = make_refcalc_fixture(
    "run_from_ase_refcalc_fails",
    _make_refcalc_fail_transforms()
    )


class TestFailingRefcalc:
    """Tests for a "reference calculation" run with successful outcome."""

    @pytest.fixture(autouse=True, name="run_refcalc")
    def fixture_run_refcalc(self, run_from_ase_refcalc_fails):
        """Run the ref-calc once for the whole class."""
        # pylint: disable=attribute-defined-outside-init
        # Could in principle add an __init__ for these, then have a
        # custom test collection, e.g., as suggested in
        # https://github.com/pytest-dev/pytest/issues/7033
        (self.refcalc_results,
         self.exec_path,
         self.ase_atoms) = run_from_ase_refcalc_fails

    @pytest.fixture(autouse=True)
    def read_theobeams_from_results(self, run_refcalc):
        """Store an (empty) list of full-dynamically calculated beams."""
        # pylint: disable=attribute-defined-outside-init
        # See note in fixture_run_refcalc
        _ = run_refcalc  # Otherwise unused-argument
        theobeams_content, *_ = self.refcalc_results
        self.theobeams = readOUTBEAMS(StringIO(theobeams_content))

    @pytest.mark.parametrize('file', ('BEAMLIST', ))
    def test_does_not_write_file(self, file):
        """Ensure that run_from_ase does not write work `file`."""
        assert not (self.exec_path / "work" / file).is_file()

    def test_output_empty(self):
        """Ensure that run_from_ase ref-calc returns an empty output."""
        theobeams_content, *_ = self.refcalc_results
        assert not theobeams_content

    def test_no_beams(self):
        """Ensure ref-calc did not produce any beam."""
        assert not self.theobeams
