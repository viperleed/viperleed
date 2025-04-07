"""Tests for functionality in viperleed.calc.from_ase.

Define fixtures and test cases appropriate for the functionality
available in the from_ase module of viperleed.calc.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-02-23'
__license__ = 'GPLv3+'

from io import StringIO

import numpy as np
import pytest
from pytest_cases import fixture, parametrize_with_cases

from viperleed.calc import from_ase as vpr_ase
from viperleed.calc.classes.slab import Slab
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.files import poscar
from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.lib.math_utils import angle

from ..helpers import TEST_DATA
from . import cases_ase


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  IMPORTANT NOTICE: all the fixtures below are class-scoped. This means the  #
#  calculations will only run once per class. This also means that if new     #
#  tests are added that modify the objects, each of the test sets working     #
#  with one modified object should be collected into its own class.           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


ASE_DATA = TEST_DATA / 'from_ase'


def _make_refcalc_ok_transforms():
    """Yield slab transforms (or sequences thereof) and names (for ids)."""
    # A single transform-nothing, on its own and as a sequence:
    no_transform = vpr_ase.SlabTransform(cut_cell_c_fraction=0.)
    yield no_transform, 'no_transform'
    yield (no_transform,), 'no_transform, sequence'

    # A simple cut
    cut_only = vpr_ase.SlabTransform(cut_cell_c_fraction=0.2)
    yield cut_only, 'cut only'

    # A 90-degrees in-plane rotation, no cutting
    rot_90_z = vpr_ase.SlabTransform(
        orthogonal_matrix=vpr_ase.rot_mat_z(90),
        cut_cell_c_fraction=0.
        )
    yield rot_90_z, '90deg z rotation'


def _make_refcalc_fail_transforms():
    """Yield slab transforms and names that cause the refcalc to fail."""
    # A single transform-nothing.
    # Fails because too few bulk layers.
    cut_half = vpr_ase.SlabTransform(cut_cell_c_fraction=0.5)
    yield cut_half, 'cut half'

    # A 90-degrees in-plane rotation and cut default.
    # Fails because of too few bulk layers.
    rot_90_z_and_cut = vpr_ase.SlabTransform(
        orthogonal_matrix=vpr_ase.rot_mat_z(90),
        cut_cell_c_fraction=0.5
        )
    yield rot_90_z_and_cut, '90deg z rotation, cut half'


@parametrize_with_cases('ase_atoms,info', cases=cases_ase)
def test_ase_n_atoms(ase_atoms, info):
    """Make sure `ase_atoms` has `n_atoms` atoms."""
    assert len(ase_atoms.positions) == info.n_atoms


def slab_from_ase(ase_atoms):
    """Return a Slab from an ase.Atoms object."""
    return Slab.from_ase(ase_atoms)


@parametrize_with_cases('case', cases=cases_ase)
def test_n_atoms_from_ase(case):
    """Make sure the number of atoms in Slab match those in ase.Atoms."""
    ase_atoms, *_ = case
    slab = slab_from_ase(ase_atoms)
    assert len(ase_atoms.positions) == slab.n_atoms


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
    @pytest.mark.parametrize('theta', _angles)
    def test_rot_axis_x(theta):
        """Test that rotation around [1,0,0] is the same as around x."""
        assert np.allclose(vpr_ase.rot_mat_axis([1, 0, 0], theta),
                           vpr_ase.rot_mat_x(theta))

    @staticmethod
    @pytest.mark.parametrize('theta', _angles)
    def test_rot_axis_z(theta):
        """Test that rotation around [1,0,0] is the same as around x."""
        assert np.allclose(vpr_ase.rot_mat_axis([0, 0, 3], theta),
                           vpr_ase.rot_mat_z(theta))

    @staticmethod
    @pytest.mark.parametrize('theta', _angles)
    def test_apply_twice_rot_x(theta):
        """Ensure that rotating twice around x is the same as 2*theta."""
        _rot_theta = vpr_ase.rot_mat_x(theta)
        _rot_2theta = vpr_ase.rot_mat_x(2*theta)
        assert np.allclose(_rot_theta.dot(_rot_theta), _rot_2theta)

    @staticmethod
    @pytest.mark.parametrize('theta', _angles)
    def test_apply_twice_rot_axis(theta):
        """Ensure that rotating twice around x is the same as 2*theta."""
        _rot_theta = vpr_ase.rot_mat_axis((-4, 3, 12), theta)
        _rot_2theta = vpr_ase.rot_mat_axis((-4, 3, 12), 2*theta)
        assert np.allclose(_rot_theta.dot(_rot_theta), _rot_2theta)

    @staticmethod
    @pytest.mark.parametrize('theta', _angles)
    def test_orthogonal_rot_x(theta):
        """Ensure R(theta).T == inv(R(theta)) == R(-theta)."""
        _rot_theta = vpr_ase.rot_mat_x(theta)
        _rot_minus_theta = vpr_ase.rot_mat_x(-theta)
        assert np.allclose(_rot_theta.dot(_rot_theta.T), np.identity(3))
        assert np.allclose(_rot_theta.T, _rot_minus_theta)

    @staticmethod
    @pytest.mark.parametrize('theta', _angles)
    def test_orthogonal_rot_axis(theta):
        """Ensure R(theta).T == inv(R(theta)) == R(-theta)."""
        _rand_axis = np.random.rand(3)
        _rot_theta = vpr_ase.rot_mat_axis(_rand_axis, theta)
        _rot_minus_theta = vpr_ase.rot_mat_axis(_rand_axis, -theta)
        assert np.allclose(_rot_theta.dot(_rot_theta.T), np.identity(3))
        assert np.allclose(_rot_theta.T, _rot_minus_theta)


class TestSlabTransforms:
    """Test of simple slab transformations."""

    _theta = 14.7  # degrees
    _axis = np.random.rand(3)

    @fixture(name='slab')
    @parametrize_with_cases('case', cases=cases_ase)
    def fixture_slab(self, case):
        """Return a slab from ASE."""
        ase_atoms, *_ = case
        return slab_from_ase(ase_atoms)

    @staticmethod
    def parallel(vec_1, vec_2):
        """Return whether vec_1 and vec_2 are parallel."""
        return np.allclose(np.cross(vec_1, vec_2), 0)

    @staticmethod
    def normalized(vec):
        """Return a unit vector along vec."""
        return vec / np.linalg.norm(vec)

    @staticmethod
    def angle3(uvec1, uvec2):
        """Return the angle between two 3D unit-norm vectors."""
        return np.arccos(np.clip(np.dot(uvec1, uvec2), -1.0, 1.0))

    def test_rot_mat_z(self, slab):
        """Assert correct rotation around z axis."""
        rot_mat = vpr_ase.rot_mat_z(self._theta)
        ucell_before = slab.ucell.T.copy()
        slab.apply_matrix_transformation(rot_mat)
        ucell_after = slab.ucell.T

        for _before, _after in zip(ucell_before, ucell_after):
            assert np.isclose(_before[2], _after[2])  # Same z
            assert np.isclose(np.linalg.norm(_before), np.linalg.norm(_after))
            if all(self.parallel(v, (0, 0, 1)) for v in (_before, _after)):
                # Skip next for vectors along the rotation axis
                continue
            assert np.isclose(np.degrees(angle(_before[:2], _after[:2])),
                              self._theta)

    def test_rot_mat_x(self, slab):
        """Assert correct rotation around x axis."""
        rot_mat = vpr_ase.rot_mat_x(self._theta)
        ucell_before = slab.ucell.T.copy()
        slab.apply_matrix_transformation(rot_mat)
        ucell_after = slab.ucell.T

        for _before, _after in zip(ucell_before, ucell_after):
            assert np.isclose(_before[0], _after[0])  # Same x
            assert np.isclose(np.linalg.norm(_before), np.linalg.norm(_after))
            if all(self.parallel(v, (1, 0, 0)) for v in (_before, _after)):
                # Skip next for vectors along the rotation axis
                continue
            assert np.isclose(
                np.degrees(angle(_before[1:], _after[1:])),
                self._theta
                )

    def test_rot_mat_axis(self, slab):
        """Assert correct rotation around z axis."""
        axis = self._axis
        rot_mat = vpr_ase.rot_mat_axis(axis, self._theta)
        ucell_before = slab.ucell.T.copy()
        slab.apply_matrix_transformation(rot_mat)
        ucell_after = slab.ucell.T

        for _before, _after in zip(ucell_before, ucell_after):
            assert np.isclose(_before.dot(axis), _after.dot(axis))
            assert np.isclose(np.linalg.norm(_before), np.linalg.norm(_after))
            if all(self.parallel(v, axis) for v in (_before, _after)):
                # Skip next for vectors along the rotation axis
                continue
            _perp_before = self.normalized(np.cross(_before, axis))
            _perp_after = self.normalized(np.cross(_after, axis))
            assert np.isclose(
                np.degrees(self.angle3(_perp_before, _perp_after)),
                self._theta
                )

    @staticmethod
    def todo_test_swap_axis_transform():
        """Apply a swap transform and verify correctness."""
        # Should use the _apply_transform function
        raise NotImplementedError

    @staticmethod
    def todo_test_apply_two_transforms():
        """Apply consecutive transforms, and verify correctness."""
        # Should use the _apply_transform function.
        # Ideas:
        # - rotate twice same angle same axis, verify with single
        # - swap twice same --> unchanged
        # - swap cyclically
        raise NotImplementedError


class TestRaises:
    """Container for test that check exceptions."""

    @staticmethod
    def test_non_existing_exec_path():
        """Test exception for non existing execution path."""
        missing_path = ASE_DATA / '__th_is__do_es_not_ex_is_t__'
        with pytest.raises(FileNotFoundError) as exc:
            vpr_ase.run_from_ase(missing_path, None)
        assert exc.match('exec_path')

    @staticmethod
    def test_non_existing_parameters():
        """Test exception for non-existing PARAMETERS file."""
        with pytest.raises(FileNotFoundError) as exc:
            vpr_ase.run_from_ase(ASE_DATA, None)
        assert exc.match('PARAMETERS')

    @parametrize_with_cases('args', cases=cases_ase)
    def test_out_of_plane_ab(self, args):
        """Test exception raised for a slab with non-xy a/b vectors."""
        transform = vpr_ase.SlabTransform(
            orthogonal_matrix=vpr_ase.rot_mat_x(20)
            )
        ase_atoms, *_ = args
        exec_path = ASE_DATA / 'initialization'  # Will not run
        with pytest.raises(ValueError) as exc:
            vpr_ase.run_from_ase(
                exec_path,
                ase_atoms,
                slab_transforms=transform
                )
        assert exc.match('z component')

    # In principle we are also raising a RuntimeError in case
    # run_calc raises any exception. In practice, this should
    # not happen as the code currently is, since run_calc
    # swallows all exceptions....


# TODO: find a nice way to generate multiple fixtures dynamically.
# See also notes in helpers.py. Difficulties:
# - make temp path dependent on the name of the first argument in
#   the signature
# - pick the right folder in 'initialization'
@fixture(name='run_from_ase_initialization', scope='class')
@parametrize_with_cases('case', cases=cases_ase, scope='class')
def fixture_run_from_ase_initialization(case, tmp_path_factory, request):
    """Return the results of an initialization run.

    Parameters
    ----------
    case : tuple
        Only the first item is used. It is the ase.Atoms object
        from which to run.
    tmp_path_factory : pytest.fixture
        The pytest default name for the temporary directory factory.
    request : pytest.fixture
        The current request fixture.

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
    ase_atoms, *_ = case
    no_cut = vpr_ase.SlabTransform(cut_cell_c_fraction=0.)
    case_id = request.param.argvalues[0]
    exec_path = tmp_path_factory.mktemp(basename=f'from_{case_id}_init',
                                        numbered=True)
    # The "initialization" folder contains only a PARAMETERS file,
    # but the required IVBEAMS or EXPBEAMS are not present, so the
    # next run_from_ase call will actually FAIL. This is fine though
    inputs_path = ASE_DATA / 'initialization'
    results = vpr_ase.run_from_ase(
        exec_path=exec_path,
        ase_object=ase_atoms,
        inputs_path=inputs_path,
        slab_transforms=no_cut,
        )
    return results, exec_path, ase_atoms


class TestFailingInitialization:
    """Tests for an 'initialization' run that will fail.

    The failure comes from the fact the the 'initialization' input
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

    @pytest.mark.parametrize('file', ('POSCAR', f'{DEFAULT_WORK}/VIBROCC'))
    def test_writes_file(self, file):
        """Ensure that run_from_ase writes `file` during initialization."""
        assert (self.exec_path / file).is_file()

    # TODO: perhaps it would be even better to store somewhere
    # checksums for files to be generated, and check them over here
    def test_writes_sensible_poscar(self):
        """Ensure that written POSCAR has the right number of atoms."""
        slab = poscar.read(self.exec_path / 'POSCAR')
        assert slab.n_atoms == len(self.ase_atoms.positions)


def make_refcalc_fixture(name, slab_transforms_and_ids,
                         accept_fails=False, **kwargs):                         # TODO: this one could be @parametrize and un-indented.
    """Return a pytest.fixture for a refcalc with name and kwargs.

    Parameters
    ----------
    name : str
        The name to give to the fixture returned.
    slab_transforms_and_ids : Iterable
        Should yield pairs of the form (transforms, id) to be
        used for tests.
    accept_fails : bool, optional
        Do not complain if the viperleed.calc execution fails.
        Default is False.
    kwargs : dict
        Keyword arguments to pass to the fixture decorator.

    Returns
    -------
    fixture : pytest.fixture
        The fixture that can be used for running a refcalc.
        The fixture is class-scoped by default. This can be
        changed by passing an appropriate keyword argument.
    """
    params_and_ids_dict = dict(zip(('params', 'ids'),
                                   zip(*slab_transforms_and_ids)))
    params_and_ids_dict.update(kwargs)

    @pytest.fixture(name=name, scope='class', **params_and_ids_dict)            # TODO: somehow here @parametrize_with_cases fails. May be due to the weird placement of this fixture.
    def _fixture(tmp_path_factory, request):
        """Return the results of a reference calculation run.

        Parameters
        ----------
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
        ase_atoms, *_ = cases_ase.case_ase_ni_100_1x1_cell()
        exec_path = tmp_path_factory.mktemp(basename='from_ase_Ni_100_refcalc',
                                            numbered=True)
        inputs_path = ASE_DATA / 'refcalc'
        kwargs = {
            'exec_path': exec_path,
            'ase_object': ase_atoms,
            'inputs_path': inputs_path,
            'slab_transforms': request.param,
            }
        try:
            results = vpr_ase.run_from_ase(**kwargs)
        except RuntimeError as exc:
            calc_failed = 'ViPErLEED calculation failed' in str(exc)
            if not accept_fails or not calc_failed:
                raise
            results = '', '', '', 0
        return results, exec_path, ase_atoms
    return _fixture


fixture_run_from_ase_refcalc = make_refcalc_fixture(
    'run_from_ase_refcalc',
    _make_refcalc_ok_transforms()
    )


class TestSuccessfulRefcalc:
    """Tests for a "reference calculation" run with successful outcome."""

    @pytest.fixture(autouse=True, name='run_refcalc')
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
    @pytest.mark.usefixtures('run_refcalc')
    def read_theobeams_from_results(self):
        """Store a list of full-dynamically calculated beams."""
        theobeams_content, *_ = self.refcalc_results
        # pylint: disable-next=attribute-defined-outside-init
        self.theobeams = readOUTBEAMS(StringIO(theobeams_content))

    @pytest.mark.parametrize('file', ('BEAMLIST', 'VIBROCC', 'IVBEAMS'))
    def test_writes_file(self, file):
        """Ensure that run_from_ase writes work `file` during refcalc."""
        assert (self.exec_path / DEFAULT_WORK / file).is_file()

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
    'run_from_ase_refcalc_fails',
    _make_refcalc_fail_transforms(),
    accept_fails=True,
    )


class TestFailingRefcalc:
    """Tests for a "reference calculation" run with successful outcome."""

    @pytest.fixture(autouse=True, name='run_refcalc')
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
    @pytest.mark.usefixtures('run_refcalc')
    def read_theobeams_from_results(self):
        """Store an (empty) list of full-dynamically calculated beams."""
        theobeams_content, *_ = self.refcalc_results
        # pylint: disable-next=attribute-defined-outside-init
        self.theobeams = readOUTBEAMS(StringIO(theobeams_content))

    @pytest.mark.parametrize('file', ('BEAMLIST', ))
    def test_does_not_write_file(self, file):
        """Ensure that run_from_ase does not write work `file`."""
        assert not (self.exec_path / DEFAULT_WORK / file).is_file()

    def test_output_empty(self):
        """Ensure that run_from_ase ref-calc returns an empty output."""
        theobeams_content, *_ = self.refcalc_results
        assert not theobeams_content

    def test_no_beams(self):
        """Ensure ref-calc did not produce any beam."""
        assert not self.theobeams
