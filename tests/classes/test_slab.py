"""Tests for viperleed.tleedmlib.classes.slab.

Created on 2023-07-28

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

from copy import deepcopy
from pathlib import Path
import sys

import numpy as np
import pytest
from pytest_cases import parametrize

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

from .. import poscar_slabs
from ..poscar_slabs import POSCARS_WITHOUT_INFO, AG_100, SLAB_36C_cm

# pylint: disable=wrong-import-position
# Cannot do anything about it until we make viperleed installable
#from viperleed.tleedmlib.base import pairwise
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.slab import Slab, SymPlane
# pylint: enable=wrong-import-position


class TestAtomTransforms:
    """Test simple transformations of the atoms of a slab."""

    def test_mirror(self, manual_slab_3_atoms):
        """Test the expected outcome of mirroring atoms of a simple slab."""
        slab = manual_slab_3_atoms
        mirrored_slab = deepcopy(slab)
        symplane = SymPlane((0, 0), (0, 1), abt=slab.surface_vectors)
        mirrored_slab.mirror(symplane)
        mirrored_slab.collapseCartesianCoordinates()
        assert all(
            at.is_same_xy(mir_at)
            for at, mir_at in zip(slab.atlist, reversed(mirrored_slab.atlist))
            )

    def test_mirror_on_slanted_cell(self, make_poscar):
        """Test the expected outcome of mirroring atoms of a slanted slab."""
        slab = make_poscar(SLAB_36C_cm)[0]
        sym_plane = SymPlane((0, 0), (0, 1), abt=slab.surface_vectors)
        assert slab.isMirrorSymmetric(sym_plane, eps=1e-6)  # will be is_mirror_symmetric on refactor branch

    def test_180_rotation(self, manual_slab_3_atoms):
        """Test the expected outcome of rotating atoms of a simple slab."""
        slab = manual_slab_3_atoms
        rotated_slab = deepcopy(slab)
        rotated_slab.rotateAtoms((0, 0), order=2)
        rotated_slab.collapseCartesianCoordinates()
        assert all(
            at.is_same_xy(rot_at)
            for at, rot_at in zip(slab.atlist, reversed(rotated_slab.atlist))
            )

    def test_180_rotation_on_slanted_cell(self, make_poscar):
        """Test the expected outcome of rotating atoms of a slanted slab."""
        slab = make_poscar(SLAB_36C_cm)[0]
        assert slab.isRotationSymmetric(np.array([0,0]), order=2, eps=1e-6)  # will be is_rotation_symmetric on refactor branch


class TestUnitCellTransforms:
    """Test simple transformations of the unit cell of a slab."""

    def test_rotation_on_trigonal_slab(self, manual_slab_1_atom_trigonal):
        """Test application of a rotation to a trigonal slab."""
        rot_15 = [[0.96592583, -0.25881905,  0.],
                  [0.25881905,  0.96592583,  0.],
                  [0.,  0.,  1.]]
        expected_cell = [[0.96592583,  0.25881905,  0.],
                         [-2.99808654, 2.30249368,  0.],
                         [0.44828774,  2.1906707,   3.]]
        expected_atom_cartpos = [-1.86064664,  1.88257645]
        slab = manual_slab_1_atom_trigonal
        slab.apply_matrix_transformation(rot_15)
        assert np.allclose(slab.ucell.T, expected_cell)
        assert np.allclose(slab.atlist[0].cartpos[:2], expected_atom_cartpos)


# This class currently contains NO TESTS, as the previous tests were
# somewhat faulty. The purpose of restore_ori_state is to:
#    (i) convert the current positions, vibrations and occupations
#        into VIBROCC offsets, and
#   (ii) fully clear the displacements of all atoms
# The tests were originally set up by @amimre, who found a bug that
# concerns runs where multiple refcalc-search pairs exists in RUN.
# This bug is documented in Issue                                               # TODO @amimre: refer to relevant issue here
class TestRestoreOristate:
    """Collection of tests for reverting a slab to its ref-calc state."""

    @pytest.fixture(name='slab_and_copy')
    def fixture_slab_and_copy(self, ag100_with_displacements_and_offsets):
        """Return a Ag(100) slab and its deepcopy."""
        slab, *_ = ag100_with_displacements_and_offsets
        return slab, deepcopy(slab)

    @staticmethod
    def check_displacements_equal(atom1, atom2, which, element, subtests):
        """Check that two atoms have the same requested displacements."""
        _test_str = f'{atom1}, {element}, '
        displ_1 = atom1.displacements[which]
        displ_2 = atom2.displacements[which]
        with subtests.test(_test_str + 'total'):
            assert displ_1[element] == pytest.approx(displ_2[element])

        offs_1 = displ_1.vibrocc_offset
        offs_2 = displ_2.vibrocc_offset
        with subtests.test(_test_str + 'vibrocc'):
            assert offs_1[element] == pytest.approx(offs_2[element])

        if which == 'geo':  # Check also offset
            offs_1 = displ_1.displacements_offset
            offs_2 = displ_2.displacements_offset
            with subtests.test(_test_str + 'geo offset'):
                assert offs_1[element] == pytest.approx(offs_2[element])


class TestAtomsAndElements:
    """Collection of tests for atom additions/removals."""

    def test_empty_slab(self):
        """Check that an empty slab has no atoms, layers, etc..."""
        slab = Slab()
        assert slab.atlist == []
        assert slab.elements == ()
        assert slab.layers == []
        assert slab.planegroup == 'unknown'

    @pytest.mark.xfail(reason='updateElementCounts is buggy', strict=True)
    def test_add_one_atom_n_elements(self):
        """Check that adding one atom to a slab updates elements correctly."""
        slab = Slab()
        new_atom = Atom('C', (0, 0, 0), 1, slab)
        slab.atlist.append(new_atom)
        slab.update_element_count()
        assert new_atom.el in slab.elements
        assert slab.n_per_elem[new_atom.el] == 1

    @pytest.mark.xfail(reason='updateElementCounts is buggy', strict=True)
    def test_remove_one_atom_n_elements(self, make_poscar):
        """Check that removing one atom updates elements correctly."""
        slab, *_ = make_poscar(poscar_slabs.AG_100)
        n_ag_atoms = slab.n_per_elem['Ag']
        slab.atlist.pop()
        slab.update_element_count()
        assert slab.n_per_elem['Ag'] == n_ag_atoms - 1


    @pytest.mark.skip(reason='to be implemented')
    def test_atlist_is_not_list(self):                                          # TODO: Should consider various situations to make sure that no Slab method messes with the atlist
        """TODO"""

    def test_update_atom_numbers(self, make_poscar):
        slab, *_ = make_poscar(AG_100)
        n_ag_atoms = slab.n_per_elem['Ag']
        slab.atlist.pop()
        slab.updateElementCount()
        assert slab.n_per_elem['Ag'] == n_ag_atoms - 1

    @pytest.mark.skip(reason='to be implemented')
    def test_chemelem_upon_element_mix_changed(self):
        """Check that chemelem are changed when ELEMENT_MIX is."""


@pytest.mark.skip(reason='to be implemented')
class TestBulk3DOperations:
    """Tests for 3D symmetry operations of BulkSlab objects."""

    def test_get_candidate_layer_periods(self):
        """TODO"""

    def test_bulk_screw_symmetric(self):                                        # TODO: especially consider the cases where an should be on an axis but it is away, and those where it should be at an n-fold position but it isn't.
        """TODO"""

    def test_bulk_glide_symmetric(self):                                        # TODO: especially consider the cases where an should be on an plane but it is away, and those where it should be at a glide-symmetric position but it isn't.
        """TODO"""


@pytest.mark.skip(reason='to be implemented')
class TestBulkDetectAndExtraBulk:
    """Collection of tests for adding bulk units to slabs."""

    def test_detect_bulk(self):                                                 # TODO: also check that rp and sl are unchanged if it fails
        """TODO"""

    def test_with_extra_bulk_units(self):                                       # TODO: also check the number of bulk layers
        """TODO"""

    def test_with_double_thickness_once(self):
        """TODO"""

    def test_with_double_thickness_twice(self):
        """TODO"""


@pytest.mark.skip(reason='to be implemented')
class TestBulkUcell:
    """Tests concerning reduction of bulk unit cell and C vector."""

    def test_get_min_c(self):                                                   # TODO: also various raises
        """TODO"""

    def test_ensure_min_c(self):
        """TODO"""

    def test_apply_bulk_ucell_reduction(self):                                  # TODO: separately ab only and c_vec. Both with recenter==True/False
        """TODO"""

    def test_minimal_bulk_ab(self):                                             # TODO: SurfaceSlab
        """TODO"""


class TestCoordinates:
    """Collection of tests for Cartesian/fractional coordinates."""
    def test_cartesian_from_fractional(self, manual_slab_3_atoms):
        """Check correct update of Cartesian atom coordinates."""
        slab = manual_slab_3_atoms
        atom = slab.atlist[0]
        atom.pos = np.array([0.1, 0.2, 0.3])
        slab.getCartesianCoordinates()
        assert atom.cartpos == pytest.approx([0.3, 0.8, -1.5])

    def test_cartesian_from_fractional_with_origin(self, manual_slab_3_atoms):
        """Check correct update of Cartesian atom coordinates with updating origin."""
        slab = manual_slab_3_atoms
        atom = slab.atlist[0]
        atom.pos = np.array([0.1, 0.2, 0.3])
        slab.getCartesianCoordinates(updateOrigin=True)
        assert atom.pos == pytest.approx([0.1, 0.2, 0.3])
        assert atom.cartpos == pytest.approx([0.3, 0.8, 0])

    def test_fractional_from_cartesian(self, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        slab.atlist[0].cartpos = np.array([0.3, 0.8, -1.5])
        slab.getFractionalCoordinates()
        assert slab.atlist[0].pos == pytest.approx([0.1, 0.2, 0.3])

    @pytest.mark.xfail(reason='collapseCartesianCoordinates is buggy??')
    def test_collapse_cartesian(self, manual_slab_3_atoms):
        """Check that cartesian coordinates are correctly collapsed."""
        slab = manual_slab_3_atoms
        slab.atlist[0].cartpos = np.array([5, 6, 3])
        slab.collapseCartesianCoordinates()
        assert slab.atlist[0].cartpos == pytest.approx([2, 2, 3])

    def test_collapse_fractional(self, manual_slab_3_atoms):
        """Check that fractional coordinates are correctly collapsed."""
        slab = manual_slab_3_atoms
        slab.atlist[0].pos = np.array([1.1, 2.2, -3.3])
        slab.collapseFractionalCoordinates()
        assert slab.atlist[0].pos == pytest.approx([0.1, 0.2, 0.7])

    def test_collapse_fractional_small_eps(self, manual_slab_3_atoms):
        """Check that fractional coordinates are correctly collapsed even with small eps."""
        slab = manual_slab_3_atoms
        slab.atlist[0].pos = np.array([1.0 - 1e-9, 2.0 + 1e-9, -3.0 -1e-15])
        slab.collapseFractionalCoordinates()
        assert slab.atlist[0].pos == pytest.approx([1, 0.0, 1.0], abs = 1e-8)


@pytest.mark.skip(reason='to be implemented')
class TestDuplicateAtoms:
    """Tests for checking detection and removal of duplicate atoms."""

    @pytest.mark.parametrize('info', poscar_slabs.WITH_DUPLICATE_ATOMS)
    def test_with_duplicate_atoms(self, info, make_poscar):
        """Check that POSCARs with duplicates are handled correctly."""
        slab, *_ = make_poscar(info)
        with pytest.raises(AtomsTooCloseError):
            slab.check_atom_collisions()

    def test_without_duplicates(self, make_poscar):
        """Check that POSCARs without duplicates are handled correctly."""
        slab, *_ = make_poscar(poscar_slabs.AG_100)
        with not_raises(AtomsTooCloseError):
            slab.check_atom_collisions()

    @pytest.mark.parametrize('info', poscar_slabs.WITH_DUPLICATE_ATOMS)
    def test_remove_duplicates(self, info):                                     # TODO: n_atoms, raises, others? check method
        """TODO"""
        slab, rpars, *_ = make_poscar(info)


@pytest.mark.skip(reason='to be implemented')
class TestSuperAndSubCell:
    """Collection of tests for creation of larger and smaller slab versions."""

    def test_supercell(self):                                                   # TODO: diagonal and non-diagonal (for some weird basis cell?). I think the old version was failing under some non-diagonal situations. Explicitly test the two removed update_origin.
        """TODO"""

    def test_subcell(self):                                                     # TODO: this is the inverse of the previous one.
        """TODO"""


class TestSlabLayers:
    """Collection of tests concerning slab (sub)layers."""

    @pytest.mark.skip(reason='to be implemented')
    def test_bulk_layers(self):
        """TODO"""

    @pytest.mark.skip(reason='to be implemented')
    def test_create_layers(self):
        """Check that layers are created correctly."""

    @pytest.mark.skip(reason='to be implemented')
    def test_create_sublayers(self):                                            # TODO: also test if this works fine excluding the second sort-by-element run
        """Check that sublayers are created correctly."""

    @pytest.mark.skip(reason='to be implemented')
    def test_full_update_with_layers_defined(self):                             # TODO: test correct behaviour for (i) coords initially outside the unit cell, and (ii) no topat_ori_z available
        """TODO"""

    @pytest.mark.skip(reason='to be implemented')
    def test_interlayer_spacing(self):                                          # TODO: also raises.
        """TODO"""

    @pytest.mark.skip(reason='to be implemented')
    def test_slab_lowocc_sublayer(self):
        """TODO"""


class TestSorting:
    """Collection of tests for slab sorting."""

    @pytest.mark.skip(reason='to be implemented')
    def test_element_sort(self):
        """Check correct element-based sorting of a Slab."""

    @pytest.mark.skip(reason='to be implemented')
    def test_element_sort_raises_with_outdated_elements(self):
        """Ensure sorting complains when elements are outdated."""

    @pytest.mark.skip(reason='to be implemented')
    def test_sort_original(self):
        """Check correct sorting of a Slab to original atom numbers."""

    @parametrize(info=POSCARS_WITHOUT_INFO)
    def test_simple_sort_by_z(self, info, make_poscar):
        slab, *_ = make_poscar(info)
        slab.sort_by_z()
        assert all(at1.pos[2] <= at2.pos[2]
                   for at1, at2 in zip(slab.atlist, slab.atlist[1:]))

    @pytest.mark.skip(reason='functions not implemented on master')
    @parametrize(info=poscar_slabs.POSCARS_WITHOUT_INFO)
    @parametrize(bottom_to_top=(True, False))
    def test_z_sort(self, info, bottom_to_top, make_poscar):
        """Check successful sorting of atoms by out-of-plane position."""
        slab, *_ = make_poscar(info)
        slab.sort_by_z(bottom_to_top=bottom_to_top)
        _ordered = operator.ge if bottom_to_top else operator.le                # TODO: swap when flipping .cartpos
        assert all(_ordered(at1.cartpos[2], at2.cartpos[2])
                   for at1, at2 in pairwise(slab))


class TestUnitCellTransforms:
    """Test simple transformations of the unit cell of a slab."""

    def test_ucell_array(self, manual_slab_3_atoms):
        """Check that the array shape of the unit cell is as expected."""
        assert manual_slab_3_atoms.ucell.shape == (3, 3)

    def test_angle_between_ucell_and_coord_sys_zero(self, manual_slab_3_atoms):
        """Check correct identification of the rotation of the a vector."""
        slab = manual_slab_3_atoms
        assert slab.angle_between_ucell_and_coord_sys == pytest.approx(0)

    def test_angle_between_ucell_and_coord_sys_30(self):
        """Check correct identification of the rotation of the a vector."""
        slab = Slab()
        slab.ucell = np.array([[np.cos(np.pi/6), np.sin(np.pi/6), 0],
                               [np.cos(np.pi/3*4), np.sin(np.pi/3*4), 0],
                               [0, 0, 1]]).T
        assert slab.angle_between_ucell_and_coord_sys == pytest.approx(30)

    _scalings = {  # Applied, expected result
        'scalar': ((2,),  np.diag((6, 8, 10))),
        'vector': ((1/3, 3.14, 0.1), np.diag((1, 12.56, 0.5)))
        }

    @parametrize('scaling,expected', _scalings.values(), ids=_scalings)
    def test_apply_scaling(self, scaling, expected, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        slab.apply_scaling(*scaling)
        assert slab.ucell == pytest.approx(expected)

    @pytest.mark.skip(reason='to be implemented')
    def test_project_c_to_z(self):
        """TODO"""

    @pytest.mark.skip(reason='to be implemented')
    def test_revert_unit_cell(self):                                            # TODO: Probably best to pick a few random operations and make sure that reverting one+rest, a few+rest, or all of them at once gives the same result. This should include unit cell as well as all atom frac and cart coordinates. Also test raises RuntimeError.
        """TODO"""

    @pytest.mark.skip(reason='to be implemented')
    def test_revert_unit_cell_undo_nothing(self):                               # TODO: both by having nothing to undo, and by passing as many as there are operations. Check especially by manually translating atoms out of the base cell.
        """TODO"""

    def test_rotation_on_trigonal_slab(self, manual_slab_1_atom_trigonal):
        """Test application of a rotation to a trigonal slab."""
        rot_15 = [[0.96592583, -0.25881905,  0.],
                  [0.25881905,  0.96592583,  0.],
                  [0.,  0.,  1.]]
        expected_cell = [[0.96592583,  0.25881905,  0.],
                         [-2.99808654, 2.30249368,  0.],
                         [0.44828774,  2.1906707,   3.]]
        expected_atom_cartpos = [-1.86064664,  1.88257645]
        slab = manual_slab_1_atom_trigonal
        slab.apply_matrix_transformation(rot_15)
        assert np.allclose(slab.ucell.T, expected_cell)
        assert np.allclose(slab.atlist[0].cartpos[:2], expected_atom_cartpos)

class TestSlabProperties:
    def test_slab_thickness(self, make_poscar):
        slab, *_ = make_poscar(AG_100)
        assert slab.thickness == pytest.approx(10.18233, abs=1e-4)

    def test_slab_vacuum_gap(self, make_poscar):
        slab, *_ = make_poscar(AG_100)
        assert slab.vacuum_gap == pytest.approx(10.18233, abs=1e-4)

@pytest.mark.skip(reason='to be implemented')
class TestUnitCellReduction:
    """Tests for minimization of various bits of the unit cell."""

    def test_ab_cell_minimization(self):                                        # TODO: Use also "Sb on Si(111)" case from Max Buchta
        """TODO"""

    def test_ab_cell_already_minimal(self):
        """TODO"""


@pytest.mark.skip(reason='to be implemented')
def test_contains():
    """TODO"""


@pytest.mark.skip(reason='to be implemented')
def test_check_ab_in_plane():
    """TODO"""


@pytest.mark.skip(reason='to be implemented')
def test_nearest_neighbors():
    """TODO"""


@pytest.mark.skip(reason='to be implemented')
def test_ucell_ori_after_reset_symmetry():                                      # TODO: ucell_ori should not change when transforming the unit cell after resetSymmetry was called
    """TODO"""


@pytest.mark.skip(reason='to be implemented')
def test_translation_symmetry_different_species():                              # TODO: could use something like an MgO slab and ascertain that Mg->O is not a valid translation
    """TODO"""


@pytest.mark.skip(reason='to be implemented')
def test_slab_equivalence():                                                    # TODO: check also cases covered by TODOs
    """TODO"""


@pytest.mark.skip(reason='to be implemented')
def test_identify_bulk_repeat():
    """TODO"""


@pytest.mark.skip(reason='to be implemented')
def test_get_bulk_repeat():
    """TODO"""


@pytest.mark.skip(reason='to be implemented')
def test_make_bulk_slab():
    """TODO"""
