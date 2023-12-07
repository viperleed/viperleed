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

from ..poscar_slabs import POSCARS_WITHOUT_INFO, AG_100, SLAB_36C_cm

# pylint: disable=wrong-import-position
# Cannot do anything about it until we make viperleed installable
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


@pytest.mark.xfail(reason='updateElementCounts is buggy', strict=True)
def test_add_one_atom_n_elements():
    """Check that adding one atom to a slab updates elements correctly."""
    slab = Slab()
    new_atom = Atom('C', (0, 0, 0), 1, slab)
    slab.atlist.append(new_atom)
    slab.updateElementCounts()
    assert new_atom.el in slab.elements
    assert slab.n_per_elem[new_atom.el] == 1


class TestAtomsAndElements:
    """Collection of tests for atom additions/removals."""

    def test_empty_slab(self):
        """Check that an empty slab has no atoms, layers, etc..."""
        slab = Slab()
        assert slab.atlist == []
        assert slab.elements == ()
        assert slab.layers == []
        assert slab.planegroup == 'unknown'

    def test_add_one_atom_n_elements(self):
        """Check that adding one atom to a slab updates elements correctly."""
        slab = Slab()
        new_atom = Atom('C', (0, 0, 0), 1, slab)
        slab.atlist.append(new_atom)
        slab.update_element_count()
        assert new_atom.el in slab.elements
        assert slab.n_per_elem[new_atom.el] == 1

    def test_remove_one_atom_n_elements(self, make_poscar):
        """Check that removing one atom updates elements correctly."""
        slab, *_ = make_poscar(poscar_slabs.AG_100)
        n_ag_atoms = slab.n_per_elem['Ag']
        slab.atlist.pop()
        slab.update_element_count()
        assert slab.n_per_elem['Ag'] == n_ag_atoms - 1

    def test_slab_thickness(self, make_poscar):
        slab, *_ = make_poscar(AG_100)
        assert slab.thickness == pytest.approx(10.18233, abs=1e-4)

    def test_slab_vacuum_gap(self, make_poscar):
        slab, *_ = make_poscar(AG_100)
        assert slab.vacuum_gap == pytest.approx(10.18233, abs=1e-4)

    @parametrize(info=POSCARS_WITHOUT_INFO)
    def test_slab_sort_by_z(self, info, make_poscar):
        slab, *_ = make_poscar(info)
        slab.sort_by_z()
        assert all(at1.pos[2] <= at2.pos[2]
                   for at1, at2 in zip(slab.atlist, slab.atlist[1:]))

    def test_updateElementCount(self, make_poscar):
        slab, *_ = make_poscar(AG_100)
        n_ag_atoms = slab.n_per_elem['Ag']
        slab.atlist.pop()
        slab.updateElementCount()
        assert slab.n_per_elem['Ag'] == n_ag_atoms - 1


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


class TestSlabUcell:
    """Test for the Slab.ucell property."""
    def test_slab_ucell(self, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        assert slab.ucell.shape == (3, 3)

    def test_apply_scaling_scalar(self, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        slab.apply_scaling(2)
        assert slab.ucell == pytest.approx(np.array([[6, 0, 0],
                                                     [0, 8, 0],
                                                     [0, 0, 10]]))

    def test_apply_scaling_vector(self, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        slab.apply_scaling(1/3, 3.14, 0.1)
        assert slab.ucell == pytest.approx(np.array([[1, 0, 0],
                                                     [0, 12.56, 0],
                                                     [0, 0, 0.5]]))

    def test_angle_between_ucell_and_coord_sys_0(self, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        assert slab.angle_between_ucell_and_coord_sys == pytest.approx(0)

    def test_angle_between_ucell_and_coord_sys_30(self, manual_slab_3_atoms):
        slab = Slab()
        slab.ucell = np.array([[np.cos(np.pi/6), np.sin(np.pi/6), 0],
                               [np.cos(np.pi/3*4), np.sin(np.pi/3*4), 0],
                               [0, 0, 1]]).T
        assert slab.angle_between_ucell_and_coord_sys == pytest.approx(30)