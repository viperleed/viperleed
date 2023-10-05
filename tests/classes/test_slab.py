"""Tests for viperleed.tleedmlib.classes.slab.

Created on 2023-07-28

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

from copy import deepcopy
import operator
from pathlib import Path
import sys

import numpy as np
import pytest
from pytest_cases import parametrize

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

from ..poscar_slabs import POSCARS_WITHOUT_INFO, AG_100

# pylint: disable=wrong-import-position
# Cannot do anything about it until we make viperleed installable
from viperleed.tleedmlib.base import pairwise
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.slab import Slab
from viperleed.tleedmlib.classes.sym_entity import SymPlane
# pylint: enable=wrong-import-position


class TestAtomTransforms:                                                       # TODO: add test for correct rotation and mirror on slanted cell. Probably enough to is_rotation_symmetric and is_mirror_symmetric with appropriate slabs. Could use "POSCAR_36C_cm" or the bulk of Fe3O4.
    """Test simple transformations of the atoms of a slab."""

    def test_mirror(self, manual_slab_3_atoms):
        """Test the expected outcome of mirroring atoms of a simple slab."""
        slab = manual_slab_3_atoms
        mirrored_slab = deepcopy(slab)
        symplane = SymPlane((0, 0), (0, 1), abt=slab.ab_cell.T)
        mirrored_slab.mirror_atoms(symplane)
        assert all(
            at.isSameXY(mir_at.cartpos[:2])
            for at, mir_at in zip(slab, reversed(mirrored_slab))
            )

    def test_180_rotation(self, manual_slab_3_atoms):
        """Test the expected outcome of rotating atoms of a simple slab."""
        slab = manual_slab_3_atoms
        rotated_slab = deepcopy(slab)
        rotated_slab.rotate_atoms(order=2)
        assert all(
            at.isSameXY(rot_at.cartpos[:2])
            for at, rot_at in zip(slab, reversed(rotated_slab))
            )


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
        slab, *_ = make_poscar(AG_100)
        n_ag_atoms = slab.n_per_elem['Ag']
        slab.atlist.pop()
        slab.update_element_count()
        assert slab.n_per_elem['Ag'] == n_ag_atoms - 1


class TestCoordinates:
    """Collection of tests for Cartesian/fractional coordinates."""

    def test_cartesian_from_fractional(self, manual_slab_3_atoms):              # TODO: also update_origin
        """Check correct update of Cartesian atom coordinates."""
        slab = manual_slab_3_atoms
        atom = slab.atlist[0]
        atom.pos = np.array([0.1, 0.2, 0.3])
        slab.update_cartesian_from_fractional()
        assert atom.cartpos == pytest.approx([0.3, 0.8, -1.5])


    def test_collapse_fractional(self, manual_slab_3_atoms):                    # TODO: use both methods. Find especially cases that are 'problematic' with %1.0: e.g., 1-1e-8, 1-1e-9, 1-1e-15
        """Check that fractional coordinates are correctly collapsed."""
        slab = manual_slab_3_atoms
        atom = slab.atlist[0]
        atom.pos = np.array([1.1, 2.2, -3.3])
        slab.collapse_fractional_coordinates()
        assert atom.pos == pytest.approx([0.1, 0.2, 0.7])

    def test_fractional_from_cartesian(self, manual_slab_3_atoms):
        """Check correct update of fractional atom coordinates."""
        slab = manual_slab_3_atoms
        atom = slab.atlist[0]
        atom.cartpos = np.array([0.3, 0.8, -1.5])
        slab.update_fractional_from_cartesian()
        assert atom.pos == pytest.approx([0.1, 0.2, 0.3])


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



class TestSorting:
    """Collection of tests for slab sorting."""

    @parametrize(info=POSCARS_WITHOUT_INFO)
    @parametrize(bottom_to_top=(True, False))
    def test_z_sort(self, info, bottom_to_top, make_poscar):
        """Check successful sorting of atoms by out-of-plane position."""
        slab, *_ = make_poscar(info)
        slab.sort_by_z(bottom_to_top=bottom_to_top)
        _ordered = operator.ge if bottom_to_top else operator.le                    # TODO: swap when flipping .cartpos
        assert all(_ordered(at1.cartpos[2], at2.cartpos[2])
                   for at1, at2 in pairwise(slab))

    @pytest.mark.skip(reason='to be implemented')
    def test_element_sort(self):
        """Check correct element-based sorting of a Slab."""








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

