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

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Cannot do anything about it until we make viperleed installable
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
        mirrored_slab.mirror(symplane)
        mirrored_slab.collapseCartesianCoordinates()
        assert all(
            at.isSameXY(mir_at.cartpos[:2])
            for at, mir_at in zip(slab, reversed(mirrored_slab))
            )

    def test_180_rotation(self, manual_slab_3_atoms):
        """Test the expected outcome of rotating atoms of a simple slab."""
        slab = manual_slab_3_atoms
        rotated_slab = deepcopy(slab)
        rotated_slab.rotateAtoms((0, 0), order=2)
        rotated_slab.collapseCartesianCoordinates()
        assert all(
            at.isSameXY(rot_at.cartpos[:2])
            for at, rot_at in zip(slab, reversed(rotated_slab))
            )


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


def test_add_one_atom_n_elements():
    """Check that adding one atom to a slab updates elements correctly."""
    slab = Slab()
    new_atom = Atom('C', (0, 0, 0), 1, slab)
    slab.atlist.append(new_atom)
    slab.update_element_count()
    assert new_atom.el in slab.elements
    assert slab.n_per_elem[new_atom.el] == 1
