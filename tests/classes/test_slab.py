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
from viperleed.tleedmlib.classes.slab import SymPlane
# pylint: enable=wrong-import-position


class TestAtomTransforms:                                                       # TODO: add test for correct rotation and mirror on slanted cell. Probably enough to is_rotation_symmetric and is_mirror_symmetric with appropriate slabs. Could use "POSCAR_36C_cm" or the bulk of Fe3O4.
    """Test simple transformations of the atoms of a slab."""

    def test_mirror(self, manual_slab_3_atoms):
        """Test the expected outcome of mirroring atoms of a simple slab."""
        slab = manual_slab_3_atoms
        mirrored_slab = deepcopy(slab)
        symplane = SymPlane((0, 0), (0, 1), abt=slab.surface_vectors)
        mirrored_slab.mirror(symplane)
        mirrored_slab.collapseCartesianCoordinates()
        assert all(
            at.isSameXY(mir_at.cartpos[:2])
            for at, mir_at in zip(slab.atlist, reversed(mirrored_slab.atlist))
            )

    def test_180_rotation(self, manual_slab_3_atoms):
        """Test the expected outcome of rotating atoms of a simple slab."""
        slab = manual_slab_3_atoms
        rotated_slab = deepcopy(slab)
        rotated_slab.rotateAtoms((0, 0), order=2)
        rotated_slab.collapseCartesianCoordinates()
        assert all(
            at.isSameXY(rot_at.cartpos[:2])
            for at, rot_at in zip(slab.atlist, reversed(rotated_slab.atlist))
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


class Test_restore_oristate:
    @pytest.mark.xfail(reason="Awaiting new implementation of displacements.")
    def test_save_restore_oristate_geo(self, ag100_slab_with_displacements_and_offsets):
            slab, param = ag100_slab_with_displacements_and_offsets
            slab_copy = deepcopy(slab)
            for at in slab.atlist:
                at.disp_geo_offset['all'] = np.array([0.1, 0.0, 0.0])
                at.offset_geo['all'] = np.array([0.0, 0.0, 0.1])
                at.mergeDisp(at.el)

            slab.restoreOriState()
            for (at_rest, at_orig) in zip(slab.atlist, slab_copy.atlist):
                assert np.allclose(at_rest.disp_geo['all'], at_orig.disp_geo['all'])

    @pytest.mark.xfail(reason="Awaiting new implementation of displacements.")
    def test_save_restore_oristate_vib(self, ag100_slab_with_displacements_and_offsets):
            slab, param = ag100_slab_with_displacements_and_offsets
            slab_copy = deepcopy(slab)
            for at in slab.atlist:
                at.offset_vib['all'] = 0.1
                at.mergeDisp(at.el)

            slab.restoreOriState()
            for (at_rest, at_orig) in zip(slab.atlist, slab_copy.atlist):
                assert np.allclose(at_rest.disp_vib['all'], at_orig.disp_vib['all'])

    @pytest.mark.xfail(reason="Awaiting new implementation of displacements.")
    def test_save_restore_oristate_occ(self, ag100_slab_with_displacements_and_offsets):
            slab, param = ag100_slab_with_displacements_and_offsets
            slab_copy = deepcopy(slab)
            for at in slab.atlist:
                at.offset_occ[at.el] = 0.1
                at.mergeDisp(at.el)

            slab.restoreOriState()
            for (at_rest, at_orig) in zip(slab.atlist, slab_copy.atlist):
                assert np.allclose(at_rest.disp_occ[at_rest.el], at_orig.disp_occ[at_rest.el])