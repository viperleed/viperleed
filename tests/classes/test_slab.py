"""classes/test_slab.py

Created on 2023-07-28

@author: Alexander M. Imre
"""

import pytest
import sys
import os
from pathlib import Path
from copy import deepcopy
import numpy as np

from viperleed.lib.classes.slab import SymPlane
from viperleed.lib.symmetry import findBulkSymmetry


class TestSlabTransforms:
    def test_mirror(self, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        mirrored_slab = deepcopy(slab)
        symplane = SymPlane(pos=(0,0),
                            dr=np.array([0,1]),
                            abt=slab.ucell.T[:2,:2])
        mirrored_slab.mirror(symplane)
        mirrored_slab.collapseCartesianCoordinates()
        assert all(at.isSameXY(mir_at.cartpos[:2])
                for at, mir_at in
                zip(slab.atlist, reversed(mirrored_slab.atlist)))

    def test_180_rotation(self, manual_slab_3_atoms):
        slab = manual_slab_3_atoms
        rotated_slab = deepcopy(slab)
        rotated_slab.rotateAtoms((0,0), order=2)
        rotated_slab.collapseCartesianCoordinates()
        assert all(at.isSameXY(mir_at.cartpos[:2])
                for at, mir_at in
                zip(slab.atlist, reversed(rotated_slab.atlist)))


# Slab Matrix operations
def test_rotation_on_trigonal_slab(manual_slab_1_atom_trigonal):
    rot_15 = np.array([[ 0.96592583, -0.25881905,  0.        ],
                       [ 0.25881905,  0.96592583,  0.        ],
                       [ 0.        ,  0.        ,  1.        ]])
    expected_cell = np.array([[ 0.44828774, -2.1906707 ,  1.        ],
                              [ 0.77645714,  2.89777748,  2.        ],
                              [ 0.        ,  0.        ,  3.        ]])
    expected_atom_cartpos = [0.63317754, 1.5903101]
    slab = manual_slab_1_atom_trigonal
    slab.apply_matrix_transformation(rot_15)
    assert np.allclose(slab.ucell.T, expected_cell)
    assert np.allclose(slab.atlist[0].cartpos[:2], expected_atom_cartpos)


@pytest.mark.parametrize('fixture', ('fe3o4_bulk_slab', 'fe3o4_thick_bulk_slab'))
def test_bulk_symmetry_thin(fixture, request):
    _, bulk, param = request.getfixturevalue(fixture)
    findBulkSymmetry(bulk, param)
    assert bulk.bulk_screws == [4]
    assert len(bulk.bulk_glides) == 2


class Test_restore_oristate:
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

    def test_save_restore_oristate_vib(self, ag100_slab_with_displacements_and_offsets):
            slab, param = ag100_slab_with_displacements_and_offsets
            slab_copy = deepcopy(slab)
            for at in slab.atlist:
                at.offset_vib['all'] = 0.1
                at.mergeDisp(at.el)

            slab.restoreOriState()
            for (at_rest, at_orig) in zip(slab.atlist, slab_copy.atlist):
                assert np.allclose(at_rest.disp_vib['all'], at_orig.disp_vib['all'])

    def test_save_restore_oristate_occ(self, ag100_slab_with_displacements_and_offsets):
            slab, param = ag100_slab_with_displacements_and_offsets
            slab_copy = deepcopy(slab)
            for at in slab.atlist:
                at.offset_occ[at.el] = 0.1
                at.mergeDisp(at.el)

            slab.restoreOriState()
            for (at_rest, at_orig) in zip(slab.atlist, slab_copy.atlist):
                assert np.allclose(at_rest.disp_occ[at_rest.el], at_orig.disp_occ[at_rest.el])