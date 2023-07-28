"""files/test_vibrocc.py

Created on 2023-07-28

@author: Alexander M. Imre
"""

import pytest
import numpy as np


def test_read_VIBROCC_offset_occ(ag100_slab_with_displacements_and_offsets):
    slab, param = ag100_slab_with_displacements_and_offsets
    assert np.allclose(slab.atlist[0].offset_occ['Ag'], -0.1)


def test_interpret_VIBROCC_offset_allowed(ag100_slab_with_displacements_and_offsets):
    slab, param = ag100_slab_with_displacements_and_offsets
    for atom in slab.atlist:
        atom.mergeDisp(atom.el)
    assert np.allclose(slab.atlist[0].disp_occ['Ag'], [0.4, 0.5, 0.6, 0.7, 0.8, 0.9])


def test_interpret_VIBROCC_offset_not_allowed(ag100_slab_with_displacements_and_offsets):
    slab, param = ag100_slab_with_displacements_and_offsets
    atom = slab.atlist[0]
    atom.offset_occ[atom.el] = +0.2
    atom.mergeDisp(atom.el)
    assert np.allclose(atom.disp_occ['Ag'], [0.5, 0.6, 0.7, 0.8, 0.9, 1.0])