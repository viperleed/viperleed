"""files/test_displacements.py

Created on 2023-07-28

@author: Alexander M. Imre
"""

import pytest
import numpy as np


class Test_readDISPLACEMENTS:
    def test_read_DISPLACEMENTS_geo(self, ag100_slab_with_displacements_and_offsets):
        slab, param = ag100_slab_with_displacements_and_offsets
        assert np.allclose(np.array(slab.atlist[0].disp_geo['all'])[:, 2], [0.2, 0.1, 0.0, -0.1, -0.2])

    def test_read_DISPLACEMENTS_vib(self, ag100_slab_with_displacements_and_offsets):
        slab, param = ag100_slab_with_displacements_and_offsets
        assert np.allclose(slab.atlist[0].disp_vib['all'], [-0.1, 0.0, 0.1])

    def test_read_DISPLACEMENTS_occ(self, ag100_slab_with_displacements_and_offsets):
        slab, param = ag100_slab_with_displacements_and_offsets
        assert np.allclose(slab.atlist[0].disp_occ['Ag'], [0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
