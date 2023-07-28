"""test_symmetry.py

Created on 2023-07-28

@author: Alexander M. Imre
"""

import pytest
import sys
import os
from pathlib import Path
from copy import deepcopy
import numpy as np

vpr_path = str(Path(__file__).parent.parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))


from viperleed.tleedmlib.files.displacements import readDISPLACEMENTS, readDISPLACEMENTS_block
from viperleed.tleedmlib.files.poscar import readPOSCAR
from viperleed.tleedmlib.files.vibrocc import readVIBROCC
from viperleed.tleedmlib.symmetry import findSymmetry, enforceSymmetry
from viperleed.tleedmlib.psgen import runPhaseshiftGen_old
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.classes.slab import Slab

class TestSymmetry():
    def test_any_pg_found(self, slab_pg_rp):
        _, slab_pg, _ = slab_pg_rp
        assert slab_pg != 'unknown'

    def test_pg_correct(self, slab_and_expectations, slab_pg_rp):
        _, _, expected_pg, _ = slab_and_expectations
        _, slab_pg, _ = slab_pg_rp
        assert slab_pg == expected_pg

    @pytest.mark.parametrize("displacement", [(4, (np.array([0.2, 0, 0]),)),
                                            (4, (np.array([0, 0.2, 0]),)),
                                            (4, (np.array([0, 0, 0.2]),)),
                                            ])
    def test_preserve_symmetry_with_displacement(self, displacement, slab_and_expectations, slab_pg_rp):
        slab, _, expected_pg, offset_at = slab_and_expectations
        _, _, rp = slab_pg_rp
        sl_copy = deepcopy(slab)
        
        # manually assign displacements
        sl_copy.atlist[offset_at].assignDisp(*displacement)

        for at in sl_copy.atlist:
            disp = at.disp_geo_offset['all'][0]
            at.cartpos += disp
        sl_copy.getFractionalCoordinates()

        assert findSymmetry(sl_copy, rp) == expected_pg
