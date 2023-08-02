"""sections/conftest.py

Created on 2023-07-26

@author: Alexander M. Imre
"""

import pytest
import sys
import os
from pathlib import Path
from copy import deepcopy
import numpy as np

from viperleed.lib.files.poscar import readPOSCAR
from viperleed.lib.classes.rparams import Rparams


@pytest.fixture()
def fe3o4_bulk_slab(poscar_path):
    file_name = "POSCAR_Fe3O4_(001)_cod1010369"
    file_path = poscar_path(file_name)
    slab = readPOSCAR(str(file_path))
    param = Rparams()
    param.LAYER_CUTS = [0.1, 0.2, '<', 'dz(1.0)']
    param.N_BULK_LAYERS = 2
    param.SYMMETRY_EPS =0.3
    param.SYMMETRY_EPS_Z = 0.3
    param.BULK_REPEAT = np.array([-0.0, -4.19199991, 4.19199991])
    slab.fullUpdate(param)
    bulk_slab = slab.makeBulkSlab(param)
    return slab, bulk_slab, param


@pytest.fixture()
def fe3o4_thick_bulk_slab(fe3o4_bulk_slab):
    slab, thin_bulk, param = fe3o4_bulk_slab
    thick_bulk = thin_bulk.doubleBulkSlab()
    return slab, thick_bulk, param
