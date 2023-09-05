"""classes/conftest.py

Created on 2023-07-26

@author: Alexander M. Imre
"""

import sys
import os
from pathlib import Path

import numpy as np
import pytest

vpr_path = str(Path(__file__).parent.parent.parent.parent)
if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

from viperleed.tleedmlib.classes import rparams
from viperleed.tleedmlib.files import poscar


@pytest.fixture(scope="function")
def atom_with_disp_and_offset(poscars_path):
    slab = poscar.read(poscars_path / "POSCAR_STO(100)-4x1")
    atom = slab.atlist[0]
    el = atom.el
    atom.disp_geo[el] = [-0.2, 0.0, 0.2]
    atom.disp_vib[el] = [-0.1, 0.0, 0.1]
    atom.disp_occ[el] = [0.7, 0.8, 0.9, 1.0]
    return atom

@pytest.fixture()
def fe3o4_bulk_slab(poscars_path):
    file_name = "POSCAR_Fe3O4_(001)_cod1010369"
    slab = poscar.read(poscars_path / file_name)
    param = rparams.Rparams()
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
