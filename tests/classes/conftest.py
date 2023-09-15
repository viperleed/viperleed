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
    slab = poscar.read(poscars_path / "POSCAR_STO(110)-4x1")
    atom = slab.atlist[0]
    el = atom.el
    atom.disp_geo[el] = [-0.2, 0.0, 0.2]
    atom.disp_vib[el] = [-0.1, 0.0, 0.1]
    atom.disp_occ[el] = [0.7, 0.8, 0.9, 1.0]
    return atom

