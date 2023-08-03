"""classes/conftest.py

Created on 2023-07-26

@author: Alexander M. Imre
"""

import pytest
import sys
import os
from pathlib import Path

from viperleed.calc.files.poscar import readPOSCAR


@pytest.fixture(scope="function")
def atom_with_disp_and_offset(poscars_path):
    slab = readPOSCAR(poscars_path / "POSCAR_STO(100)-4x1")
    atom = slab.atlist[0]
    el = atom.el
    atom.disp_geo[el] = [-0.2, 0.0, 0.2]
    atom.disp_vib[el] = [-0.1, 0.0, 0.1]
    atom.disp_occ[el] = [0.7, 0.8, 0.9, 1.0]
    return atom
