"""Test configuration for viperleed.tests.classes.

Created on 2023-07-26

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)

Fixtures
--------
manually_displaced_atom
    One atom with displacements assigned manually to its element.
manual_slab_1_atom_trigonal
    A simple trigonal Slab with a single Atom.
manual_slab_3_atoms
    An orthorhombic slab with three atoms.
"""

from pathlib import Path
import sys

from pytest_cases import fixture
import numpy as np

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.classes.slab import Slab
# pylint: enable=wrong-import-position


@fixture
def manually_displaced_atom(ag100):
    """Return one atom with displacements assigned to its element."""
    slab, *_ = ag100
    atom = slab.atlist[0]
    element = atom.el
    atom.disp_geo[element] = np.array([[-0.2, 0, 0], [0.0, 0, 0], [0.2, 0, 0]])
    atom.disp_vib[element] = [-0.1, 0.0, 0.1]
    atom.disp_occ[element] = [0.7, 0.8, 0.9, 1.0]
    return atom


# TODO: the two slabs below are probably better considered as CASES

@fixture
def manual_slab_1_atom_trigonal():
    """Return a simple trigonal slab with a single carbon atom."""
    slab = Slab()
    slab.ucell = np.array([[1, 0, 0],
                           [-2.3, 3, 0],
                           [1, 2, 3]]).T
    slab.atlist.append(Atom('C', np.array([0.2, 0.7, 0.1]), 1, slab))
    slab.full_update(Rparams())
    return slab


@fixture
def manual_slab_3_atoms():
    """Return an orthorhombic slab with three carbon atoms."""
    slab = Slab()
    slab.ucell = np.diag([3., 4., 5.])
    positions = (-0.25, 0, 0), (0, 0, 0), (0.25, 0, 0)
    slab.atlist.extend(Atom('C', np.array(pos), i+1, slab)
                       for i, pos in enumerate(positions))
    slab.full_update(Rparams())
    return slab
