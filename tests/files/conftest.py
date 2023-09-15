"""Test configuration for viperleed.tests.files.

Created on 2023-09-06

@author: Michele Riva (@michele-riva)

Fixtures
--------
poscar_with_group
    A Slab from POSCAR with known symmetry, an Rparams and a TestInfo.
    Parametrized with all known (non-bulk) POSCAR files.
tmp_poscar
    Fresh path to a 'POSCAR' file.
"""

from pathlib import Path
import sys

from pytest_cases import fixture, parametrize_with_cases

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Will be fixed in installable version
from viperleed.tleedmlib import symmetry

from ..poscar_slabs import CasePOSCARSlabs
# pylint: enable=wrong-import-position


@fixture
def tmp_poscar(tmp_path):
    """Return a fresh path to a POSCAR file."""
    return tmp_path / 'POSCAR'


@fixture
@parametrize_with_cases('args', cases=CasePOSCARSlabs)
def poscar_with_group(args):
    """Return slab, Rparams and TestInfo after symmetry-finding."""
    slab, rpars, info = args
    symmetry.findSymmetry(slab, rpars)
    return slab, rpars, info
