"""Test configuration for viperleed.tests.calc.files.

Fixtures
--------
displaced_atom
    An Atom with displacements assigned (but not applied).
poscar_with_group
    A Slab from POSCAR with known symmetry, an Rparams and a TestInfo.
    Parametrized with all known (non-bulk) POSCAR files.
tmp_poscar
    Fresh path to a 'POSCAR' file.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-09-06'
__license__ = 'GPLv3+'

from pytest_cases import fixture, parametrize_with_cases

from viperleed.calc import symmetry

from ..poscar_slabs import CasePOSCARSlabs


@fixture
def displaced_atom(ag100_with_displacements_and_offsets):
    """Return an Atom with displacements assigned (but not applied)."""
    slab, *_ = ag100_with_displacements_and_offsets
    return slab.atlist[0]


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
