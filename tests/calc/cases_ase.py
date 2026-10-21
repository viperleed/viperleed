"""Test cases for test_viperleed_from_ase."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-09-15'
__license__ = 'GPLv3+'

from dataclasses import dataclass

import ase.build

from ..helpers import InfoBase


@dataclass(repr=False)
class ASEInfo(InfoBase):
    """Collection of test-related information for an ASE case."""
    n_atoms: int = 0


def case_ase_ni_100_1x1_cell():
    """Return an ase.Atoms Ni(100)-1x1 with 6 layers, and a its info."""
    element = 'Ni'
    _atoms = ase.build.fcc110(element, size=(1, 1, 6), vacuum=3)
    info = ASEInfo()
    info.n_atoms = 6
    return _atoms, info
