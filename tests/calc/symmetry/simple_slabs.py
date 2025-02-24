"""Module simple_slabs of viperleed.tests.calc.symmetry.

This module defines two classes as containers of pytest cases. The
cases consist of hand-made slabs with all relevant combinations of
plane groups, including atoms at all Wyckoff positions.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-04-13'
__license__ = 'GPLv3+'

from enum import Enum

import numpy as np
from pytest_cases import case, parametrize

from viperleed.calc.classes.atom import Atom
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.slab import Slab

from ..tags import CaseTag as Tag
from ..testinfo import DisplacementInfo
from ..testinfo import TestInfo


def remove_atoms(slab, param, *atom_nrs):
    """Remove atoms with the given numbers from slab."""
    for atom_nr in atom_nrs:
        slab.atlist.remove(slab.atlist.get(atom_nr))
    slab.update_element_count()
    slab.full_update(param)


def tilt_c_axis(slab, direction):
    """Modify slab's c vector to tilt in a fractional direction."""
    slab.c_vector[:2] = 0.3 * np.dot(direction, slab.ab_cell.T)
    slab.ucell_ori = slab.ucell.copy()
    slab.update_cartesian_from_fractional()


def set_displacement_infos(test_info, *displaced_atoms):
    """Assign DisplacementInfos from a sequence.

    Parameters
    ----------
    test_info : TestInfo
        The object that will be updated.
    *displaced_atoms : Sequence
        Each item is a pair (int, bool) specifying which atom
        number (1-based) is to be displaced, and whether it is
        on a rotation axis.

    Returns
    -------
    None.
    """
    test_info.displacements.extend(DisplacementInfo(*at_info)
                                   for at_info in displaced_atoms)


class UnitCell(Enum):
    """An enumeration of unit cells."""

    OBLIQUE = (1, 0, 0), (0.2, 1, 0), (0, 0, 1)
    RECTANGULAR = (2., 0, 0), (0, 1, 0), (0, 0, 1)
    RHOMBIC = (1, 0.5, 0), (-1, 0.5, 0), (0, 0, 1)
    SQUARE = (1., 0, 0), (0, 1, 0), (0, 0, 1)
    HEXAGONAL = (1, 0, 0), (-0.5, (3**0.5)/2, 0), (0, 0, 1)


class CaseSimpleSlabs:  # pylint: disable=too-many-public-methods
    """A collection of hand-made slabs with select symmetry."""

    @staticmethod
    def add_atoms(slab, atoms):
        """Add atoms to a slab."""
        slab.atlist.extend(atoms)
        slab.update_element_count()

    @staticmethod
    def make_slab(ucell):
        """Return an empty slab with a given unit cell."""
        if isinstance(ucell, str):
            ucell = UnitCell[ucell.upper()]
        if isinstance(ucell, UnitCell):
            ucell = ucell.value
        slab = Slab()
        slab.ucell = 8*np.array(ucell, dtype=float).T
        slab.ucell_ori = slab.ucell.copy()
        return slab

    @staticmethod
    def update(slab, param=None):
        """Update a slab from param, and return the latter."""
        if not param:
            param = Rparams()
        slab.full_update(param)
        return param

    def case_p1(self):
        """Return a slab with p1 plane group."""
        slab = self.make_slab('oblique')
        x, y = 0.1, 0.15  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0, 0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set((), (), {1: {1}, 2: {2}}, 'p1')
        return slab, param, info

    def case_p2(self):
        """Return a p2-symmetric slab."""
        slab = self.make_slab('oblique')
        x, y = 0.1, 0.15  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (-x, -y, 0.1), 3, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set((), (1,), {1: {1}, 2: {2, 3}}, 'p2')
        set_displacement_infos(info, (2, False))
        return slab, param, info

    def case_pm_10(self):
        """Return a pm[1 0]-symmetric slab."""
        slab = self.make_slab('rectangular')
        x, y = 0.1, 0.15                 # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (x, -y, 0.1), 3, slab),
                              Atom('Ni', (2*x, 0.5, 0.3), 4, slab),))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set((1, 4), (), {1: {1}, 2: {2, 3}, 4: {4}}, 'pm')
        set_displacement_infos(info, (2, False))
        return slab, param, info

    def case_pg_10(self):
        """Return a pg[1 0]-symmetric slab."""
        slab = self.make_slab('rectangular')
        x, y = 0.1, 0.15  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Ca', (0.5, 0.0, 0.0), 2, slab),
                              Atom('Fe', (x, -y, 0.1), 3, slab),
                              Atom('Fe', (0.5+x, y, 0.1), 4, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set((), (), {1: {1, 2}, 3: {3, 4}}, 'pg')
        set_displacement_infos(info, (3, False))
        return slab, param, info

    def case_cm_1m1(self):
        """Return a cm[1 -1]-symmetric slab."""
        slab = self.make_slab('rhombic')
        x, y = 0.1, -0.2  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.5, 0.5, 0.0), 1, slab),
                              Atom('Fe', (0.5+x, 0.5+y, 0.1), 2, slab),
                              Atom('Fe', (0.5-y, 0.5-x, 0.1), 3, slab),
                              Atom('Ni', (0.5-x, 0.5+x, 0.2), 4, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set((1, 4), (), {1: {1}, 2: {2, 3}, 4: {4}}, 'cm')
        set_displacement_infos(info, (2, False))
        return slab, param, info

    @parametrize(ucell=('rectangular', 'square'))
    def case_rcm(self, ucell):
        """Return an rcm[1 0]-symmetric slab."""
        slab = self.make_slab(ucell)
        x, y = 0.1, 0.15  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Ca', (0.5, 0.5, 0.0), 2, slab),
                              Atom('Fe', (x, y, 0.1), 3, slab),
                              Atom('Fe', (x, -y, 0.1), 4, slab),
                              Atom('Fe', (0.5+x, 0.5+y, 0.1), 5, slab),
                              Atom('Fe', (0.5+x, 0.5-y, 0.1), 6, slab),
                              Atom('Ni', (3*x, 0.0, 0.2), 7, slab),
                              Atom('Ni', (0.5+3*x, 0.5, 0.2), 8, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set((1, 2, 7, 8), (),
                          {1: {1, 2}, 3: {3, 4, 5, 6}, 7: {7, 8}}, 'rcm')
        set_displacement_infos(info, (3, False))
        return slab, param, info

    @parametrize(ucell=('rectangular', 'square'))
    def case_pmm(self, ucell):
        """Return a pmm-symmetric slab."""
        slab = self.make_slab(ucell)
        x, y = 0.1, 0.15  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (x, -y, 0.1), 3, slab),
                              Atom('Fe', (-x, y, 0.1), 4, slab),
                              Atom('Fe', (-x, -y, 0.1), 5, slab),
                              Atom('Ni', (0.5, 0.5, 0.2), 6, slab),
                              Atom('Ti', (0.5, 0.0, 0.3), 7, slab),
                              Atom('Cr', (0.0, 0.5, 0.4), 8, slab),
                              Atom('Mn', (x, 0.0, 0.5), 9, slab),
                              Atom('Mn', (-x, 0.0, 0.5), 10, slab),
                              Atom('Co', (0.0, y, 0.6), 11, slab),
                              Atom('Co', (0.0, -y, 0.6), 12, slab),
                              Atom('Cu', (0.5, y, 0.7), 13, slab),
                              Atom('Cu', (0.5, -y, 0.7), 14, slab),
                              Atom('Zn', (x, 0.5, 0.8), 15, slab),
                              Atom('Zn', (-x, 0.5, 0.8), 16, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            (9, 10, 11, 12, 13, 14, 15, 16),
            (1, 6, 7, 8),
            {1: {1}, 2: {2, 3, 4, 5}, 6: {6}, 7: {7}, 8: {8}, 9: {9, 10},
             11: {11, 12}, 13: {13, 14}, 15: {15, 16}},
            'pmm'
            )
        set_displacement_infos(info, (2, False))
        return slab, param, info

    @parametrize(ucell=('rectangular', 'square'))
    def case_pmg(self, ucell):
        """Return a pmg[1 0]-symmetric slab."""
        slab = self.make_slab(ucell)
        x, y = 0.1, -0.15  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Ca', (0.5, 0.0, 0.0), 2, slab),
                              Atom('Fe', (0.5+x, -y, 0.1), 3, slab),
                              Atom('Fe', (x, y, 0.1), 4, slab),
                              Atom('Fe', (-x, -y, 0.1), 5, slab),
                              Atom('Fe', (0.5-x, y, 0.1), 6, slab),
                              Atom('Ni', (0.5, 0.5, 0.2), 7, slab),
                              Atom('Ni', (0.0, 0.5, 0.2), 8, slab),
                              Atom('Ti', (0.25, 0.5-y, 0.4), 9, slab),
                              Atom('Ti', (0.75, 0.5+y, 0.4), 10, slab),
                              Atom('Mn', (0.25, x, 0.3), 11, slab),
                              Atom('Mn', (0.75, -x, 0.3), 12, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            (9, 10, 11, 12),
            (1, 2, 7, 8),
            {1: {1, 2}, 3: {3, 4, 5, 6}, 7: {7, 8}, 9: {9, 10}, 11: {11, 12}},
            'pmg'
            )
        set_displacement_infos(info, (3, False))
        return slab, param, info

    @parametrize(ucell=('rectangular', 'square'))
    def case_pgg(self, ucell):
        """Return a pgg-symmetric slab."""
        slab = self.make_slab(ucell)
        x, y = 0.1, 0.15  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Ca', (0.5, 0.5, 0.0), 2, slab),
                              Atom('Fe', (x, y, 0.1), 4, slab),
                              Atom('Fe', (-x, -y, 0.1), 5, slab),
                              Atom('Fe', (0.5-x, 0.5+y, 0.1), 6, slab),
                              Atom('Fe', (0.5+x, 0.5-y, 0.1), 3, slab),
                              Atom('Ni', (0.0, 0.5, 0.2), 7, slab),
                              Atom('Ni', (0.5, 0.0, 0.2), 8, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set((), (1, 2, 7, 8),
                          {1: {1, 2}, 3: {3, 4, 5, 6}, 7: {7, 8}}, 'pgg')
        set_displacement_infos(info, (4, False))
        return slab, param, info

    def case_cmm(self):
        """Return a cmm-symmetric slab."""
        slab = self.make_slab('rhombic')
        x, y = 0.1, -0.2  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.5, 0.5, 0.0), 1, slab),
                              Atom('Fe', (0.5+x, 0.5+y, 0.1), 2, slab),
                              Atom('Fe', (0.5+y, 0.5+x, 0.1), 3, slab),
                              Atom('Fe', (0.5-y, 0.5-x, 0.1), 4, slab),
                              Atom('Fe', (0.5-x, 0.5-y, 0.1), 5, slab),
                              Atom('Ti', (0.0, 0.0, 0.2), 6, slab),
                              Atom('Ni', (0.5-x, 0.5+x, 0.3), 7, slab),
                              Atom('Ni', (0.5+x, 0.5-x, 0.3), 8, slab),
                              Atom('Cu', (0.5+x, 0.5+x, 0.4), 9, slab),
                              Atom('Cu', (0.5-x, 0.5-x, 0.4), 10, slab),
                              Atom('Mn', (0.0, 0.5, 0.5), 11, slab),
                              Atom('Mn', (0.5, 0.0, 0.5), 12, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            (7, 8, 9, 10),
            (1, 6, 11, 12),
            {1: {1}, 2: {2, 3, 4, 5}, 6: {6}, 7: {7, 8}, 9: {9, 10},
             11: {11, 12}},
            'cmm'
            )
        set_displacement_infos(info, (2, False))
        return slab, param, info

    @parametrize(ucell=('rectangular', 'square'))
    def case_rcmm(self, ucell):
        """Return an rcmm-symmetric slab."""
        slab = self.make_slab(ucell)
        x, y = 0.1, 0.15  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Ca', (0.5, 0.5, 0.0), 2, slab),
                              Atom('Fe', (x, y, 0.1), 3, slab),
                              Atom('Fe', (x, -y, 0.1), 4, slab),
                              Atom('Fe', (-x, y, 0.1), 5, slab),
                              Atom('Fe', (-x, -y, 0.1), 6, slab),
                              Atom('Fe', (0.5+x, 0.5+y, 0.1), 7, slab),
                              Atom('Fe', (0.5+x, 0.5-y, 0.1), 8, slab),
                              Atom('Fe', (0.5-x, 0.5+y, 0.1), 9, slab),
                              Atom('Fe', (0.5-x, 0.5-y, 0.1), 10, slab),
                              Atom('Ni', (x, 0.0, 0.2), 11, slab),
                              Atom('Ni', (-x, 0.0, 0.2), 12, slab),
                              Atom('Ni', (0.5+x, 0.5, 0.2), 13, slab),
                              Atom('Ni', (0.5-x, 0.5, 0.2), 14, slab),
                              Atom('Cu', (0.0, y, 0.3), 15, slab),
                              Atom('Cu', (0.0, -y, 0.3), 16, slab),
                              Atom('Cu', (0.5, 0.5+y, 0.3), 17, slab),
                              Atom('Cu', (0.5, 0.5-y, 0.3), 18, slab),
                              Atom('Mn', (0.25, 0.25, 0.4), 19, slab),
                              Atom('Mn', (0.25, 0.75, 0.4), 20, slab),
                              Atom('Mn', (0.75, 0.25, 0.4), 21, slab),
                              Atom('Mn', (0.75, 0.75, 0.4), 22, slab),
                              Atom('Ti', (0.5, 0.0, 0.5), 23, slab),
                              Atom('Ti', (0.0, 0.5, 0.5), 24, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            (11, 12, 13, 14, 15, 16, 17, 18),
            (1, 2, 19, 20, 21, 22, 23, 24),
            {1: {1, 2}, 3: {3, 4, 5, 6, 7, 8, 9, 10}, 11: {11, 12, 13, 14},
             15: {15, 16, 17, 18}, 19: {19, 20, 21, 22}, 23: {23, 24}},
            'rcmm'
            )
        set_displacement_infos(info, (3, False))
        return slab, param, info

    def case_p4(self):
        """Return a p4-symmetric slab."""
        slab = self.make_slab('square')
        x, y = 0.1, 0.2  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (y, -x, 0.1), 3, slab),
                              Atom('Fe', (-y, x, 0.1), 4, slab),
                              Atom('Fe', (-x, -y, 0.1), 5, slab),
                              Atom('Ni', (0.5, 0.5, 0.2), 6, slab),
                              Atom('Cu', (0.5, 0.0, 0.3), 7, slab),
                              Atom('Cu', (0.0, 0.5, 0.3), 8, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set((), (1, 6, 7, 8),
                          {1: {1}, 2: {2, 3, 4, 5}, 6: {6}, 7: {7, 8}}, 'p4')
        set_displacement_infos(info, (2, False))
        return slab, param, info

    def case_p4m(self):
        """Return a p4m-symmetric slab."""
        slab = self.make_slab('square')
        x, y = 0.1, 0.2  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (y, x, 0.1), 3, slab),
                              Atom('Fe', (y, -x, 0.1), 4, slab),
                              Atom('Fe', (-y, x, 0.1), 5, slab),
                              Atom('Fe', (x, -y, 0.1), 6, slab),
                              Atom('Fe', (-x, y, 0.1), 7, slab),
                              Atom('Fe', (-x, -y, 0.1), 8, slab),
                              Atom('Fe', (-y, -x, 0.1), 9, slab),
                              Atom('Ti', (x, 0.0, 0.2), 10, slab),
                              Atom('Ti', (0.0, x, 0.2), 11, slab),
                              Atom('Ti', (-x, 0.0, 0.2), 12, slab),
                              Atom('Ti', (0.0, -x, 0.2), 13, slab),
                              Atom('Co', (y, y, 0.3), 14, slab),
                              Atom('Co', (y, -y, 0.3), 15, slab),
                              Atom('Co', (-y, y, 0.3), 16, slab),
                              Atom('Co', (-y, -y, 0.3), 17, slab),
                              Atom('Zn', (0.5, x, 0.4), 18, slab),
                              Atom('Zn', (x, 0.5, 0.4), 19, slab),
                              Atom('Zn', (0.5, -x, 0.4), 20, slab),
                              Atom('Zn', (-x, 0.5, 0.4), 21, slab),
                              Atom('Ni', (0.5, 0.5, 0.5), 22, slab),
                              Atom('Cu', (0.5, 0.0, 0.6), 23, slab),
                              Atom('Cu', (0.0, 0.5, 0.6), 24, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            tuple(range(10, 22)),
            (1, 22, 23, 24),
            {1: {1}, 2: {2, 3, 4, 5, 6, 7, 8, 9}, 10: {10, 11, 12, 13},
             14: {14, 15, 16, 17}, 18: {18, 19, 20, 21}, 22: {22},
             23: {23, 24}},
            'p4m'
            )
        set_displacement_infos(info, (2, False))
        return slab, param, info

    def case_p4g(self):
        """Return a p4g-symmetric slab."""
        slab = self.make_slab('square')
        x, y = 0.1, 0.2  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (y, -x, 0.1), 3, slab),
                              Atom('Fe', (-y, x, 0.1), 4, slab),
                              Atom('Fe', (-x, -y, 0.1), 5, slab),
                              Atom('Fe', (0.5+y, 0.5+x, 0.1), 6, slab),
                              Atom('Fe', (0.5+x, 0.5-y, 0.1), 7, slab),
                              Atom('Fe', (0.5-y, 0.5-x, 0.1), 8, slab),
                              Atom('Fe', (0.5-x, 0.5+y, 0.1), 9, slab),
                              Atom('Ti', (x, 0.5-x, 0.2), 10, slab),
                              Atom('Ti', (0.5-x, -x, 0.2), 11, slab),
                              Atom('Ti', (0.5+x, x, 0.2), 12, slab),
                              Atom('Ti', (-x, 0.5+x, 0.2), 13, slab),
                              Atom('Ca', (0.5, 0.5, 0.0), 14, slab),
                              Atom('Cu', (0.5, 0.0, 0.6), 15, slab),
                              Atom('Cu', (0.0, 0.5, 0.6), 16, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            (10, 11, 12, 13),
            (1, 14, 15, 16),
            {1: {1, 14}, 2: {2, 3, 4, 5, 6, 7, 8, 9},
             10: {10, 11, 12, 13}, 15: {15, 16}},
            'p4g'
            )
        set_displacement_infos(info, (2, False))
        return slab, param, info

    def case_p3(self):
        """Return a p3-symmetric slab."""
        slab = self.make_slab('hexagonal')
        x, y = 0.1, 0.3  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (-y, x-y, 0.1), 3, slab),
                              Atom('Fe', (y-x, -x, 0.1), 4, slab),
                              Atom('Ti', (2/3, 1/3, 0.2), 5, slab),
                              Atom('Cu', (1/3, 2/3, 0.3), 6, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set((), (1, 5, 6),
                          {1: {1}, 2: {2, 3, 4}, 5: {5}, 6: {6}}, 'p3')
        set_displacement_infos(info, (2, False))
        return slab, param, info

    def case_p3m1(self):
        """Return a p3m1-symmetric slab."""
        slab = self.make_slab('hexagonal')
        x, y = 0.15, 0.45  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (-y, x-y, 0.1), 3, slab),
                              Atom('Fe', (y-x, -x, 0.1), 4, slab),
                              Atom('Fe', (-y, -x, 0.1), 5, slab),
                              Atom('Fe', (y-x, y, 0.1), 6, slab),
                              Atom('Fe', (x, x-y, 0.1), 7, slab),
                              Atom('Ti', (2/3, 1/3, 0.2), 8, slab),
                              Atom('Cu', (1/3, 2/3, 0.3), 9, slab),
                              Atom('Mn', (x, -x, 0.4), 10, slab),
                              Atom('Mn', (x, 2*x, 0.4), 11, slab),
                              Atom('Mn', (-2*x, -x, 0.4), 12, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            (10, 11, 12),
            (1, 8, 9),
            {1: {1}, 2: {2, 3, 4, 5, 6, 7}, 8: {8}, 9: {9}, 10: {10, 11, 12}},
            'p3m1'
            )
        set_displacement_infos(info, (2, False))
        return slab, param, info

    def case_p31m(self):
        """Return a p31m-symmetric slab."""
        slab = self.make_slab('hexagonal')
        x, y = 0.15, 0.45  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (-y, x-y, 0.1), 3, slab),
                              Atom('Fe', (y-x, -x, 0.1), 4, slab),
                              Atom('Fe', (y, x, 0.1), 5, slab),
                              Atom('Fe', (x-y, -y, 0.1), 6, slab),
                              Atom('Fe', (-x, y-x, 0.1), 7, slab),
                              Atom('Ti', (2/3, 1/3, 0.2), 8, slab),
                              Atom('Ti', (1/3, 2/3, 0.2), 9, slab),
                              Atom('Mn', (x, 0, 0.4), 10, slab),
                              Atom('Mn', (0, x, 0.4), 11, slab),
                              Atom('Mn', (-x, -x, 0.4), 12, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            (10, 11, 12),
            (1, 8, 9),
            {1: {1}, 2: {2, 3, 4, 5, 6, 7}, 8: {8, 9}, 10: {10, 11, 12}},
            'p31m'
            )
        set_displacement_infos(info, (2, False))
        return slab, param, info

    def case_p6(self):
        """Return a p6-symmetric slab."""
        slab = self.make_slab('hexagonal')
        x, y = 0.15, 0.45  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (-y, x-y, 0.1), 3, slab),
                              Atom('Fe', (y-x, -x, 0.1), 4, slab),
                              Atom('Fe', (-x, -y, 0.1), 5, slab),
                              Atom('Fe', (y, y-x, 0.1), 6, slab),
                              Atom('Fe', (x-y, x, 0.1), 7, slab),
                              Atom('Ti', (2/3, 1/3, 0.2), 8, slab),
                              Atom('Ti', (1/3, 2/3, 0.2), 9, slab),
                              Atom('Cu', (0.5, 0, 0.3), 10, slab),
                              Atom('Cu', (0, 0.5, 0.3), 11, slab),
                              Atom('Cu', (0.5, 0.5, 0.3), 12, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            (),
            (1, 8, 9, 10, 11, 12),
            {1: {1}, 2: {2, 3, 4, 5, 6, 7}, 8: {8, 9}, 10: {10, 11, 12}},
            'p6'
            )
        set_displacement_infos(info, (2, False))
        return slab, param, info

    def case_p6m(self):
        """Return a p6m-symmetric slab."""
        slab = self.make_slab('hexagonal')
        x, y = 0.15, 0.45  # pylint: disable=invalid-name
        self.add_atoms(slab, (Atom('Ca', (0.0, 0.0, 0.0), 1, slab),
                              Atom('Fe', (x, y, 0.1), 2, slab),
                              Atom('Fe', (-y, x-y, 0.1), 3, slab),
                              Atom('Fe', (y-x, -x, 0.1), 4, slab),
                              Atom('Fe', (-y, -x, 0.1), 5, slab),
                              Atom('Fe', (y-x, y, 0.1), 6, slab),
                              Atom('Fe', (x, x-y, 0.1), 7, slab),
                              Atom('Fe', (y, x, 0.1), 8, slab),
                              Atom('Fe', (x-y, -y, 0.1), 9, slab),
                              Atom('Fe', (-x, y-x, 0.1), 10, slab),
                              Atom('Fe', (-x, -y, 0.1), 11, slab),
                              Atom('Fe', (y, y-x, 0.1), 12, slab),
                              Atom('Fe', (x-y, x, 0.1), 13, slab),
                              Atom('Ti', (2/3, 1/3, 0.2), 14, slab),
                              Atom('Ti', (1/3, 2/3, 0.2), 15, slab),
                              Atom('Cu', (0.5, 0, 0.3), 16, slab),
                              Atom('Cu', (0, 0.5, 0.3), 17, slab),
                              Atom('Cu', (0.5, 0.5, 0.3), 18, slab),
                              Atom('Mn', (x, -x, 0.4), 19, slab),
                              Atom('Mn', (x, 2*x, 0.4), 20, slab),
                              Atom('Mn', (2*x, x, 0.4), 21, slab),
                              Atom('Mn', (-2*x, -x, 0.4), 22, slab),
                              Atom('Mn', (-x, x, 0.4), 23, slab),
                              Atom('Mn', (-x, -2*x, 0.4), 24, slab),
                              Atom('Co', (y, 0.0, 0.5), 25, slab),
                              Atom('Co', (0.0, y, 0.5), 26, slab),
                              Atom('Co', (y, y, 0.5), 27, slab),
                              Atom('Co', (-y, 0.0, 0.5), 28, slab),
                              Atom('Co', (0.0, -y, 0.5), 29, slab),
                              Atom('Co', (-y, -y, 0.5), 30, slab)))
        param = self.update(slab)
        info = TestInfo()
        info.symmetry.set(
            tuple(range(19, 31)),
            (1, 14, 15, 16, 17, 18),
            {1: {1}, 2: {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
             14: {14, 15}, 16: {16, 17, 18}, 19: {19, 20, 21, 22, 23, 24},
             25: {25, 26, 27, 28, 29, 30}},
            'p6m'
            )
        set_displacement_infos(info, (2, False))
        return slab, param, info


class CaseSimpleHexagonalSlabs:
    """Collection of simple hexagonal slabs without the full p6m symmetry."""

    @staticmethod
    def make_p6m():
        """Return a p6m-symmetric slab."""
        return CaseSimpleSlabs().case_p6m()

    # The cm slabs are obtained by tilting the c axis in the right direction
    def case_hex_cm_11(self):
        """Return an hexagonal slab with cm[1 1] plane group."""
        slab, param, info = self.make_p6m()
        tilt_c_axis(slab, (1, 1))
        info.symmetry.hermann = 'cm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 18, 27, 30)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 8}, 3: {3, 9}, 4: {4, 10}, 5: {5, 11}, 6: {6, 12},
            7: {7, 13}, 14: {14, 15}, 16: {16, 17}, 18: {18}, 19: {19, 23},
            20: {20, 21}, 22: {22, 24}, 25: {25, 26}, 27: {27}, 28: {28, 29},
            30: {30}
            }
        return slab, param, info

    def case_hex_cm_1m1(self):
        """Return an hexagonal slab with cm[1 -1] plane group."""
        slab, param, info = self.make_p6m()
        tilt_c_axis(slab, (1, -1))
        info.symmetry.hermann = 'cm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 14, 15, 18, 19, 23)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 5}, 3: {3, 6}, 4: {4, 7}, 8: {8, 11}, 9: {9, 12},
            10: {10, 13}, 14: {14}, 15: {15}, 16: {16, 17}, 18: {18}, 19: {19},
            20: {20, 22}, 21: {21, 24}, 23: {23}, 25: {25, 29}, 26: {26, 28},
            27: {27, 30}
            }
        return slab, param, info

    @case(tags=(Tag.NEED_ROTATION,))
    def case_hex_cm_10(self):
        """Return an hexagonal slab with cm[1 0] plane group."""
        slab, param, info = self.make_p6m()
        tilt_c_axis(slab, (1, 0))
        info.symmetry.hermann = 'cm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 16, 25, 28)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 9}, 3: {3, 10}, 4: {4, 8}, 5: {5, 13}, 6: {6, 11},
            7: {7, 12}, 14: {14, 15}, 16: {16}, 17: {17, 18}, 19: {19, 21},
            20: {20, 24}, 22: {22, 23}, 25: {25}, 26: {26, 30}, 27: {27, 29},
            28: {28}
            }
        # Rotation: how far to rotate one of the current unit
        # vectors to get the one after the plane group has been
        # detected, or reduced from p6m/p31m
        info.symmetry.rotation = -60  # becomes [1 1]
        return slab, param, info

    @case(tags=(Tag.NEED_ROTATION,))
    def case_hex_cm_01(self):
        """Return an hexagonal slab with cm[0 1] plane group."""
        slab, param, info = self.make_p6m()
        tilt_c_axis(slab, (0, 1))
        info.symmetry.hermann = 'cm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 17, 26, 29)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 10}, 3: {3, 8}, 4: {4, 9}, 5: {5, 12}, 6: {6, 13},
            7: {7, 11}, 14: {14, 15}, 16: {16, 18}, 17: {17}, 19: {19, 24},
            20: {20, 23}, 21: {21, 22}, 25: {25, 30}, 26: {26}, 27: {27, 28},
            29: {29}
            }
        # Rotation: how far to rotate one of the current unit
        # vectors to get the one after the plane group has been
        # detected, or reduced from p6m/p31m
        info.symmetry.rotation = 60  # becomes [1 1]
        return slab, param, info

    @case(tags=(Tag.NEED_ROTATION,))
    def case_hex_cm_21(self):
        """Return an hexagonal slab with cm[2 1] plane group."""
        slab, param, info = self.make_p6m()
        tilt_c_axis(slab, (2, 1))
        info.symmetry.hermann = 'cm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 14, 15, 17, 21, 22)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 7}, 3: {3, 5}, 4: {4, 6}, 8: {8, 12}, 9: {9, 13},
            10: {10, 11}, 14: {14}, 15: {15}, 16: {16, 18}, 17: {17},
            19: {19, 20}, 21: {21}, 22: {22}, 23: {23, 24}, 25: {25, 27},
            26: {26, 29}, 28: {28, 30}
            }
        # Rotation: how far to rotate one of the current unit
        # vectors to get the one after the plane group has been
        # detected, or reduced from p6m/p3m1
        info.symmetry.rotation = 60  # becomes [1 -1]
        return slab, param, info

    @case(tags=(Tag.NEED_ROTATION,))
    def case_hex_cm_12(self):
        """Return an hexagonal slab with cm[1 2] plane group."""
        slab, param, info = self.make_p6m()
        tilt_c_axis(slab, (1, 2))
        info.symmetry.hermann = 'cm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 14, 15, 16, 20, 24)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 6}, 3: {3, 7}, 4: {4, 5}, 8: {8, 13}, 9: {9, 11},
            10: {10, 12}, 14: {14}, 15: {15}, 16: {16}, 17: {17, 18},
            19: {19, 22}, 20: {20}, 21: {21, 23}, 24: {24}, 25: {25, 28},
            26: {26, 27}, 29: {29, 30}
            }
        # Rotation: how far to rotate one of the current unit
        # vectors to get the one after the plane group has been
        # detected, or reduced from p6m/p3m1
        info.symmetry.rotation = -60  # becomes [1 -1]
        return slab, param, info

    # The cmm slabs, are instead obtained by removing a pair of atoms
    # across the sixfold axis. This removed the 3-fold rotation.
    # Which pair determines which mirror planes are left
    def case_hex_cmm_11(self):  # Also [1 -1]
        """Return an hexagonal slab with cmm[1 1] plane group."""
        slab, param, info = self.make_p6m()
        remove_atoms(slab, param, 19, 23)
        info.symmetry.hermann = 'cmm'
        info.symmetry.on_axes = (1, 16, 17, 18)
        info.symmetry.on_planes = (14, 15, 27, 30)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 5, 8, 11}, 3: {3, 6, 9, 12}, 4: {4, 7, 10, 13},
            14: {14, 15}, 16: {16, 17}, 18: {18}, 20: {20, 21, 22, 24},
            25: {25, 26, 28, 29}, 27: {27, 30}
            }
        return slab, param, info

    @case(tags=(Tag.NEED_ROTATION,))
    def case_hex_cmm_10(self):  # Also [1 2]
        """Return an hexagonal slab with cmm[1 0] plane group."""
        slab, param, info = self.make_p6m()
        remove_atoms(slab, param, 20, 24)
        info.symmetry.hermann = 'cmm'
        info.symmetry.on_axes = (1, 16, 17, 18)
        info.symmetry.on_planes = (14, 15, 25, 28)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 6, 9, 11}, 3: {3, 7, 10, 12}, 4: {4, 5, 8, 13},
            14: {14, 15}, 16: {16}, 17: {17, 18}, 19: {19, 21, 22, 23},
            25: {25, 28}, 26: {26, 27, 29, 30}
            }
        info.symmetry.invalid_direction = '[1 1]'
        # Rotation: how far to rotate one of the current unit
        # vectors to get the one after the plane group has been
        # detected, or reduced from p6m
        info.symmetry.rotation = -60
        return slab, param, info

    @case(tags=(Tag.NEED_ROTATION,))
    def case_hex_cmm_01(self):  # Also [2 1]
        """Return an hexagonal slab with cmm[0 1] plane group."""
        slab, param, info = self.make_p6m()
        remove_atoms(slab, param, 21, 22)
        info.symmetry.hermann = 'cmm'
        info.symmetry.on_axes = (1, 16, 17, 18)
        info.symmetry.on_planes = (14, 15, 26, 29)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 7, 10, 11}, 3: {3, 5, 8, 12}, 4: {4, 6, 9, 13},
            14: {14, 15}, 16: {16, 18}, 17: {17}, 19: {19, 20, 23, 24},
            25: {25, 27, 28, 30}, 26: {26, 29}
            }
        info.symmetry.invalid_direction = '[1 1]'
        # Rotation: how far to rotate one of the current unit
        # vectors to get the one after the plane group has been
        # detected, or reduced from p6m
        info.symmetry.rotation = 60
        return slab, param, info

    def case_hex_p2(self):
        """Return an hexagonal slab with p2 plane group."""
        # For this one, we remove another two atoms to get
        # rid of all the mirror planes
        slab, param, info = self.make_p6m()
        remove_atoms(slab, param, 6, 9, 19, 23)
        info.symmetry.hermann = 'p2'
        info.symmetry.on_axes = (1, 16, 17, 18)
        info.symmetry.on_planes = ()
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 11}, 3: {3, 12}, 5: {5, 8}, 4: {4, 13}, 7: {7, 10},
            14: {14, 15}, 16: {16}, 17: {17}, 18: {18}, 20: {20, 24},
            21: {21, 22}, 25: {25, 28}, 26: {26, 29}, 27: {27, 30}
            }
        return slab, param, info


class CaseSimpleSquareSlabs:
    """Collection of simple square slabs without the full p4m/g symmetry."""

    @staticmethod
    def make_p4g():
        """Return a p4g-symmetric slab."""
        return CaseSimpleSlabs().case_p4g()

    @staticmethod
    def make_p4m():
        """Return a p4m-symmetric slab."""
        return CaseSimpleSlabs().case_p4m()

    # As for the hex slabs, we tilt c on a p4g to get pg
    def case_square_pg_10(self):
        """Return a square slab with pg[1 0] plane group."""
        slab, param, info = self.make_p4g()
        tilt_c_axis(slab, (1, 0))
        info.symmetry.hermann = 'pg'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = ()
        info.symmetry.link_groups = {
            1: {1, 14}, 2: {2, 7}, 3: {3, 6}, 4: {4, 8}, 5: {5, 9},
            10: {10, 12}, 11: {11, 13}, 15: {15, 16}
            }
        return slab, param, info

    def case_square_pg_01(self):
        """Return a square slab with pg[0 1] plane group."""
        slab, param, info = self.make_p4g()
        tilt_c_axis(slab, (0, 1))
        info.symmetry.hermann = 'pg'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = ()
        info.symmetry.link_groups = {
            1: {1, 14}, 2: {2, 9}, 3: {3, 8}, 4: {4, 6}, 5: {5, 7},
            10: {10, 11}, 12: {12, 13}, 15: {15, 16}
            }
        return slab, param, info

    def case_square_pmm(self):
        """Return a square slab with pmm plane group."""
        slab, param, info = self.make_p4m()
        remove_atoms(slab, param, 11, 13)
        info.symmetry.hermann = 'pmm'
        info.symmetry.on_axes = (1, 22, 23, 24)
        info.symmetry.on_planes = (10, 12, 18, 19, 20, 21)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 6, 7, 8}, 3: {3, 4, 5, 9}, 10: {10, 12},
            14: {14, 15, 16, 17}, 18: {18, 20}, 19: {19, 21}, 22: {22},
            23: {23}, 24: {24}
            }
        return slab, param, info

    def case_square_pm_10(self):
        """Return a square slab with pm[1 0] plane group."""
        slab, param, info = self.case_square_pmm()
        tilt_c_axis(slab, (1, 0))
        info.symmetry.hermann = 'pm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 10, 12, 19, 21, 22, 23, 24)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 6}, 3: {3, 4}, 5: {5, 9}, 7: {7, 8}, 10: {10},
            12: {12}, 14: {14, 15}, 16: {16, 17}, 18: {18, 20}, 19: {19},
            21: {21}, 22: {22}, 23: {23}, 24: {24}
            }
        return slab, param, info

    def case_square_pm_01(self):
        """Return a square slab with pm[0 1] plane group."""
        slab, param, info = self.case_square_pmm()
        tilt_c_axis(slab, (0, 1))
        info.symmetry.hermann = 'pm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 18, 20, 22, 23, 24)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 7}, 3: {3, 5}, 4: {4, 9}, 6: {6, 8}, 10: {10, 12},
            14: {14, 16}, 15: {15, 17}, 18: {18}, 19: {19, 21}, 20: {20},
            22: {22}, 23: {23}, 24: {24}
            }
        return slab, param, info

    def case_square_cmm(self):
        """Return a square slab with cmm plane group."""
        # We could do the cmm also from p4g, but the origin would
        # need to be shifted to get the mirrors in the right spot
        slab, param, info = self.make_p4m()
        remove_atoms(slab, param, 4, 5, 6, 7)
        info.symmetry.hermann = 'cmm'
        info.symmetry.on_axes = (1, 22, 23, 24)
        info.symmetry.on_planes = (14, 15, 16, 17)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 3, 8, 9}, 10: {10, 11, 12, 13}, 14: {14, 17},
            15: {15, 16}, 18: {18, 19, 20, 21}, 22: {22}, 23: {23, 24}
            }
        return slab, param, info

    def case_square_cm_11(self):
        """Return a square slab with cm[1 1] plane group."""
        slab, param, info = self.case_square_cmm()
        tilt_c_axis(slab, (1, 1))
        info.symmetry.hermann = 'cm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 14, 17, 22)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 3}, 8: {8, 9}, 10: {10, 11}, 12: {12, 13}, 14: {14},
            15: {15, 16}, 17: {17}, 18: {18, 19}, 20: {20, 21}, 22: {22},
            23: {23, 24}
            }
        return slab, param, info

    def case_square_cm_1m1(self):
        """Return a square slab with cm[1 -1] plane group."""
        slab, param, info = self.case_square_cmm()
        tilt_c_axis(slab, (1, -1))
        info.symmetry.hermann = 'cm'
        info.symmetry.on_axes = ()
        info.symmetry.on_planes = (1, 15, 16, 22)
        info.symmetry.link_groups = {
            1: {1}, 2: {2, 9}, 3: {3, 8}, 10: {10, 13}, 11: {11, 12},
            14: {14, 17}, 15: {15}, 16: {16}, 18: {18, 21}, 19: {19, 20},
            22: {22}, 23: {23, 24}
            }
        return slab, param, info
