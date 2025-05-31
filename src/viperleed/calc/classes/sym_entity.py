"""Module sym_entity of viperleed.calc.classes.

Part of the functionality defined here used to be in the slab module.
Defines the SymPlane class, which represent a mirror/glide operations
with a 2D location.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-02-16'
__license__ = 'GPLv3+'

import numpy as np

from viperleed.calc.lib.base import dist_from_line
from viperleed.calc.lib.coordinates import add_edges_and_corners


class SymPlane:
    """Candidate plane for a symmetry operation. 'ty' pre-defines a type
    (mirror or glide), 'index2' allows the (1,2) and (2,1) directions if True,
    and collapse moves pos into the (0,0) unit cell if True."""

    def __init__(self, pos, dr, abt, ty="none", index2=False, collapse=True):
        if collapse:  # collapse to (0,0) cell
            self.pos = np.dot(abt.T, (np.dot(np.linalg.inv(abt.T), pos) % 1.0))
        else:
            self.pos = pos
        self.dir = dr/np.linalg.norm(dr)
        # normalized vector perpendicular to pos = in-plane
        self.type = ty
        self.par = []
        optionlist = [(1, 0), (0, 1), (1, 1), (1, -1)]
        if index2:
            optionlist.extend([(2, 1), (1, 2)])
        for (i, j) in optionlist:
            if abs((abs(np.dot(self.dir, (i*abt[0]+j*abt[1])))
                    / (np.linalg.norm(self.dir)
                    * np.linalg.norm(i*abt[0]+j*abt[1])))-1.0) < 0.001:
                self.par = np.array([i, j])

    def __str__(self):
        """Return a string representation of this SymPlane."""
        return f'SymPlane(pos={self.pos}, par={self.par})'

    @property
    def is_glide(self):
        """Return whether this is a mirror plane."""
        return self.type == 'glide'

    @property
    def is_mirror(self):
        """Return whether this is a mirror plane."""
        return self.type == 'mirror'

    @property
    def normal(self):
        """Return a unit vector normal to this plane (even without perp)."""
        return np.array((self.dir[1], -self.dir[0]))

    def distanceFromOrigin(self, abt):
        pointlist = [(0, 0), (1, 0), (0, 1), (1, 1)]
        return min([dist_from_line(self.pos, self.pos+self.dir,
                                   p[0]*abt[0]+p[1]*abt[1])
                    for p in pointlist])

    def isEquivalent(self, pl2, abt, eps=0.001):
        """Checks whether two symmetry planes have the same position and
        direction (including duplicates in next unit cell)"""
        if not np.array_equal(self.par, pl2.par):
            return False

        # If we're close to an edge or corner, also check translations
        fpos = np.dot(np.linalg.inv(np.transpose(abt)), self.pos) % 1.0
        releps = eps / np.linalg.norm(abt, axis=1)
        complist, _ = add_edges_and_corners([self.pos], (fpos,), releps, abt)

        for p in complist:
            if dist_from_line(pl2.pos, pl2.pos+pl2.dir, p) < eps:
                return True
        return False

    def point_operation(self, n_dim=2):
        """Return a matrix for mirroring across this plane."""
        # See https://en.wikipedia.org/wiki/Householder_transformation
        normal = self.normal
        mirror = np.identity(n_dim)
        mirror[:2, :2] -= 2*np.outer(normal, normal)
        return mirror
