"""Module sym_entity of viperleed.tleedmlib.classes.

Created on 2023-02-16

@author: Michele Riva (@michele-riva)
@author: Florian Kraushofer (@fkraushofer)

Part of the functionality defined here used to be in the slab module.
Defines SymEntity and its subclasses SymPlane and RotationAxis, which
represent symmetry operations with a 2D location.
"""

import numpy as np

from viperleed.tleedmlib.base import dist_from_line


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

    def distanceFromOrigin(self, abt):
        pointlist = [(0, 0), (1, 0), (0, 1), (1, 1)]
        return min([dist_from_line(self.pos, self.pos+self.dir,
                                   p[0]*abt[0]+p[1]*abt[1])
                    for p in pointlist])

    def __str__(self):
        return ("SymPlane(pos = {}, par = {})".format(self.pos, self.par))

    def isEquivalent(self, pl2, abt, eps=0.001):
        """Checks whether two symmetry planes have the same position and
        direction (including duplicates in next unit cell)"""
        if not np.array_equal(self.par, pl2.par):
            return False
        complist = [self.pos]
        fpos = np.dot(np.linalg.inv(np.transpose(abt)), self.pos) % 1.0
        # if we're close to an edge or corner, also check translations
        for i in range(0, 2):
            releps = eps / np.linalg.norm(abt[i])
            if abs(fpos[i]) < releps:
                complist.append(self.pos+abt[i])
            if abs(fpos[i]-1) < releps:
                complist.append(self.pos-abt[i])
        if len(complist) == 3:  # coner - add the diagonally opposed one
            complist.append(complist[1]+complist[2]-complist[0])

        for p in complist:
            if dist_from_line(pl2.pos, pl2.pos+pl2.dir, p) < eps:
                return True
        return False