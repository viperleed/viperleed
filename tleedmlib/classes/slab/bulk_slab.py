# -*- coding: utf-8 -*-
"""Module bulk_slab of viperleed.tleedmlib.classes.slab.

Created on 2023-02-21, originally Jun 13 2019

@author: Florian Kraushofer (@fkraushofer)
@author: Michele Riva (@michele-riva)

Defines the BulkSlab class, a BaseSlab describing a 3D-periodic crystal.
This module was created as part of the refactoring of slab.py.
"""

import copy
import itertools

import numpy as np
import scipy.spatial as sps

from viperleed.tleedmlib.base import angle, rotation_matrix
from viperleed.tleedmlib.base import rotation_matrix_order

from .base_slab import BaseSlab


class BulkSlab(BaseSlab):
    """A class representing an infinite solid, period in 3D.

    Contains unit cell, element information and atom coordinates.
    Also has a variety of convenience functions for manipulating
    and updating the atoms.

    Attributes
    ----------
    ucell : np.array
        The unit cell as vectors a, b, c (columns)
    poscar_scaling : float
        The original scaling factor from POSCAR
    elements : tuple of str
        Element labels as read from POSCAR
    chemelem : set of str
        Chemical elements in the slab, including from `ELEMENT_MIX`
    n_per_elem : dict {str: int}
        The number of atoms per POSCAR element.
    atlist : list of Atom
        List of all atoms in the slab.
    layers : list of Layer
        List of Layer objects, where each `layer` is a composite
        of sublayers, as in TensErLEED
    sublayers : list of Layer
        List of Layer objects, each containing atoms of equal
        element and Z coordinate
    sitelist : list of Sitetype
        List of distinct sites as Sitetype, storing information
        on vibration and concentration
    ucell_mod : list of tuples (str, numpy.ndarray)
        Stored modifications made to the unit cell; each is a tuple
        of (type, array), where type is 'lmul', 'rmul', or 'add'
    topat_ori_z : float
        Stores the original position of the topmost atom in Cartesian
        coordinates
    celltype : str                                                              # TODO: would be nicer with an Enum
        Unit-cell shape as string. Values: 'oblique', 'rhombic',
        'rectangular', 'square', 'hexagonal'
    planegroup : str
        Symmetry group of the slab. May be reduced by the user
        relative to `foundplanegroup`.
    foundplanegroup : str
        Highest symmetry found. Doesn't get modified when user
        reduces symmetry manually.
    orisymplane : SymPlane
        Only stored if the `planegroup` is ambiguous as to which unit
        vector the symmetry plane at the origin is parallel to
    linklists : list of list of Atom
        List of lists of atoms which are linked by a symmetry operation
    symbaseslab : Slab or None                                                  # TODO: do we need one??
        Slab with the smallest in-plane unit-cell area that shows
        the full symmetry of the slab.
    bulk_screws : list of int
        Integer list of rotation orders present in the bulk.
    bulk_glides : list of SymPlane
        List of glide-symmetry planes present in the bulk.
    """

    def __init__(self):
        """Initialize instance."""
        super().__init__()
        del self.bulkslab
        self.bulk_screws = []
        self.bulk_glides = []

    @property
    def is_bulk(self):
        """Return whether this is a bulk slab."""
        return True

    def doubleBulkSlab(self):
        """Returns a copy of the bulk slab which is doubled in thickness."""
        ts = copy.deepcopy(self)
        bulkc = np.transpose(np.copy(ts.ucell))[2]
        bulkc[2] *= -1
        ts.update_cartesian_from_fractional(update_origin=True)
        tmplist = ts.atlist[:]
        for at in tmplist:
            at.duplicate()
            at.cartpos = at.cartpos + bulkc
        ts.ucell[:, 2] *= 2
        ts.collapse_cartesian_coordinates(update_origin=True)
        return ts

    def getBulk3Dstr(self):
        """Returns a one-line string containing information about the bulk
        screw axes and glide planes. Only to be used for bulk slabs. Format of
        the string is 'r(2, 4), m([1,1], [ 1,-1])'. If neither screw axes nor
        glide planes exist, returns string 'None'."""
        b3ds = ""
        if self.bulk_screws:
            b3ds += "r({})".format(", ".join([str(v)
                                              for v in self.bulk_screws]))
        if self.bulk_glides:
            if b3ds:
                b3ds += ", "
            b3ds += "m({})".format(", ".join([np.array2string(gp.par,
                                                              separator=",")
                                              for gp in self.bulk_glides]))
        if not b3ds:
            return "None"
        return b3ds

    def getBulkRepeat(self, rp):
        """Based on a pre-existing definition of the bulk, tries to identify
        a repeat vector for which the bulk matches the slab above. Returns that
        vector in cartesian coordinates, or None if no match is found."""
        raise NotImplementedError

    def getCandidateLayerPeriod(self, eps=1e-4):
        """For a bulk slab, find at what offsets the sublayers repeat,
        checking only element and number of atoms. Returns a list of integer
        offsets between layer indices that are potentially equivalent."""
        if self.n_sublayers < 2:
            return([])
        cl = []     # candidate layers
        h = self.ucell[2, 2]  # cell height; periodicity cannot go beyond h/2
        l0 = self.sublayers[0]
        nl = self.n_sublayers
        l0el = l0.atlist[0].el
        l0n = len(l0.atlist)
        for i, lay in enumerate(self.sublayers[1:]):
            if abs(lay.cartbotz - l0.cartbotz) > h/2 + eps:
                break
            if lay.atlist[0].el == l0el and len(lay.atlist) == l0n:
                cl.append(i+1)
        if len(cl) == 0:
            return([])
        i = 0
        while i < len(cl):
            wrong = False
            zoff = self.sublayers[cl[i]].cartbotz - self.sublayers[0].cartbotz
            for j in range(1, int(np.ceil(self.n_sublayers/2))):
                if (self.sublayers[(j + cl[i]) % nl].atlist[0].el
                        != self.sublayers[j].atlist[0].el):
                    wrong = True
                    break
                if abs(zoff - ((self.sublayers[(j + cl[i]) % nl].cartbotz
                                - self.sublayers[j].cartbotz) % h)) > eps:
                    wrong = True
                    break
            if wrong:
                cl.pop(i)
            else:
                i += 1
        return(cl)

    def getMinC(self, rp, z_periodic=True):
        """Checks whether there is a vector c with a smaller length than
        the current one. If so, returns the minimized vector, else returns
        None."""
        eps = rp.SYMMETRY_EPS
        pcands = self.getCandidateLayerPeriod(eps)
        if len(pcands) == 0:
            return None
        ts = copy.deepcopy(self)
        ts.update_cartesian_from_fractional()
        ts.create_sublayers(eps)
        baseLayer = ts.sublayers[0]
        baseInd = ts.sublayers.index(baseLayer)
        nl = ts.n_sublayers
        ori = baseLayer.atlist[0].cartpos  # compare displacements from here
        repeatC = None
        for per in pcands:
            ind = (baseInd + per) % nl
            for at in ts.sublayers[ind].atlist:
                v = ori - at.cartpos
                if ts.isTranslationSymmetric(v, eps, z_periodic=z_periodic):
                    repeatC = at.cartpos - ori
                    break
            if repeatC is not None:
                break
        if repeatC is None:
            return None
        # optimize C vector to be close to Z, if possible
        cFracBase = np.dot(np.linalg.inv(ts.ab_cell), repeatC[:2]) % 1.0
        newC = np.append(np.dot(ts.ab_cell, cFracBase), -repeatC[2])
        for (i, j) in [(0, -1), (-1, 0), (-1, -1)]:
            v = np.dot(ts.ab_cell, cFracBase + np.array([i, j]))
            if (np.linalg.norm(np.append(v, -repeatC[2]))
                    < np.linalg.norm(newC)):
                newC[:2] = v
        return newC

    def isBulkGlideSymmetric(self, symplane, sldisp, eps):
        """Evaluates whether the bulk has a glide plane along a given
        direction, i.e. mirror at this direction, then some translation."""
        m = np.identity(3, dtype=float)
        ang = angle(symplane.dir, np.array([1, 0]))
        rotm = rotation_matrix(ang)
        m[:2, :2] = np.dot(np.linalg.inv(rotm),
                           np.dot(np.array([[1, 0], [0, -1]]), rotm))
        return self.isBulkTransformSymmetric(m, sldisp, eps)

    def isBulkScrewSymmetric(self, order, sldisp, eps):
        """Evaluates whether the slab has a screw axis of the given order when
        translated by the given number of sublayers."""
        m = rotation_matrix_order(order, dim=3)
        return self.isBulkTransformSymmetric(m, sldisp, eps)

    def isBulkTransformSymmetric(self, matrix, sldisp, eps):
        """Evalues whether the slab is self-equivalent under a given symmetry
        operation, and subsequent translation by a given number of sublayers"""
        self.check_a_b_in_plane()
        uc = self.ucell
        uct = np.transpose(uc)
        releps = [eps / np.linalg.norm(uct[j]) for j in range(0, 3)]
        # get translation vectors to check
        transVecs = []
        lowocclayer = self.fewest_atoms_sublayer
        baseInd = self.sublayers.index(lowocclayer)
        ori = lowocclayer.atlist[0].cartpos
        for at in self.sublayers[(baseInd + sldisp) % self.n_sublayers].atlist:
            transVecs.append((at.cartpos - np.dot(matrix, ori)).reshape(3, 1))
        for (i, sl) in enumerate(self.sublayers):
            coordlist = [at.cartpos for at in sl.atlist]
            oricm = np.array(coordlist)  # original cartesian coordinate matrix
            oripm = np.dot(np.linalg.inv(uc), oricm.transpose()) % 1.0
            # collapse (relative) coordinates to base unit cell
            oricm = np.dot(uc, oripm).transpose()
            # original cartesian coordinates collapsed to base unit cell
            transcoords = np.copy(oricm).transpose()
            transcoords = np.dot(matrix, transcoords)
            # now get coordinates of the sublayer to compare to
            sl2 = self.sublayers[(i + sldisp) % self.n_sublayers]
            oricm2 = np.array([at.cartpos for at in sl2.atlist])
            oripm2 = np.dot(np.linalg.inv(uc), oricm2.transpose()) % 1.0
            oricm2 = np.dot(uc, oripm2).transpose()
            # for every point in matrix, check whether is equal:
            for (i, p) in enumerate(oripm2.transpose()):
                # get extended comparison list for edges/corners:
                addlist = []
                for j in range(0, 3):
                    if abs(p[j]) < releps[j]:
                        addlist.append(oricm2[i]+uct[j])
                    if abs(p[j]-1) < releps[j]:
                        addlist.append(oricm2[i]-uct[j])
                if len(addlist) == 2:
                    # 2D coner - add the diagonally opposed point
                    addlist.append(addlist[0]+addlist[1]-oricm2[i])
                elif len(addlist) == 3:
                    # 3D corner - add all diagonally opposed points
                    addlist.extend([(p1 + p2 - oricm2[i]) for (p1, p2) in
                                    itertools.combinations(addlist, 2)])
                    addlist.append(addlist[0] + addlist[1] + addlist[2]
                                   - 2*oricm2[i])
                for v in addlist:
                    oricm2 = np.concatenate((oricm2, v.reshape(1, 3)))
            j = 0
            while j < len(transVecs):
                v = transVecs[j]
                shiftm = np.tile(v, len(coordlist))
                tmpcoords = transcoords + shiftm
                tmpcoords = np.dot(uc, (np.dot(np.linalg.inv(uc), tmpcoords)
                                        % 1.0))
                distances = sps.distance.cdist(tmpcoords.transpose(), oricm2,
                                               'euclidean')
                mismatch = False
                for sublist in distances:
                    if min(sublist) > eps:
                        mismatch = True
                        break
                if mismatch:
                    transVecs.pop(j)
                else:
                    j += 1
            if len(transVecs) == 0:
                return False
        return True