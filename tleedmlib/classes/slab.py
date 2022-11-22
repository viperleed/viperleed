# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Class accumulating atoms and layers, listing their elements and various other
properties. Includes functions for manipulation of those properties.
"""

import logging
import copy
import re
import itertools
from numbers import Real
from typing import Sequence

import numpy as np
import scipy.spatial as sps
from scipy.spatial import KDTree

try:
    import ase
    has_ase = True
except ImportError:
    has_ase = False

from viperleed.tleedmlib.base import (angle, rotation_matrix_order,
                                      rotation_matrix, dist_from_line,
                                      make_unique_list)
from viperleed.tleedmlib.classes.atom import Atom
import viperleed.tleedmlib as tl

logger = logging.getLogger("tleedm.slab")


class SymPlane:
    """Candidate plane for a symmetry operation."""

    def __init__(self, pos, dr, abt, ty="none", index2=True, collapse=True):
        """Initialize instance.
        
        Parameters
        ----------
        pos : np.ndarray
            Cartesian position of this symmetry plane. Can later be accessed
            as the .pos attribute. If collapse=True, self.pos contains the
            collapsed version of pos.
        dr : np.ndarray                                                                 # TODO: given the use in symmetry.py, this could become fractional
            In-plane vector identifying this plane, in cartesian coordinates.
            A unit vector parallel to dr can be accessed via .dir; a version
            in fractional coordinates is available in .par, instead.
        abt : np.ndarray
            2D unit cell of the slab for which this plane is a candidate
            symmetry plane. Unit vectors are rows, i.e., a, b = abt.                    # TODO: is this correct? Could change this name to "ucell" or similar
        ty : {'none', 'mirror', 'glide'}, optional
            Pre-defines a type for this symmetry plane. This property can
            later be accessed via the .type attribute. Default is 'none'.
        index2 : bool, optional
            Allows also (2,1) and (1,2) directions if True. When False,
            only (1,0), (0,1), (1,1) and (1,-1) directions are considered.
            Default is True.
        collapse : bool, optional
            Moves pos into the (0,0) unit cell if True. Default is True.

        Returns
        -------
        None.
        """
        self.dir = dr/np.linalg.norm(dr)
        self.type = ty
        self.pos = pos
        if collapse:  # collapse to (0,0) cell
            self.collapse(abt)

        # Identify fractional directions .par and .perp                                 # TODO: this could be done as (i, j) = np.dot(vec, inv(abt))
        self.par = []                                                                   # and would give general frac directions, not necessarily ints
        self.perp = []                                                                  # The current implementation can give errors should this identification
        optionlist = [(1, 0), (0, 1), (1, 1), (1, -1)]                                  # fail as then .par/.perp are '[]' (e.g., get_implied fails). Giving dr
        if index2:                                                                      # as fractional would solve at least the 'par' part.
            optionlist.extend([(2, 1), (1, 2)])

        perp_dir = np.dot(rotation_matrix(np.pi/2), self.dir)
        for (i, j) in optionlist:
            cart_dir = i * abt[0] + j * abt[1]
            cart_dir /= np.linalg.norm(cart_dir)
            if abs(abs(np.dot(self.dir, cart_dir)) - 1.0) < 0.001:
                self.par = np.array([i, j])
            if abs(abs(np.dot(perp_dir, cart_dir)) - 1.0) < 0.001:
                self.perp = np.array([i, j])

    def __str__(self):
        return (f"SymPlane(pos={self.pos}, par={self.par}, "
                f"type={self.type})")

    def collapse(self, abt):
        """Collapse self.pos to the (0, 0) unit cell."""
        # We need only collapse the component of pos perpendicular to the
        # plane, as the parallel one does not matter, and can be set to 0
        perp_dir = np.dot(rotation_matrix(np.pi/2), self.dir)
        self.pos = np.dot(self.pos, perp_dir) * perp_dir
        self.pos = np.dot(abt.T,
                          (np.dot(np.linalg.inv(abt.T), self.pos) % 1.0))

    def get_implied(self, abt, collapse=False):
        """Return the second symmetry plane that is implied by this one."""
        return SymPlane(pos=self.pos + np.dot(self.perp, abt)/2,
                        dr=self.dir, abt=abt, ty=self.type, index2=True,
                        collapse=collapse)

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
            if tl.base.dist_from_line(pl2.pos, pl2.pos+pl2.dir, p) < eps:
                return True
        return False


class Slab:
    """
    Contains unit cell, element information and atom coordinates. Also has a
    variety of convenience functions for manipulating and updating the atoms.
    Slabs can be created from an ase.Atoms object by passing an ase.Atoms object.

    Attributes
    ----------
    ucell : np.array
        The unit cell as vectors a, b, c (columns)
    poscar_scaling : float
        The original scaling factor from POSCAR
    elements : list of str
        Element labels as read from POSCAR
    chemelem : list of str
        Chemical elements in the slab, including from ELEMENT_MIX
    n_per_elem : dict {str: int}
        The number of atoms per element.
    last_element_mix : list of str
        Stores the last value of the ELEMENT_MIX parameter that was applied.
    atlist : list of Atom
        List of all atoms in the slab.
    layers : list of Layer
        List of Layer objects, where "layer" is a composite of sublayers, as
        in TensErLEED
    sublayers : list of Layer
        List of Layer objects, each containing atoms of equal element and Z
        coordinate
    sitelist : list of Sitetype
        List of distinct sites as Sitetype, storing information on vibration
        and concentration
    ucell_mod : list of tuples (str, np.array)
        Stored modifications made to the unit cell; each is a tuple of
        (type, matrix), where type is lmul, rmul, or add
    topat_ori_z : float
        Stores the original position of the topmost atom in cartesian
        coordinates
    celltype : str
        Unit cell type as string. Values: oblique, rhombic, rectangular,
        square, hexagonal
    planegroup : str
        Symmetry group of the slab, as string
    foundplanegroup : str
        Highest symmetry found, doesn't get modified when user reduces
        symmetry manually
    orisymplane : SymPlane
        Only stored if the planegroup is ambigious as to which unit vector the
        symmetry plane at the origin is parallel to
    linklists : list of list of Atom
        List of lists of atoms which are linked by a symmetry operation
    displists : list of list of Atom
        List of lists of atoms which are displaced together. This differs from
        'linklists' in that while linklists stores the linking based on the
        current symmetry, the 'displists' store actual displacement linking
        based on the symmetry set at the time the displacement is assigned.
    sites_initialized : bool
        Set by self.initSites
    layers_initialized : bool
        Set by self.createLayers
    deltas_initialized : bool
        Set by Rparams.generateSearchPars
    preprocessed : bool
        True if the POSCAR that his slab was read from had the
        'Plane group = XY' comment in the header, indicating that it was
        processed by tleedm before.
    symbaseslab : Slab
        Slab object collapsed to the base unit cell
    bulkslab : Slab
        Slab object containing only bulk layers
    bulk_screws : list of int
        Only assigned to the bulkslab object. Integer list of rotation orders
        present in the bulk.
    bulk_glides : list of SymPlane
        Only assigned to the bulkslab object. List of symmetry planes present
        in the bulk.
    """

    def __init__(self, ase_atoms=None):
        global has_ase

        self.ucell = np.array([])
        self.poscar_scaling = 1.
        self.elements = []
        self.chemelem = []
        self.n_per_elem = {}
        self.last_element_mix = None
        self.atlist = []
        self.layers = []
        self.sublayers = []
        self.sitelist = []
        self.ucell_mod = []
        self.ucell_ori = np.array([])
        self.topat_ori_z = None
        self.celltype = "unknown"
        self.planegroup = "unknown"
        self.foundplanegroup = "unknown"
        self.orisymplane = None
        self.linklists = []
        self.displists = []
        self.sites_initialized = False
        self.layers_initialized = False
        self.preprocessed = False
        self.deltas_initialized = False
        self.symbaseslab = None
        self.bulkslab = None
        self.bulk_screws = []
        self.bulk_glides = []

        if ase_atoms is None:
            return
        if not has_ase:
            logger.warning("Slab creation from ase.Atoms not supported: "
                           "Module ase not found.")
            return
        if type(ase_atoms) != ase.Atoms:
            raise TypeError("Slab object must be created empty or from "
                            "ase.Atoms, found unexpected type "
                            + str(type(ase_atoms)))
        # initialize from ase_atoms
        self.ucell = np.transpose(ase_atoms.cell[:])
        elems = [v.capitalize() for v in ase_atoms.get_chemical_symbols()]
        self.elements = make_unique_list(elems)
        self.n_per_elem = {k: elems.count(k) for k in self.elements}
        for i, (el, pos) in enumerate(zip(elems,
                                          ase_atoms.get_scaled_positions())):
            self.atlist.append(Atom(el, pos, i+1, self))
        self.getCartesianCoordinates()

    def resetSymmetry(self):
        """Sets all symmetry information back to default values."""
        self.ucell_mod = []
        self.ucell_ori = self.ucell
        self.celltype = "unknown"
        self.planegroup = "unknown"
        self.foundplanegroup = "unknown"
        self.orisymplane = None

    def resetAtomOriN(self):
        """Gets new 'original' numbers for atoms in the slab. If a bulkslab
        is defined, also updates the numbers there to keep the two consistent.
        """
        self.sortOriginal()
        self.sortByEl()
        bulkAtsRenumbered = []
        for (i, at) in enumerate(self.atlist):
            if self.bulkslab is not None:
                for bat in [a for a in self.bulkslab.atlist
                            if a.oriN == at.oriN
                            and a not in bulkAtsRenumbered]:
                    bat.oriN = i+1
                    bulkAtsRenumbered.append(bat)
            at.oriN = i+1

    def fullUpdate(self, rparams):
        """readPOSCAR initializes the slab with information from POSCAR;
        fullUpdate re-initializes the atom list, then uses the information
        from the parameters file to create layers, calculate cartesian
        coordinates (absolute and per layer), and to update elements and
        sites."""
        self.collapseFractionalCoordinates()
        self.getCartesianCoordinates()
        if not self.layers_initialized:
            self.createLayers(rparams)
        self.updateElements(rparams)
        if not self.sites_initialized:
            self.initSites(rparams)
        if rparams.fileLoaded["VIBROCC"]:
            for at in self.atlist:
                at.initDisp()

    def getCartesianCoordinates(self, updateOrigin=False):
        """Assigns absolute cartesian coordinates to all atoms, with x,y using
        the unit cell (top plane), while z = 0 for the topmost atom and
        positive going down through the slab. If updateOrigin is set True, the
        cartesian origin relative to the fractional origin will be updated,
        otherwise it is static."""
        al = self.atlist[:]     # temporary copy
        al.sort(key=lambda atom: atom.pos[2])
        topat = al[-1]
        topcart = np.dot(self.ucell, topat.pos)
        if updateOrigin or self.topat_ori_z is None:
            self.topat_ori_z = topcart[2]
        for atom in al:
            atom.cartpos = np.dot(self.ucell, atom.pos)
            atom.cartpos[2] = self.topat_ori_z - atom.cartpos[2]

    def getFractionalCoordinates(self):
        """Calculates fractional coordinates for all atoms from their
        cartesian coordinates, using the slab unit cell."""
        uci = np.linalg.inv(self.ucell)
        for at in self.atlist:
            tp = np.copy(at.cartpos)
            tp[2] = self.topat_ori_z-tp[2]
            at.pos = np.dot(uci, tp)

    def collapseFractionalCoordinates(self):
        """Finds atoms outside the parallelogram spanned by the unit vectors
        a and b and moves them inside."""
        for at in self.atlist:
            at.pos = at.pos % 1.0

    def collapseCartesianCoordinates(self, updateOrigin=False):
        """Finds atoms outside the parallelogram spanned by the unit vectors
        a and b and moves them inside. If keepOriZ is True, the old value of
        the top atom position will be preserved.."""
        self.getFractionalCoordinates()
        self.collapseFractionalCoordinates()
        self.getCartesianCoordinates(updateOrigin=updateOrigin)

    def updateLayerCoordinates(self):
        for layer in self.layers:
            layer.getLayerPos()

    def createLayers(self, rparams, bulk_cuts=[]):
        """Creates a list of Layer objects based on the N_BULK_LAYERS and
        LAYER_CUTS parameters in rparams. If layers were already defined,
        overwrite. The bulk_cuts kwarg allows specifically inserting
        automatically detected bulk layer cuts. Returns the cuts as a sorted
        list of floats."""
        # first interpret LAYER_CUTS parameter - can be a list of strings
        ct = []
        rgx = re.compile(r'\s*(dz|dc)\s*\(\s*(?P<cutoff>[0-9.]+)\s*\)')
        al = self.atlist[:]
        al.sort(key=lambda atom: atom.pos[2])
        for (i, s) in enumerate(rparams.LAYER_CUTS):
            if type(s) == float:
                ct.append(s)
                continue
            s = s.lower()
            if "dz" in s or "dc" in s:
                m = rgx.match(s)
                if not m:
                    logger.warning("Error parsing part of LAYER_CUTS: " + s)
                    continue
                cutoff = float(m.group('cutoff'))
                lowbound = 0.
                if bulk_cuts:
                    lowbound = max(bulk_cuts)
                highbound = 1.
                val = None
                if (i > 1) and (rparams.LAYER_CUTS[i-1] in ["<", ">"]):
                    try:
                        val = float(rparams.LAYER_CUTS[i-2])
                    except ValueError:
                        logger.warning("LAYER_CUTS: Error parsing left-hand "
                                       "boundary for " + s)
                    if val is not None:
                        if rparams.LAYER_CUTS[i-1] == "<":
                            lowbound = val
                        else:
                            highbound = val
                if i < len(rparams.LAYER_CUTS) - 2 and (rparams.LAYER_CUTS[i+1]
                                                        in ["<", ">"]):
                    try:
                        val = float(rparams.LAYER_CUTS[i+2])
                    except ValueError:
                        logger.warning("LAYER_CUTS: Error parsing right-hand "
                                       "boundary for " + s)
                    if val is not None:
                        if rparams.LAYER_CUTS[i+1] == ">":
                            lowbound = val
                        else:
                            highbound = val
                if 'dc' in s:
                    cutoff *= (self.ucell[2, 2]
                               / np.linalg.norm(self.ucell[:, 2]))
                for i in range(1, len(al)):
                    if ((abs(al[i].cartpos[2]-al[i-1].cartpos[2]) > cutoff)
                            and al[i].pos[2] > lowbound
                            and al[i].pos[2] < highbound
                            and al[i-1].pos[2] > lowbound
                            and al[i-1].pos[2] < highbound):
                        ct.append(abs((al[i].pos[2]+al[i-1].pos[2])/2))
            elif s not in ["<", ">"]:
                try:
                    ct.append(float(s))
                except ValueError:
                    logger.warning("LAYER_CUTS: Could not parse value: " + s)
                    continue
        if bulk_cuts:
            ct = [v for v in ct if v > max(bulk_cuts) + 1e-6] + bulk_cuts
        ct.sort()
        self.layers = []
        tmplist = self.atlist[:]
        self.sortByZ()
        laynum = 0
        b = True if rparams.N_BULK_LAYERS > 0 else False
        newlayer = tl.Layer(self, 0, b)
        self.layers.append(newlayer)
        for atom in self.atlist:
            # only check for new layer if we're not in the top layer already
            if laynum < len(ct):
                if atom.pos[2] > ct[laynum]:
                    # if atom is higher than the next cutoff, make a new layer
                    laynum += 1
                    b = True if rparams.N_BULK_LAYERS > laynum else False
                    newlayer = tl.Layer(self, laynum, b)
                    self.layers.append(newlayer)
                    check = True    # check for empty layer
                    while check:
                        if laynum >= len(ct):
                            check = False
                        elif atom.pos[2] <= ct[laynum]:
                            check = False
                        else:
                            laynum += 1
                            b = (True if rparams.N_BULK_LAYERS > laynum
                                 else False)
                            newlayer = tl.Layer(self, laynum, b)
                            self.layers.append(newlayer)
            atom.layer = newlayer
            newlayer.atlist.append(atom)
        dl = []
        for layer in self.layers:
            if not layer.atlist:
                logger.warning('A layer containing no atoms was found. Layer '
                               'will be deleted. Check LAYER_CUTS parameter.')
                rparams.setHaltingLevel(2)
                dl.append(layer)
        for layer in dl:
            if layer.isBulk:
                self.layers[layer.num+1].isBulk = True
            self.layers.remove(layer)
            del layer
        self.layers.reverse()
        for i, layer in enumerate(self.layers):
            layer.getLayerPos()
            layer.num = i
        self.atlist = tmplist
        self.layers_initialized = True
        return ct

    def createSublayers(self, eps=0.001):
        """Sorts the atoms in the slab into sublayers, sorted by element and Z
        coordinate."""
        self.sortByZ()
        subl = []  # will be a list of sublayers, using the Layer class
        for el in self.elements:
            sublists = [[a for a in self.atlist if a.el == el]]
            # first, split at points where two atoms are more than eps apart
            i = 0
            while i < len(sublists):
                brk = False
                if len(sublists[i]) > 1:
                    tmplist = sublists[i][:]
                    for j in range(1, len(tmplist)):
                        if (abs(tmplist[j].cartpos[2]
                                - tmplist[j-1].cartpos[2]) > eps):
                            sublists.append(tmplist[:j])
                            sublists.append(tmplist[j:])
                            sublists.pop(i)
                            brk = True
                            break
                    if not brk:
                        i += 1
                else:
                    i += 1
            # now, go through again and split sublayers at greatest interlayer
            #   distance, if they are too thick overall
            i = 0
            while i < len(sublists):
                brk = False
                if len(sublists[i]) > 1:
                    if abs(sublists[i][0].cartpos[2]
                           - sublists[i][-1].cartpos[2]) > eps:
                        maxdist = abs(sublists[i][1].cartpos[2]
                                      - sublists[i][0].cartpos[2])
                        maxdistindex = 1
                        for j in range(2, len(sublists[i])):
                            d = abs(sublists[i][j].cartpos[2]
                                    - sublists[i][j-1].cartpos[2])
                            if d > maxdist:
                                maxdist = d
                                maxdistindex = j
                        sublists.append(sublists[i][:maxdistindex])
                        sublists.append(sublists[i][maxdistindex:])
                        sublists.pop(i)
                        brk = True
                else:
                    i += 1
                if not brk:
                    i += 1
            # now, create sublayers based on sublists:
            for ls in sublists:
                newsl = tl.Layer(self, 0, sublayer=True)
                subl.append(newsl)
                newsl.atlist = ls
                newsl.cartbotz = ls[0].cartpos[2]
        self.sublayers = []
        subl.sort(key=lambda sl: -sl.cartbotz)
        while subl:
            acc = [subl.pop()]  # accumulate sublayers with same z
            while subl:
                if abs(subl[-1].cartbotz - acc[0].cartbotz) < eps:
                    acc.append(subl.pop())
                else:
                    break
            acc.sort(key=lambda sl: sl.atlist[0].el)  # sort by element
            self.sublayers.extend(acc)
        for (i, sl) in enumerate(self.sublayers):
            sl.num = i

    def getMinLayerSpacing(self):
        """Returns the minimum distance (cartesian) between two layers in the
        slab. Returns zero if there is only one layer, or none are defined."""
        if len(self.layers) < 2:
            return 0
        self.getCartesianCoordinates()
        return min([(self.layers[i].cartori[2] - self.layers[i-1].cartbotz)
                    for i in range(1, len(self.layers))])

    def updateElements(self, rp):
        """Updates elements based on the ELEMENT_MIX parameter, and warns in
        case of a naming conflict."""
        if self.last_element_mix == rp.ELEMENT_MIX:
            return     # don't update if up to date
        # update nelem
        c = 0
        oldels = self.elements[:]
        for i, pel in enumerate(oldels):
            if pel not in rp.ELEMENT_MIX:
                c += 1
            else:
                c += len(rp.ELEMENT_MIX[pel])
                # check for overlapping names:
                for el in rp.ELEMENT_MIX[pel]:
                    if el in oldels:
                        logger.warning(
                            "Element name "+el+" given in ELEMENT_MIX is also "
                            "an element name in POSCAR. It is recommended you "
                            "rename the element in the POSCAR file.")
        self.chemelem = []
        for el in self.elements:
            if el in rp.ELEMENT_MIX:
                self.chemelem.extend([e.capitalize()
                                      for e in rp.ELEMENT_MIX[el]
                                      if not e.capitalize() in self.chemelem])
            else:
                self.chemelem.append(el.capitalize())
        self.last_element_mix = rp.ELEMENT_MIX

    def updateElementCount(self):
        """Updates the number of atoms per element."""
        self.n_per_elem = {}
        del_elems = []
        for el in self.elements:
            n = len([at for at in self.atlist if at.el == el])
            if n > 0:
                self.n_per_elem[el] = n
            else:
                del_elems.append(el)
        for el in del_elems:
            self.elements.remove(el)

    def updateAtomNumbers(self):
        """Updates atom oriN - should not happen normally, but necessary if
        atoms get deleted."""
        for (i, at) in enumerate(self.atlist):
            at.oriN = i+1

    def initSites(self, rp):
        """Goes through the atom list and supplies them with appropriate
        SiteType objects, based on the SITE_DEF parameters from the supplied
        Rparams."""
        atlist = self.atlist[:]     # copy to not have any permanent changes
        atlist.sort(key=lambda atom: atom.oriN)
        sl = []
        for el in rp.SITE_DEF:
            for sitename in rp.SITE_DEF[el]:
                newsite = tl.Sitetype(el, sitename)
                sl.append(newsite)
                for i in rp.SITE_DEF[el][sitename]:
                    try:
                        if atlist[i-1].el != el:
                            logger.warning(
                                'SITE_DEF tries to assign atom number '
                                + str(i) + ' as ' + el + ', but POSCAR has it '
                                'as '+atlist[i-1].el+'. Atom will be skipped '
                                'and left as default site type!')
                            rp.setHaltingLevel(1)
                        else:
                            atlist[i-1].site = newsite
                    except IndexError:
                        logger.error('SITE_DEF: atom number out of bounds.')
                        raise
        for el in self.elements:
            newsite = tl.Sitetype(el, 'def')
            found = False
            for at in atlist:
                if at.el == el and at.site is None:
                    at.site = newsite
                    found = True
            if found:
                sl.append(newsite)
        for site in [s for s in sl if s.el in rp.ELEMENT_MIX]:
            site.mixedEls = rp.ELEMENT_MIX[site.el][:]
        self.sitelist = sl
        self.sites_initialized = True

    def sortByZ(self, botToTop=False):
        """Sorts atlist by z coordinate"""
        self.atlist.sort(key=lambda atom: atom.pos[2])
        if botToTop:
            self.atlist.reverse()

    def sortByEl(self):
        """Sorts atlist by elements, preserving the element order from the
        original POSCAR"""
        # unfortunately, simply calling the sort function by element does not
        #    preserve the element order from the POSCAR
        esortlist = sorted(self.atlist, key=lambda atom: atom.el)
        lastel = ''
        tmpElList = []
        isoLists = []
        # generate sub-lists isolated by elements
        for at in esortlist:
            if at.el != lastel:
                tmpElList.append(at.el)
                isoLists.append([])
                lastel = at.el
            isoLists[-1].append(at)
        sortedlist = []
        # going through the elements in the order they appear in POSCAR, find
        #   the corresponding index in tmpElList and append the atoms of that
        #   type to sorted list
        for el in self.elements:
            try:
                i = tmpElList.index(el)
            except ValueError:
                logger.error("Unexpected point encountered in Slab.sortByEl: "
                             "Could not find element in element list")
            else:
                sortedlist.extend(isoLists[i])
        self.atlist = sortedlist

    def sortOriginal(self):
        """Sorts atlist by original atom order from POSCAR"""
        self.atlist.sort(key=lambda atom: atom.oriN)

    def projectCToZ(self):
        """makes the c vector of the unit cell perpendicular to the surface,
        changing all atom coordinates to fit the new base"""
        if self.ucell[0, 2] != 0.0 or self.ucell[1, 2] != 0.0:
            self.getCartesianCoordinates()
            self.ucell[:, 2] = np.array([0, 0, self.ucell[2, 2]])
            self.collapseCartesianCoordinates()
            # implicitly also gets new fractional coordinates

    def restoreOriState(self, keepDisp=False):
        """Resets the atom positions and site vibrational amplitudes to the
        original state, and stores the deviations as offset instead."""
        for site in self.sitelist:
            if site.oriState is None:
                continue
            siteats = [at for at in self.atlist if at.site == site]
            for el in site.occ:
                for at in siteats:
                    o = site.vibamp[el] - site.oriState.vibamp[el]
                    if el in at.offset_vib:
                        o += at.offset_vib[el]
                    at.offset_vib[el] = o
                    if abs(at.offset_vib[el]) < 1e-6:
                        del at.offset_vib[el]
                site.vibamp[el] = site.oriState.vibamp[el]
                site.occ[el] = site.oriState.occ[el]
        uci = np.linalg.inv(self.ucell)
        self.getCartesianCoordinates()
        for at in self.atlist:
            if at.oriState is None:
                continue
            for el in at.disp_occ.keys():
                o = np.array([0., 0., 0.])
                if el in at.offset_geo:
                    o += at.offset_geo[el]
                if el in at.disp_geo_offset:
                    o += at.disp_geo_offset[el][0]
                else:
                    o += at.disp_geo_offset["all"][0]
                o[2] *= -1
                ro = np.dot(uci, o)
                nro = at.pos - at.oriState.pos + ro
                for i in range(0, len(nro)):  # minimize abs(nro)
                    while nro[i] > 0.5:
                        nro[i] -= 1
                    while nro[i] < -0.5:
                        nro[i] += 1
                no = np.dot(self.ucell, nro)
                no[2] *= -1
                at.offset_geo[el] = no
            at.pos = at.oriState.pos
            at.disp_geo_offset = {"all": [np.zeros(3)]}
        self.collapseFractionalCoordinates()
        self.getCartesianCoordinates()
        self.updateLayerCoordinates()
        if keepDisp:
            return
        for at in self.atlist:
            at.deltasGenerated = []
            at.initDisp(force=True)
            at.constraints = {1: {}, 2: {}, 3: {}}
        return

    def rotateAtoms(self, order, axis=np.array([0., 0.])):
        """Translates the atoms in the slab to have the axis in the origin,
        applies an order-fold rotation matrix to the atom positions, then
        translates back"""
        self.getCartesianCoordinates()
        m = rotation_matrix_order(order)
        for at in self.atlist:
            # translate origin to candidate point, rotate, translate back
            at.cartpos[0:2] = np.dot(m, at.cartpos[0:2] - axis) + axis
        self.getFractionalCoordinates()

    def rotateUnitCell(self, order, append_ucell_mod=True):
        """Rotates the unit cell (around the origin), leaving atom positions
        the same. Note that this rotates in the opposite direction as
        rotateAtoms."""
        self.getCartesianCoordinates()
        m = rotation_matrix_order(order)
        m3 = np.identity(3, dtype=float)
        m3[:2, :2] = m
        self.ucell = np.dot(m3, self.ucell)
        if append_ucell_mod:
            self.ucell_mod.append(('lmul', m3))
        self.getFractionalCoordinates()

    def mirror(self, symplane, glide=False):
        """Translates the atoms in the slab to have the symplane in the
        origin, applies a mirror or glide matrix, then translates back.
        Very inefficient implementation!"""
        ang = angle(symplane.dir, np.array([1, 0]))
        rotm = rotation_matrix(ang)
        rotmirm = np.dot(np.linalg.inv(rotm),
                         np.dot(np.array([[1, 0], [0, -1]]), rotm))
        # rotates to have plane in x direction, mirrors on x
        if glide:
            abt = self.ucell[:2, :2].T
            glidevec = (symplane.par[0]*abt[0]+symplane.par[1]*abt[1])/2
        else:
            glidevec = np.zeros(2)
        for at in self.atlist:
            at.cartpos[:2] -= symplane.pos     # translate to plane
            at.cartpos[:2] = np.dot(rotmirm, at.cartpos[:2])    # apply mirror
            at.cartpos[:2] += symplane.pos     # translate back
            at.cartpos[:2] += glidevec   # 0 if not glides

    def getLowOccLayer(self):
        """Finds and returns the lowest occupancy sublayer"""
        minlen = len(self.sublayers[0].atlist)
        lowocclayer = self.sublayers[0]
        for lay in self.sublayers:
            if len(lay.atlist) < minlen:
                lowocclayer = lay
                minlen = len(lay.atlist)
        return lowocclayer

    def isRotationSymmetric(self, axis, order, eps):
        """Evaluates whether the slab is equivalent to itself when rotated
        around the axis with the given rotational order"""
        m = rotation_matrix_order(order)
        ab = self.ucell[:2, :2]
        abt = ab.T
        releps = [eps / np.linalg.norm(abt[j]) for j in range(0, 2)]
        shiftv = axis.reshape(2, 1)
        for sl in self.sublayers:
            coordlist = [at.cartpos[0:2] for at in sl.atlist]
            shiftm = np.tile(shiftv, len(coordlist))
            # matrix to shift all coordinates by axis
            oricm = np.array(coordlist)  # original cartesian coordinate matrix
            oripm = np.dot(np.linalg.inv(ab), oricm.transpose()) % 1.0
            # collapse (relative) coordinates to base unit cell
            oricm = np.dot(ab, oripm).transpose()
            # original cartesian coordinates collapsed to base unit cell
            tmpcoords = np.copy(oricm).transpose()
            # copy of coordinate matrix to be rotated
            tmpcoords -= shiftm
            tmpcoords = np.dot(m, tmpcoords)
            tmpcoords += shiftm
            tmpcoords = np.dot(ab,
                               (np.dot(np.linalg.inv(ab), tmpcoords) % 1.0))
            # collapse coordinates to base unit cell
            # for every point in matrix, check whether is equal:
            for (i, p) in enumerate(oripm.transpose()):
                # get extended comparison list for edges/corners:
                addlist = []
                for j in range(0, 2):
                    if abs(p[j]) < releps[j]:
                        addlist.append(oricm[i]+abt[j])
                    if abs(p[j]-1) < releps[j]:
                        addlist.append(oricm[i]-abt[j])
                if len(addlist) == 2:
                    # coner - add the diagonally opposed one
                    addlist.append(addlist[0]+addlist[1]-oricm[i])
                for v in addlist:
                    oricm = np.concatenate((oricm, v.reshape(1, 2)))
            distances = sps.distance.cdist(tmpcoords.transpose(), oricm,
                                           'euclidean')
            for sublist in distances:
                if min(sublist) > eps:
                    return False
        return True

    def getCompareCoords(self, eps):
        """
        Generates a numpy array containing all atom coordinates, plus
        equivalent positions for atoms at the edge of the unit cell. Also
        returns a second numpy array containing the sublayer per atom, for
        filtering the coordinate array.

        Parameters
        ----------
        eps : float
            Error tolerance for positions (cartesian)

        Returns
        -------
        compare_coords : numpy.array
            2D array containing all atom coordiantes in order
        compare_sublayers : numpy.array
            1D array containing the sublayer index for each atom in the same
            order as compare_coords. For sublayer-wise comparison.

        """
        uc = self.ucell
        releps = tuple(eps / np.linalg.norm(uc.T[j, :2]) for j in range(0, 2))
        frac_coords = np.dot(np.linalg.inv(uc),
                             np.array([at.cartpos
                                       for at in self.atlist]).T) % 1.0
        compare_coords = np.dot(uc, frac_coords).T  # already collapsed
        compare_sublayers = np.array(
            [next(i for i, sl in enumerate(self.sublayers) if at in sl.atlist)
             for at in self.atlist])
        for (i, p) in enumerate(frac_coords.T):
            # get extended comparison list for edges/corners:
            addlist = []
            for j in range(0, 2):
                if abs(p[j]) < releps[j]:
                    addlist.append(compare_coords[i]+uc.T[j])
                if abs(p[j]-1) < releps[j]:
                    addlist.append(compare_coords[i]-uc.T[j])
            if len(addlist) == 2:
                # corner - add the diagonally opposed one
                addlist.append(addlist[0]+addlist[1]-compare_coords[i])
            for v in addlist:
                compare_coords = np.concatenate((compare_coords,
                                                 v.reshape(1, 3)))
                compare_sublayers = np.append(compare_sublayers,
                                              compare_sublayers[i])
        return compare_coords, compare_sublayers

    def isTranslationSymmetric_2D(self, tv, eps, compare_to=None):
        """
        Evaluates whether the slab is equivalent to itself under a given 2D
        translation. Slightly faster than isTranslationSymmetric.

        Parameters
        ----------
        tv : numpy.array
            2D or 3D translation vector; third component should be (near) zero.
        eps : float
            Error tolerance for positions (cartesian)
        compare_to : tuple, or None, optional
            when a tuple it should contain two elements:
                compare_coords : numpy.array
                    2D array containing cartesian coordiantes of all atoms,
                    including equivalent positions for atoms at the edge of a unit
                    cell
                compare_sublayers : numpy.array
                    1D array containing the sublayer index for each atom in the
                    same order as compare_coords. For sublayer-wise comparison.
            If None or not given, compare_to is self.getCompareCoords(eps).
            Default is None.

        Returns
        -------
        bool
            True if translation symmetric, else False.
        """
        if len(tv) not in (2, 3):
            raise ValueError("isTranslationSymmetric_2D: not a 2D or 3D "
                             "vector")
        uc = self.ucell
        if len(tv) == 2:
            tv = np.append(tv, 0.)
        # create 2D array of cart. coordinates, collapsed to base unit cell
        cart_coords = np.dot(
            uc, np.dot(np.linalg.inv(uc), np.array([at.cartpos for at
                                                    in self.atlist]).T)
            % 1.0).T
        shifted_coords = cart_coords.T + np.tile(tv.reshape(3, 1),
                                                 len(self.atlist))
        shifted_coords = np.dot(
            uc, (np.dot(np.linalg.inv(uc), shifted_coords) % 1.0))
        if compare_to is None:
            compare_coords, compare_sublayers = self.getCompareCoords(eps)
        else:
            compare_coords, compare_sublayers = compare_to
        sublayer_per_atom = np.array(
            [next(i for i, sl in enumerate(self.sublayers) if at in sl.atlist)
             for at in self.atlist])
        for i in range(len(self.sublayers)):  # sublayer-wise comparison
            distances = sps.distance.cdist(
                shifted_coords.T[sublayer_per_atom == i],
                compare_coords[compare_sublayers == i],
                'euclidean')
            if any(min(sublist) > eps for sublist in distances):
                return False
        return True

    def isTranslationSymmetric(self, tv, eps, z_periodic=True, z_range=None):
        """
        Evaluates whether the slab is equivalent to itself when translated
        along the given cartesian translation vector tv.

        Parameters
        ----------
        tv : numpy array
            2- or 3-dimensional translation vectors are accepted.
        eps : float
            Error tolerance for positions (cartesian)
        z_periodic : bool, optional
            True for checking periodicity of a bulk slab, in which the c vector
            is a true unit cell vector. False otherwise.
        z_range : tuple of floats, optional
            Limit check to only atoms within a given range of cartesian
            coordinates. The default is None.

        Returns
        -------
        bool
            True if translation symmetric, else False.
        """
        if len(tv) == 2:  # two-dimensional displacement. append zero for z
            tv = np.append(tv, 0.)
        uc = np.copy(self.ucell)
        uc[:, 2] *= -1   # mirror c vector down
        uct = uc.T
        releps = [eps / np.linalg.norm(uct[j]) for j in range(0, 3)]
        shiftv = tv.reshape(3, 1)
        # unlike in-plane operations, this one cannot be done sublayer-internal
        coordlist = [at.cartpos for at in self.atlist]
        shiftm = np.tile(shiftv, len(coordlist))
        oricm = np.array(coordlist)  # original cartesian coordinate matrix
        oricm[:, 2] *= -1
        shiftm[2] *= -1
        # collapse (relative) coordinates to base unit cell
        oripm = np.dot(np.linalg.inv(uc), oricm.T) % 1.0
        # original cartesian coordinates collapsed to base unit cell
        oricm = np.dot(uc, oripm).T
        # copy of coordinate matrix to be manipulated
        tmpcoords = np.copy(oricm).T
        # determine which z to check
        if z_range is None:
            min_z = np.min(tmpcoords[2]) - eps
            max_z = np.max(tmpcoords[2]) + eps
        else:
            z_range = tuple(-v for v in z_range)
            min_z, max_z = min(z_range) - eps, max(z_range) + eps
        tmpcoords += shiftm
        if not z_periodic:
            # discard atoms that moved out of range in z
            tmpcoords = tmpcoords[:, tmpcoords[2] >= min_z]
            tmpcoords = tmpcoords[:, tmpcoords[2] <= max_z]
        # collapse coordinates to base unit cell
        tmpcoords = np.dot(uc, (np.dot(np.linalg.inv(uc), tmpcoords) % 1.0))
        # for every point in matrix, check whether is equal:
        # use list here, convert to np.array later
        element_per_atom = [at.el for at in self.atlist]
        element_per_atom_extended = element_per_atom.copy()
        for (i, p) in enumerate(oripm.T):
            # get extended comparison list for edges/corners:
            addlist = []
            for j in range(0, 3):
                if abs(p[j]) < releps[j]:
                    addlist.append(oricm[i]+uct[j])
                if abs(p[j]-1) < releps[j]:
                    addlist.append(oricm[i]-uct[j])
            if len(addlist) == 2:
                # 2D coner - add the diagonally opposed point
                addlist.append(addlist[0]+addlist[1]-oricm[i])
            elif len(addlist) == 3:
                # 3D corner - add all diagonally opposed points
                addlist.extend([(p1 + p2 - oricm[i]) for (p1, p2) in
                                itertools.combinations(addlist, 2)])
                addlist.append(addlist[0] + addlist[1] + addlist[2]
                               - 2*oricm[i])
            for v in addlist:
                oricm = np.concatenate((oricm, v.reshape(1, 3)))
                element_per_atom_extended.append(element_per_atom[i])
        # need to convert to array here for use as mask
        element_per_atom = np.array(element_per_atom)
        element_per_atom_extended = np.array(element_per_atom_extended)
        for el in self.elements:  # element-wise comparison
            distances = sps.distance.cdist(
                tmpcoords.T[element_per_atom == el],
                oricm[element_per_atom_extended == el],
                'euclidean')
            if any(min(sublist) > eps for sublist in distances):
                return False
        return True

    def isBulkTransformSymmetric(self, matrix, sldisp, eps):
        """Evalues whether the slab is self-equivalent under a given symmetry
        operation, and subsequent translation by a given number of sublayers"""
        uc = self.ucell
        uct = np.transpose(uc)
        releps = [eps / np.linalg.norm(uct[j]) for j in range(0, 3)]
        # get translation vectors to check
        transVecs = []
        lowocclayer = self.getLowOccLayer()
        baseInd = self.sublayers.index(lowocclayer)
        ori = lowocclayer.atlist[0].cartpos
        for at in self.sublayers[(baseInd + sldisp)
                                 % len(self.sublayers)].atlist:
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
            sl2 = self.sublayers[(i + sldisp) % len(self.sublayers)]
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

    def isBulkScrewSymmetric(self, order, sldisp, eps):
        """Evaluates whether the slab has a screw axis of the given order when
        translated by the given number of sublayers."""
        m = rotation_matrix_order(order, dim=3)
        return self.isBulkTransformSymmetric(m, sldisp, eps)

    def isBulkGlideSymmetric(self, symplane, sldisp, eps):
        """Evaluates whether the bulk has a glide plane along a given
        direction, i.e. mirror at this direction, then some translation."""
        m = np.identity(3, dtype=float)
        ang = angle(symplane.dir, np.array([1, 0]))
        rotm = rotation_matrix(ang)
        m[:2, :2] = np.dot(np.linalg.inv(rotm),
                           np.dot(np.array([[1, 0], [0, -1]]), rotm))
        return self.isBulkTransformSymmetric(m, sldisp, eps)

    def isMirrorSymmetric(self, symplane, eps, glide=False):
        """Evaluates whether the slab is equivalent to itself when applying a
        mirror or glide operation at a given plane"""
        ang = angle(symplane.dir, np.array([1, 0]))
        rotm = rotation_matrix(ang)
        rotmirm = np.dot(np.linalg.inv(rotm),
                         np.dot(np.array([[1, 0], [0, -1]]), rotm))
        # rotates to have plane in x direction, mirrors on x
        ab = self.ucell[:2, :2]
        abt = ab.T
        releps = [eps / np.linalg.norm(abt[j]) for j in range(0, 2)]
        shiftv = symplane.pos.reshape(2, 1)
        if glide:
            glidev = ((symplane.par[0]*abt[0]+symplane.par[1]*abt[1])
                      / 2).reshape(2, 1)
        for sl in self.sublayers:
            coordlist = [at.cartpos[:2] for at in sl.atlist]
            shiftm = np.tile(shiftv, len(coordlist))  # shift all coordinates
            if glide:
                glidem = np.tile(glidev, len(coordlist))
            oricm = np.array(coordlist)  # original cartesian coordinate matrix
            oripm = np.dot(np.linalg.inv(ab), oricm.transpose()) % 1.0
            # collapse (relative) coordinates to base unit cell
            oricm = np.dot(ab, oripm).transpose()
            # original cartesian coordinates collapsed to base unit cell
            tmpcoords = np.copy(oricm).transpose()
            # copy of coordinate matrix to be rotated
            tmpcoords -= shiftm
            tmpcoords = np.dot(rotmirm, tmpcoords)
            tmpcoords += shiftm
            if glide:
                tmpcoords += glidem
            tmpcoords = np.dot(ab, (np.dot(np.linalg.inv(ab),
                                           tmpcoords) % 1.0))
            # collapse coordinates to base unit cell
            # for every point in matrix, check whether is equal:
            for (i, p) in enumerate(oripm.transpose()):
                # get extended comparison list for edges/corners:
                addlist = []
                for j in range(0, 2):
                    if abs(p[j]) < releps[j]:
                        addlist.append(oricm[i]+abt[j])
                    if abs(p[j]-1) < releps[j]:
                        addlist.append(oricm[i]-abt[j])
                if len(addlist) == 2:
                    # coner - add the diagonally opposed one
                    addlist.append(addlist[0]+addlist[1]-oricm[i])
                for v in addlist:
                    oricm = np.concatenate((oricm, v.reshape(1, 2)))
            distances = sps.distance.cdist(tmpcoords.T, oricm,
                                           'euclidean')
            for sublist in distances:
                if min(sublist) > eps:
                    return False
        return True

    def isEquivalent(self, slab, eps=0.001):
        """Compares the slab to another slab, returns True if all atom cartpos
        match (with at least one other atom, if there are duplicates), False
        if not. Both slabs are copied and collapsed to the (0,0) cell
        before."""
        slab1 = copy.deepcopy(self)
        slab2 = copy.deepcopy(slab)
        slab1.collapseCartesianCoordinates()
        slab2.collapseCartesianCoordinates()
        # reorder sublayers by Z to then compare by index
        slab1.sublayers.sort(key=lambda sl: sl.cartbotz)
        slab2.sublayers.sort(key=lambda sl: sl.cartbotz)
        ab = self.ucell[:2, :2]
        for (i, sl) in enumerate(slab1.sublayers):
            if (len(sl.atlist) != len(slab2.sublayers[i].atlist)
                    or abs(sl.cartbotz-slab2.sublayers[i].cartbotz) > eps
                    or sl.atlist[0].el != slab2.sublayers[i].atlist[0].el):
                return False
            for at1 in sl.atlist:
                complist = [at1.cartpos[0:2]]
                # if we're close to an edge or corner, also check translations
                for j in range(0, 2):
                    releps = eps / np.linalg.norm(ab[:, j])
                    if abs(at1.pos[j]) < releps:
                        complist.append(at1.cartpos[:2] + ab[:, j])
                    if abs(at1.pos[j]-1) < releps:
                        complist.append(at1.cartpos[:2] - ab[:, j])
                if len(complist) == 3:
                    # coner - add the diagonally opposed one
                    complist.append(complist[1] + complist[2] - complist[0])
                found = False
                for at2 in slab2.sublayers[i].atlist:
                    for p in complist:
                        if np.linalg.norm(p-at2.cartpos[0:2]) < eps:
                            found = True
                            break
                    if found:
                        break
                if not found:
                    return False
        return True

    def revertUnitCell(self, restoreTo=None):
        """If the unit cell in a and b was transformed earlier, restore the
        original form and coordinates. If a 'restoreTo' argument is passed,
        restore only back to the point defined by the argument."""
        if restoreTo is None:
            restoreTo = []
        if len(self.ucell_mod) > 0:
            self.getCartesianCoordinates()
            oplist = self.ucell_mod[len(restoreTo):]
            for op in list(reversed(oplist)):
                if op[0] == 'add':
                    for at in self.atlist:
                        at.cartpos[0:2] -= op[1]
                    self.collapseCartesianCoordinates()
                elif op[0] == 'lmul':
                    self.ucell = np.dot(np.linalg.inv(op[1]), self.ucell)
                    self.collapseCartesianCoordinates()
                elif op[0] == 'rmul':
                    self.ucell = np.dot(self.ucell, np.linalg.inv(op[1]))
                    self.collapseCartesianCoordinates()
            self.ucell_mod = self.ucell_mod[:len(restoreTo)]

    def getMinUnitCell(self, rp, warn_convention=False):
        """Checks whether there is a unit cell (a,b) with a smaller area than
        the current one. If so, returns True and the minimized unit cell, else
        returns False and the current unit cell."""
        eps = rp.SYMMETRY_EPS
        epsz = rp.SYMMETRY_EPS_Z
        abst = self.ucell[:2, :2].T
        # create a testslab: C projected to Z
        ts = copy.deepcopy(self)
        ts.projectCToZ()
        ts.sortByZ()
        ts.createSublayers(epsz)
        # find the lowest occupancy sublayer
        lowocclayer = ts.getLowOccLayer()
        if len(lowocclayer.atlist) < 2:
            return (False, abst)  # nothing to check
        # create list of translation vectors to test
        plist = [at.cartpos[:2] for at in lowocclayer.atlist]
        vlist = [p - plist[0] for p in plist[1:]]
        mincell = abst.copy()
        smaller = False
        compare_to = ts.getCompareCoords(eps)  # for isTranslationSymmetric_2D
        for v in (v for v in vlist if ts.isTranslationSymmetric_2D(
                v, eps, compare_to=compare_to)):
            tcell = np.array([mincell[0], v])
            if (abs(np.linalg.det(tcell)) < (abs(np.linalg.det(mincell)) - eps)
                    and abs(np.linalg.det(tcell)) > eps):
                mincell = tcell
                smaller = True
            else:
                tcell = np.array([mincell[1], v])
                if (abs(np.linalg.det(tcell)) < (abs(np.linalg.det(mincell))
                                                 - eps)
                        and abs(np.linalg.det(tcell)) > eps):
                    mincell = tcell
                    smaller = True
        if not smaller:
            return(False, abst)
        else:
            mincell, _, _ = tl.leedbase.reduceUnitCell(mincell)
            # cosmetic corrections
            if abs(mincell[0, 0]) < eps and abs(mincell[1, 1]) < eps:
                # if matrix is diagonal, swap a and b
                tmp = mincell[1].copy()
                mincell[1] = mincell[0].copy()
                mincell[0] = tmp
            if abs(mincell[1, 0]) < eps and abs(mincell[0, 1]) < eps:
                # if matrix is diagonal, make elements positive
                mincell = abs(mincell)
            # by convention, make the shorter vector the first one
            if np.linalg.norm(mincell[0]) > np.linalg.norm(mincell[1]) + eps:
                if abs(mincell[1, 0]) < eps and abs(mincell[0, 1]) < eps:
                    # if matrix is diagonal, DO NOT make it off-diagonal
                    if warn_convention:
                        logger.warning(
                            "The unit cell orientation does not follow "
                            "standard convention: to keep SUPERLATTICE matrix "
                            "diagonal, the first bulk vector must be larger "
                            "than the second. Consider swapping the unit cell "
                            "vectors.")
                else:
                    mincell = np.dot(np.array([[0, 1], [-1, 0]]), mincell)
            # finally, make sure it's right-handed
            if angle(mincell[0], mincell[1]) < 0:
                mincell = np.dot(np.array([[1, 0], [0, -1]]), mincell)
            return(True, mincell)

    def getBulkRepeat(self, rp):
        """Based on a pre-existing definition of the bulk, tries to identify
        a repeat vector for which the bulk matches the slab above. Returns that
        vector in cartesian coordinates, or None if no match is found."""
        eps = rp.SYMMETRY_EPS
        if len(self.sublayers) == 0:
            self.createSublayers(rp.SYMMETRY_EPS_Z)
        if self.bulkslab is None:
            self.makeBulkSlab(rp)
        if len(self.bulkslab.sublayers) == 0:
            self.bulkslab.createSublayers(rp.SYMMETRY_EPS_Z)
        nsub = len(self.bulkslab.sublayers)
        if len(self.sublayers) < 2*nsub:
            return None
        # nonbulk_subl = self.sublayers[:-nsub]
        z_range = (self.sublayers[-nsub].atlist[0].cartpos[2],
                   self.sublayers[-1].atlist[0].cartpos[2])
        baseLayer = self.sublayers[-1-nsub]
        ori = baseLayer.atlist[0].cartpos  # compare displacements from here
        repeat_vectors = []
        for at in self.sublayers[-1].atlist:
            v = at.cartpos - ori
            if self.isTranslationSymmetric(v, eps, z_periodic=False,
                                           z_range=z_range):
                repeat_vectors.append(-v)
        if len(repeat_vectors) == 0:
            return None
        # pick the shortest repeat vector
        cv = min(repeat_vectors, key=lambda x: np.linalg.norm(x))
        cv[2] = -cv[2]  # leed coordinates to standard
        return cv

    def getMinC(self, rp, z_periodic=True):
        """Checks whether there is a vector c with a smaller length than
        the current one. If so, returns the minimized vector, else returns
        None."""
        eps = rp.SYMMETRY_EPS
        pcands = self.getCandidateLayerPeriod(eps)
        if len(pcands) == 0:
            return None
        ts = copy.deepcopy(self)
        ts.getCartesianCoordinates()
        ts.createSublayers(eps)
        baseLayer = ts.sublayers[0]
        baseInd = ts.sublayers.index(baseLayer)
        nl = len(ts.sublayers)
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
        cFracBase = np.dot(np.linalg.inv(ts.ucell[:2, :2]), repeatC[:2]) % 1.0
        newC = np.append(np.dot(ts.ucell[:2, :2], cFracBase), -repeatC[2])
        for (i, j) in [(0, -1), (-1, 0), (-1, -1)]:
            v = np.dot(ts.ucell[:2, :2], cFracBase + np.array([i, j]))
            if (np.linalg.norm(np.append(v, -repeatC[2]))
                    < np.linalg.norm(newC)):
                newC[:2] = v
        return newC

    def getCandidateLayerPeriod(self, eps=1e-4):
        """For a bulk slab, find at what offsets the sublayers repeat,
        checking only element and number of atoms. Returns a list of integer
        offsets between layer indices that are potentially equivalent."""
        if len(self.sublayers) <= 1:
            return([])
        cl = []     # candidate layers
        h = self.ucell[2, 2]  # cell height; periodicity cannot go beyond h/2
        l0 = self.sublayers[0]
        nl = len(self.sublayers)
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
            for j in range(1, int(np.ceil(len(self.sublayers)/2))):
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

    def makeSupercell(self, transform):
        """Returns a copy of the slab with the unit cell transformed by the
        given integer-valued, (2x2) transformation matrix."""
        if np.any(abs(np.round(transform) - transform) > 1e-6):
            raise ValueError("Slab.makeSupercell: transformation matrix "
                             "contains non-integer elements")
        transform = np.round(transform).astype(int)
        transformSize = int(round(abs(np.linalg.det(transform))))
        ts = copy.deepcopy(self)
        if transformSize > 1:
            transformDiag = [1, 1]
            if np.max(transform[:, 0]) > np.max(transform[:, 1]):
                longSide = 0
            else:
                longSide = 1
            transformDiag[longSide] = np.max(transform)
            while transformSize / transformDiag[longSide] % 1 != 0:
                transformDiag[longSide] -= 1
            transformDiag[1-longSide] = int(transformSize
                                            / transformDiag[longSide])
            cpatlist = ts.atlist[:]
            for at in cpatlist:
                for i in range(0, transformDiag[0]):
                    for j in range(0, transformDiag[1]):
                        if i == j == 0:
                            continue
                        tmpat = at.duplicate() # duplicate saves duplicated atom in slab
                        tmpat.pos[0] += i
                        tmpat.pos[1] += j
        ts.resetAtomOriN()
        ts.getCartesianCoordinates(updateOrigin=True)
        tm = np.identity(3, dtype=float)
        tm[:2, :2] = transform
        ts.ucell = np.transpose(np.dot(tm, np.transpose(ts.ucell)))
        ts.getFractionalCoordinates()
        ts.getCartesianCoordinates(updateOrigin=True)
        return ts

    def changeBulkCell(self, rp, newcell):
        """Takes a unit cell (a,b), calculates a new SUPERLATTICE parameter
        from it and creates a new bulk slab with that cell. If a different
        SUPERLATTICE was defined by the user, outputs an error and returns."""
        # calculate new SUPERLATTICE matrix
        abst = np.transpose(self.ucell[:2, :2])
        newSL = np.dot(abst, np.linalg.inv(newcell))
        if not np.all(abs(newSL - newSL.round()) < 5e-2):
            logger.error(
                "Automatically detected bulk SUPERLATTICE is "
                "not integer-valued: \n"+str(newSL)+"\n"
                "# User-defined SUPERLATTICE will be used instead.")
            rp.setHaltingLevel(2)
            return
        else:
            newSL = newSL.round()
        if rp.superlattice_defined:
            logger.warning(
                "Automatically detected minimum-area bulk unit cell differs "
                "from the cell defined by the SUPERLATTICE parameter. "
                "Consider changing the SUPERLATTICE parameter. Found matrix: "
                "\n" + str(newSL.astype(int)))
            rp.setHaltingLevel(2)
            return
        else:
            rp.SUPERLATTICE = newSL
            self.bulkslab = self.makeBulkSlab(rp)

    def addBulkLayers(self, rp, n=1):
        """Returns a copy of the slab with n bulk units appended at the
        bottom, and a list of the new atoms that were added."""
        ts = copy.deepcopy(self) # temporary slab
        newbulkats = []
        duplicated = []
        zdiff = 0.
        for _ in range(n):
            blayers = [lay for lay in ts.layers if lay.isBulk]
            if isinstance(rp.BULK_REPEAT, np.ndarray):
                bulkc = np.copy(rp.BULK_REPEAT)
                if bulkc[2] < 0:
                    # perhaps vector given from surface to bulk instead of reverse...
                    bulkc = -bulkc
            else:
                cvec = ts.ucell[:, 2]
                if zdiff == 0. and rp.BULK_REPEAT is None:
                    # assume that interlayer vector from bottom non-bulk to top
                    # bulk layer is the same as between bulk units
                    zdiff = (blayers[-1].cartbotz
                             - ts.layers[blayers[0].num-1].cartbotz)
                elif zdiff == 0. and isinstance(rp.BULK_REPEAT,
                                                (float, np.floating)):
                    zdiff = rp.BULK_REPEAT
                bulkc = cvec * zdiff / cvec[2]
            ts.getCartesianCoordinates()
            cfact = (ts.ucell[2, 2] + abs(bulkc[2])) / ts.ucell[2, 2]
            ts.ucell[:, 2] = ts.ucell[:, 2] * cfact
            bulkc[2] = -bulkc[2]
            original_atoms = ts.atlist[:] # all atoms before adding layers
            
            # split bulkc into parts parallel and orthogonal to unit cell c
            # this allows to keep the same ucell and shift correctly the new bulk layers
            c_direction = ts.ucell[:, 2] / np.dot(ts.ucell[:, 2], ts.ucell[:, 2])
            bulkc_project_to_c = np.dot(bulkc, ts.ucell[:, 2]) * c_direction
            bulkc_perp_to_c = bulkc - bulkc_project_to_c
            added_this_loop = []
            for at in original_atoms:
                if at.layer.isBulk and at not in duplicated:
                    new_atom = at.duplicate()
                    newbulkats.append(new_atom)
                    duplicated.append(at)
                    added_this_loop.append(new_atom)
                    new_atom.oriN = len(ts.atlist)

                # old atoms get shifted up along ucell c
                at.cartpos += bulkc_project_to_c
            for at in added_this_loop:
                # new atoms get shifted perpendicular to ucell c
                at.cartpos -= bulkc_perp_to_c
            # TODO: could be done outside loop?
            ts.collapseCartesianCoordinates(updateOrigin=True)
            ts.sortOriginal()
        return ts, newbulkats

    def doubleBulkSlab(self):
        """Returns a copy of the bulk slab which is doubled in thickness."""
        ts = copy.deepcopy(self)
        bulkc = np.transpose(np.copy(ts.ucell))[2]
        bulkc[2] *= -1
        ts.getCartesianCoordinates(updateOrigin=True)
        tmplist = ts.atlist[:]
        for at in tmplist:
            at.duplicate()
            at.cartpos = at.cartpos + bulkc
        ts.ucell[:, 2] *= 2
        ts.collapseCartesianCoordinates(updateOrigin=True)
        return ts

    def detectBulk(self, rp, second_cut_min_spacing=1.2):
        """
        Determine the minimal bulk repeat vector from BULK_LIKE_BELOW.

        Parameters
        ----------
        rp : Rparams
            Run parameters object. Contains BULK_LIKE_BELOW, not modified.
        second_cut_min_spacing : float, optional
            Minimal z spacing in Angstrom between sublayers of the detected
            bulk unit to make an additional cut through the bulk. If no two
            sublayers are below this threshold, the bulk will be treated as
            a single layer.

        Raises
        ------
        RuntimeError
            Raised if no repeat vector is found.

        Returns
        -------
        bulk_repeat : np.array
            The new repeat vector.
        slab_cuts : list of float
            Layer cuts needed for the bulk.

        """
        tsl = copy.deepcopy(self)   # temporary copy to mess up layers in
        rp_dummy = copy.deepcopy(rp)
        rp_dummy.LAYER_CUTS = [rp.BULK_LIKE_BELOW]
        rp_dummy.N_BULK_LAYERS = 1
        rp_dummy.BULK_LIKE_BELOW = 0.
        tsl.createLayers(rp_dummy)
        # create a pseudo-bulkslab
        tsl.bulkslab = tsl.makeBulkSlab(rp_dummy)
        bsl = tsl.bulkslab
        bsl.createSublayers(rp.SYMMETRY_EPS_Z)
        # reduce unit cell in xy
        changecell, mincell = bsl.getMinUnitCell(rp)
        if changecell:
            tsl.changeBulkCell(rp, mincell)
            bsl = tsl.bulkslab
        bsl.getCartesianCoordinates()
        # detect new C vector
        newC = bsl.getMinC(rp, z_periodic=False)
        if newC is None:
            logger.error("Automatic bulk detection failed: Found no bulk "
                         "repeat vector below the specified cutoff.")
            raise RuntimeError("Failed to detect bulk repeat vector.")

        # reduce the bulk slab
        rp_dummy.BULK_REPEAT = -newC
        rp_dummy.SUPERLATTICE = np.eye(2)
        tsl.bulkslab = tsl.makeBulkSlab(rp_dummy)
        bsl = tsl.bulkslab
        bsl.createSublayers(rp.SYMMETRY_EPS_Z)

        # calculate cut plane
        frac_bulk_thickness = abs(newC[2]) / abs(self.ucell[2, 2])
        frac_lowest_pos = min([at.pos[2] for at in self.atlist])
        frac_bulk_onset = (frac_lowest_pos + frac_bulk_thickness
                           - (rp.SYMMETRY_EPS_Z / self.ucell[2, 2]))
        slab_cuts = [(max([at.pos[2] for at in self.atlist
                          if at.pos[2] < frac_bulk_onset])
                     + min([at.pos[2] for at in self.atlist
                            if at.pos[2] > frac_bulk_onset]))
                     / 2]

        # now look for potential second cut
        if len(bsl.sublayers) > 1:
            maxdist = abs(bsl.sublayers[1].cartbotz
                          - bsl.sublayers[0].cartbotz)
            cutlayer = 0
            for i in range(1, len(bsl.sublayers) - 1):
                d = abs(bsl.sublayers[i+1].cartbotz
                        - bsl.sublayers[i].cartbotz)
                if d > maxdist:
                    maxdist = d
                    cutlayer = i
            if maxdist >= second_cut_min_spacing:
                bulkcut_frac_from_lowest = (
                    bsl.sublayers[cutlayer].atlist[0].pos[2]
                    - (maxdist / (2 * abs(bsl.ucell[2, 2])))
                    - min([at.pos[2] for at in bsl.atlist]))
                slab_cuts.append(frac_lowest_pos
                                 + (bulkcut_frac_from_lowest
                                    * bsl.ucell[2, 2] / self.ucell[2, 2]))
        return (-newC, slab_cuts)

    def makeBulkSlab(self, rp):
        """Copies self to create a bulk slab, in which everything apart from the
        bulk layers is deleted. Returns that bulk slab."""
        if  len(self.layers) < 2:
            err_txt = "Less than two layers detected. Check POSCAR and consider modifying LAYER_CUTS."
            logger.error(err_txt)
            raise RuntimeError(err_txt)
            
        # construct bulk slab
        bsl = copy.deepcopy(self)
        bsl.resetSymmetry()
        bsl.atlist = [at for at in bsl.atlist if at.layer.isBulk]
        bsl.layers = [lay for lay in bsl.layers if lay.isBulk]
        bsl.getCartesianCoordinates()
        al = bsl.atlist[:]     # temporary copy
        al.sort(key=lambda atom: atom.pos[2])
        topat = al[-1]
        if type(rp.BULK_REPEAT) == np.ndarray:
            bulkc = np.copy(rp.BULK_REPEAT)
            if bulkc[2] < 0:
                bulkc = -bulkc
        else:
            cvec = bsl.ucell[:, 2]
            zdiff = 0
            if rp.BULK_REPEAT is None:
                # assume that interlayer vector from bottom non-bulk to top
                #  bulk layer is the same as between bulk units
                zdiff = (bsl.layers[-1].cartbotz
                         - self.layers[bsl.layers[0].num-1].cartbotz)
            elif isinstance(rp.BULK_REPEAT, (float, np.floating)):
                zdiff = rp.BULK_REPEAT
            if  abs(zdiff) < 1e-5:
                err_txt = "Unable to detect bulk interlayer vector. Check POSCAR and consider explicitly setting BULK_REPEAT."
                logger.error(err_txt)
                raise RuntimeError(err_txt)   
            bulkc = cvec * zdiff / cvec[2]
        bsl.ucell[:, 2] = bulkc
        # reduce dimensions in xy
        superlattice = np.identity(3, dtype=float)
        superlattice[:2, :2] = rp.SUPERLATTICE.T
        bsl.ucell = np.dot(bsl.ucell, np.linalg.inv(superlattice))
        if (rp.superlattice_defined and np.linalg.norm(bsl.ucell[:2, 0]) >
                np.linalg.norm(bsl.ucell[:2, 1]) + 1e-4):
            logger.warning(
                "The bulk unit cell defined by SUPERLATTICE does "
                "not follow standard convention: the first vector is larger "
                "than the second. Make sure beams are labelled correctly.")
            rp.setHaltingLevel(1)
        bsl.collapseCartesianCoordinates(updateOrigin=True)
        # if self.ucell_mod is not empty, don't drag that into the bulk slab.
        bsl.ucell_mod = []
        # position in c is now random; make topmost bulk atom top
        topc = topat.pos[2]
        for at in bsl.atlist:
            at.pos[2] = (at.pos[2] - topc + 0.9999) % 1.
        # now center around cell midpoint
        cl = [at.pos[2] for at in bsl.atlist]
        midpos = (max(cl)+min(cl))/2
        for at in bsl.atlist:
            at.pos[2] = (at.pos[2] - midpos + 0.5) % 1.
        bsl.getCartesianCoordinates(updateOrigin=True)
        bsl.updateElementCount()   # update the number of atoms per element
        # remove duplicates
        bsl.createSublayers(rp.SYMMETRY_EPS_Z)
        newatlist = []
        for subl in bsl.sublayers:
            i = 0
            while i < len(subl.atlist):
                j = i+1
                while j < len(subl.atlist):
                    if subl.atlist[i].isSameXY(subl.atlist[j].cartpos[:2],
                                               eps=rp.SYMMETRY_EPS):
                        subl.atlist.pop(j)
                    else:
                        j += 1
                i += 1
            newatlist.extend(subl.atlist)
        bsl.atlist = newatlist
        bsl.updateElementCount()   # update number of atoms per element again
        # update the layers. Don't use Slab.createLayers here to keep it
        #   consistent with the slab layers
        for i, layer in enumerate(bsl.layers):
            layer.slab = bsl
            layer.getLayerPos()
            layer.num = i
            layer.atlist = [at for at in layer.atlist if at in bsl.atlist]
        return bsl

    def makeSymBaseSlab(self, rp, transform=None):
        """Copies self to create a symmetry base slab by collapsing to the
        cell defined by rp.SYMMETRY_CELL_TRANSFORM, then removing duplicates.
        Also assigns the duplicateOf variable for all atoms in self.atlist.
        By default, the transformation matrix will be taken from rp, but a
        different matrix can also be passed."""
        ssl = copy.deepcopy(self)
        ssl.resetSymmetry()
        ssl.getCartesianCoordinates()
        # reduce dimensions in xy
        transform3 = np.identity(3, dtype=float)
        if transform is not None:
            transform3[:2, :2] = transform
        else:
            transform3[:2, :2] = rp.SYMMETRY_CELL_TRANSFORM
        ssl.ucell = np.dot(ssl.ucell, np.linalg.inv(np.transpose(transform3)))
        ssl.collapseCartesianCoordinates(updateOrigin=True)
        ssl.ucell_mod = []
        # if self.ucell_mod is not empty, don't drag that into the new slab.
        # remove duplicates
        ssl.createSublayers(rp.SYMMETRY_EPS_Z)
        newatlist = []
        for subl in ssl.sublayers:
            i = 0
            while i < len(subl.atlist):
                j = i+1
                baseat = [a for a in self.atlist
                          if a.oriN == subl.atlist[i].oriN][0]
                while j < len(subl.atlist):
                    if subl.atlist[i].isSameXY(subl.atlist[j].cartpos[:2],
                                               eps=rp.SYMMETRY_EPS):
                        for a in [a for a in self.atlist
                                  if a.oriN == subl.atlist[j].oriN]:
                            a.duplicateOf = baseat
                        subl.atlist.pop(j)
                    else:
                        j += 1
                i += 1
            newatlist.extend(subl.atlist)
        ssl.atlist = newatlist
        ssl.updateElementCount()   # update number of atoms per element again
        # update the layers. Don't use Slab.createLayers here to keep it
        #   consistent with the slab layers
        for i, layer in enumerate(ssl.layers):
            layer.slab = ssl
            layer.getLayerPos()
            layer.num = i
            layer.atlist = [at for at in layer.atlist if at in ssl.atlist]
        return ssl

    def getSurfaceAtoms(self, rp):
        """Checks which atoms are 'at the surface', returns them as a set."""

        _PTL = set(el.lower() for el in tl.leedbase.PERIODIC_TABLE)
        _RADII = tl.leedbase.COVALENT_RADIUS

        atoms = copy.deepcopy(self.atlist)
        # run from top to bottom of slab
        atoms.sort(key=lambda atom: -atom.pos[2])
        covered = set()
        surfats = set()
        for atom in atoms:
            if atom.el in rp.ELEMENT_MIX:
                # radius as weighted average
                totalocc = 0.0
                r = 0.0
                for chemel in rp.ELEMENT_MIX[atom.el]:
                    if chemel.lower() in _PTL:
                        if chemel in atom.site.occ:
                            r += (_RADII[chemel.capitalize()]
                                  * atom.site.occ[chemel])
                            totalocc += atom.site.occ[chemel]
                    else:
                        logger.error(
                            "Error identifying surface atoms: Could "
                            f"not identify {chemel} as a chemical element.")
                        rp.setHaltingLevel(2)
                        return []
                if totalocc == 0:
                    logger.error("Unexpected point encountered in "
                                 "generateSearchInput: GS001")
                else:
                    r /= totalocc
            else:
                if atom.el.lower() in _PTL:
                    r = _RADII[atom.el.capitalize()]
                else:
                    logger.error("Error identifying surface atoms: Could not "
                                 f"identify {atom.el} as a chemical element.")
                    rp.setHaltingLevel(2)
                    return []
            r *= 1.2    # TODO: !!! test if this is enough
            surfats.update(a for a in self.atlist
                           if (a.pos[2] >= atom.pos[2]
                               and a not in covered))
            covered.update(a for a in self.atlist
                           if (a.pos[2] < atom.pos[2]
                               and a.isSameXY(atom.cartpos[:2], eps=r)))
            if len(covered) + len(surfats) >= len(atoms):
                break   # that's all of them
        return surfats

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

    def getNearestNeigbours(self):
        """Returns a list listing the nearest neighbor distance for all atoms
        in the slab taking periodic boundary conditions into account. For this
        calculation, the cell is internally expanded into a supercell."""

        # unit vectors
        a = self.ucell[:, 0]  # vector a
        b = self.ucell[:, 1]  # vector b

        # Compare unit vector lengths and decide based on this how many cells to add around
        # A minimum 3x3 supercell is constructed for nearest neighbor query, but may be exteneded if vector lengths
        # are very different
        max_length = max(np.linalg.norm(a), np.linalg.norm(b))
        i = np.ceil(max_length/np.linalg.norm(a))
        j = np.ceil(max_length/np.linalg.norm(b))

        # Makes supercell minimum size 3x3 original
        transform = np.array([[2*i+1, 0],
                              [0, 2*j+1]])
        supercell = self.makeSupercell(transform)


        atom_coords = [atom.cartpos for atom in supercell.atlist] # Atom coordinates in supercell
        # For NN query use KDTree from scipy.spacial
        tree = KDTree(atom_coords)

        NN_dict = {}  # Dict containing Atom and NN will be returned

        # Now query atoms in center cell for NN distances and save to dict
        for atom in self.atlist:
            coord = atom.cartpos
            coord += (i+1)*a + (j+1)*b  # central cell

            dists, _ = tree.query(coord, k=2)  # second argument irrelevant; would be index of NN atoms (supercell, not original!)
            NN_dict[atom] = dists[1]  # element 0 is distance to atom itself (< 1e-15)

        return NN_dict

    def apply_matrix_transformation(self, trafo_matrix):
        """Applies an orthogonal transformation to the unit cell and all atoms.

        The transformation is described by an orthogonal transfomation matrix
        (O) which is applied to both the unit cell and all atomic positions.

        This transformation is essentially equivalent to a change of basis.
        The transformation of the unit cell (U) is applied as U' = OUO^T.
        Fractional and cartesian atomic positions (v) are transformed as Ov.
        This method can be used to switch indices, rotate the slab, etc..

        Parameters
        ----------
        trafo_matrix : Sequence
                trafo_matrix must be an orthogonal 3-by-3 matrix.
                Contains the transformation matrix (O) describing
                the applied transformation.

        Raises
        ------
        ValueError
            If trafo_matrix is not 3-by-3 or not orthogonal.

        Examples
        --------
        >>> theta = np.pi/2
        >>> rot_mat = [[np.cos(theta), -np.sin(theta), 0],
                       [np.sin(theta),  np.cos(theta), 0],
                       [0, 0, 1]]
        >>> slab.apply_matrix_transformation(rot_mat)

            Applies a rotation by 90 deg around the z axis to the unit cell.
            (In positive direction, i.e. clockwise looking along z)

        >>> swap_b_c = [[1, 0, 0],
                        [0, 0, 1],
                        [0, 1, 0]]
        >>> slab.apply_matrix_transformation(swap_b_c)

            Applies a transformation that switches the unit cell vectors b
            and c in a non-handedness-conserving manner.
        """
        trafo_matrix = np.asarray(trafo_matrix)
        if trafo_matrix.shape != (3, 3):
            raise ValueError(
                "apply_matrix_transformation: not a 3-by-3 matrix"
            )

        if not np.allclose(np.linalg.inv(trafo_matrix), trafo_matrix.T):
            raise ValueError("apply_matrix_transformation: matrix is not "
                             "orthogonal. Consider unsing apply_scaling.")

        self.ucell = trafo_matrix.dot(self.ucell).dot(trafo_matrix.T)
        for atom in self.atlist:
            atom.pos = trafo_matrix.dot(atom.pos)
        self.getCartesianCoordinates(updateOrigin=True)
        self.collapseFractionalCoordinates()

    def apply_scaling(self, scaling):
        """Applies a scaling along the unit cell vectors.

        This can be used to apply an isotropic scaling in all directions
        or to strech/compress along unit cell vectors in order to change
        lattice constants in some direction. To apply other (orthogonal)
        transformations (e.g. rotation, flipping), use
        apply_matrix_transformation.

        Parameters
        ----------
        scaling : Real, Sequence
            If the scaling is given as a number or a list containing one item, an
            isotropic scaling will be applied to the unit cell and atom positions.
            If a list with 3 entries is provided, the scaling will be applied along
            the unit cell vector in the given order.
        
        Raises
        ------
        TypeError
            If scaling is not a real number nor a sequence
        ValueError
            If scaling is a Sequence with length different than 1 or 3,
            or if scaling would make the unit cell singular (i.e., reduce
            the length of a unit vecotr to zero).
        
        Examples
        ----------
        >>> slab.apply_scaling(2)
        
            Streches the unit cell by a factor of two along a, b and c. This doubles
            the lattice constant and increases the volume 8-fold. 
        >>> slab.apply_scaling([1,1,1/3])
        
            Compresses the unit cell by a factor of 3 along c.
        """
        if isinstance(scaling, Real):
            scaling = (scaling,)*3
        elif isinstance(scaling, (np.ndarray, Sequence)):
            if (len(scaling) not in (1, 3)
                    or len(np.shape(scaling)) != 1):
                raise ValueError(
                    "Slab.apply_scaling: only 1-, or 3-sequences allowed"
                    )
            if len(scaling) == 1:
                scaling = (scaling[0],)*3
        else:
            raise TypeError(
                f"Slab.apply_scaling: invalid type {type(scaling)}. "
                f"Expected number or Sequence."
            )

        if any(abs(s) < 1e-5 for s in scaling):
            raise ValueError("Slab.apply_scaling: cannot reduce unit"
                             "vector(s) to zero length")
        scaling_mat = np.diag(scaling)

        # Apply to unit cell (basis). Notice the inverted order, because the
        # unit cell is stored with unit vectors as columns (i.e., a = ucell[:, 0])
        self.ucell = self.ucell.dot(scaling_mat)

        self.getCartesianCoordinates(updateOrigin=True)  # update cartesian values
