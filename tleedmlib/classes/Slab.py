# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Class accumulating atoms and layers, listing their elements and various other
properties. Includes functions for manipulation of those properties.
"""

#from timeit import default_timer as timer

import logging
import numpy as np
import copy
import re
import scipy.spatial as sps
import itertools

from tleedmlib.base import angle, rotMatrix
import tleedmlib as tl
# from tleedmlib import DEFAULT

logger = logging.getLogger("tleedm.slab")

class SymPlane:
    """Candidate plane for a symmetry operation. 'ty' pre-defines a type
    (mirror or glide), 'index2' allows the (1,2) and (2,1) directions if True,
    and collapse moves pos into the (0,0) unit cell if True."""
    def __init__(self, pos, dr, abt, ty="none", index2=False, collapse=True):
        if collapse:  #collapse to (0,0) cell
            self.pos = np.dot(np.transpose(abt),
                         (np.dot(np.linalg.inv(np.transpose(abt)), pos) % 1.0))
        else:
            self.pos = pos
        self.dir = dr/np.linalg.norm(dr)
                            #normalized vector perpendicular to pos = in-plane
        self.type = ty
        self.par = []
        optionlist = [(1,0),(0,1),(1,1),(1,-1)]
        if index2: 
            optionlist.extend([(2,1),(1,2)])
        for (i,j) in optionlist:
            if abs((abs(np.dot(self.dir,(i*abt[0]+j*abt[1])))
                  / (np.linalg.norm(self.dir)
                     * np.linalg.norm(i*abt[0]+j*abt[1])))-1.0) < 0.001:
                self.par = np.array([i,j])

    def isEquivalent(self,pl2,abt,eps=0.001):
        """Checks whether two symmetry planes have the same position and
        direction (including duplicates in next unit cell)"""
        if not np.array_equal(self.par,pl2.par): return False
        complist = [self.pos]
        fpos = np.dot(np.linalg.inv(np.transpose(abt)), self.pos) % 1.0
        # if we're close to an edge or corner, also check translations
        for i in range(0,2):
            releps = eps / np.linalg.norm(abt[i])
            if abs(fpos[i]) < releps:
                complist.append(self.pos+abt[i])
            if abs(fpos[i]-1) < releps:
                complist.append(self.pos-abt[i])
        if len(complist) == 3:  # coner - add the diagonally opposed one
            complist.append(complist[1]+complist[2]-complist[0])

        for p in complist:
            if tl.base.distanceLineThroughPointsFromPoint(pl2.pos,
                                                     pl2.pos+pl2.dir,p) < eps:
                return True
        return False

class Slab:
    """Contains unit cell, element information and atom coordinates. After
    feeding raw data into atpos, call initAtomList to create a list of Atom
    objects."""
    def __init__(self):
        self.ucell = np.array([])         # unit cell
        self.atpos = []	        # list of atom positions as read from POSCAR
                                #   (list of arrays)
        self.oriels = []        # element labels as read from POSCAR
        self.elements = []      # element labels as read from POSCAR, with
                                #   potential modifications after ELEMENT_MIX
        self.chemelem = []      # chemical elements, including from ELEMENT_MIX
        self.nelem = 0          # number of different elements, including from
                                #   ELEMENT_MIX
        self.nperelem = []      # number of atoms per element
        self.lastupdateelmix = None    # keep track of whether ELEMENT_MIX
                                #  changed since updateElements was last called
        self.atlist = []        # list of Atom objects
        self.layers = []        # list of Layer objects, where "layer" is a
                                #   composite of sublayers, as in TensErLEED
        self.sublayers = []     # list of Layer objects, each containing atoms
                                #   of equal element and Z coordinate
        self.sitelist = []      # list of distinct sites as Sitetype elements
        self.uCellMod = []      # stored modifications made to the unit cell;
                                #   each is a tuple of (type, matrix), where
                                #   type is lmul, rmul, or add
        self.uCellOri = np.array([])
        self.topatOriZ = None     # stores the original position of the topmost
                                  #   atom in cartesian coordinates
        self.celltype = "unknown"       # unit cell type as string
        self.planegroup = "unknown"     # symmetry group of the slab, as string
        self.foundplanegroup = "unknown" # highest symmetry found, doesn't get
                                #  modified when user reduces symmetry manually
        self.orisymplane = None      # only stored if the planegroup is
                                #   ambigious as to which unit vector the
                                #   symmetry plane at the origin is parallel to
        self.linklists = []     # list of lists of atoms which are linked by a
                                #   symmetry operation
        self.displists = []     # list of lists of atoms which are displaced
                                #   together

        self.sitesInitialized = False
        self.layersInitialized = False
        self.preprocessed = False    # True if the POSCAR it was read from had
                                     #   the 'Plane group = XY' comment
        self.deltasInitialized = False

        self.bulkslab = None       # Slab object containing only bulk layers
        self.bulkScrews = []       # only assigned to the bulkslab object!
                    # Integer list of rotation orders present in the bulk
        self.bulkGlides = []       # only assigned to the bulkslab object!
                    # List of symplanes present in the bulk


    def fullUpdate(self,rparams):
        """readPOSCAR initializes the slab with information from POSCAR;
        fullUpdate re-initializes the atom list, then uses the information
        from the parameters file to create layers, calculate cartesian
        coordinates (absolute and per layer), and to update elements and
        sites."""
        #self.initAtomList()
        self.collapseFractionalCoordinates()
        self.getCartesianCoordinates()
        if not self.layersInitialized:
            self.createLayers(rparams)
            self.layersInitialized = True
        self.updateElements(rparams)
        if not self.sitesInitialized:
            self.initSites(rparams)
            self.sitesInitialized = True
        if rparams.fileLoaded["VIBROCC"]:
            for at in self.atlist:
                at.initDisp()

    def updateNperElem(self):
        """Updates the nperelem variable based on the atoms that are actually
        in the slab atlist."""
        for i in range(0,len(self.elements)):
            self.nperelem[i] = 0
            for at in self.atlist:
                if at.el == self.elements[i]: self.nperelem[i] += 1

    def initAtomList(self):
        """Creates a list of Atom objects based on the data read previously"""
        n = 0
        for nat in self.nperelem:
            n += nat
        self.atlist = []
        elnum = 0
        atcount = 0
        for i in range(0,n):
            self.atlist.append(tl.Atom(self.elements[elnum],self.atpos[i],
                                       i+1,self))
            atcount += 1
            if atcount == self.nperelem[elnum]:
                elnum += 1
                atcount = 0
        self.getCartesianCoordinates()

    def getCartesianCoordinates(self, updateOrigin=False):
        """Assigns absolute cartesian coordinates to all atoms, with x,y using
        the unit cell (top plane), while z = 0 for the topmost atom and
        positive going down through the slab. If updateOrigin is set True, the 
        cartesian origin relative to the fractional origin will be updated, 
        otherwise it is static."""
        al = self.atlist[:]     #temporary copy
        al.sort(key=lambda atom: atom.pos[2])
        topat = al[-1]
        topcart = np.dot(self.ucell, topat.pos)
        if updateOrigin or self.topatOriZ is None:
            self.topatOriZ = topcart[2]
        for atom in al:
            atom.cartpos = np.dot(self.ucell, atom.pos)
            atom.cartpos[2] = self.topatOriZ - atom.cartpos[2]

    def getFractionalCoordinates(self):
        """Calculates fractional coordinates for all atoms from their
        cartesian coordinates, using the slab unit cell."""
        uci = np.linalg.inv(self.ucell)
        for at in self.atlist:
            tp = np.copy(at.cartpos)
            tp[2] = self.topatOriZ-tp[2]
            at.pos = np.dot(uci,tp)

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

    def createLayers(self,rparams):
        """Creates a list of Layer objects based on the N_BULK_LAYERS and
        LAYER_CUTS parameters in rparams. If layers were already defined,
        overwrite."""
        # first interpret LAYER_CUTS parameter - is a list of strings
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
                highbound = 1.
                val = None
                if (i > 1) and (rparams.LAYER_CUTS[i-1] in ["<", ">"]):
                    try:
                        val = float(rparams.LAYER_CUTS[i-2])
                    except:
                        logger.warning("LAYER_CUTS: Error parsing left-hand "
                                        "boundary for "+ s)
                    if val is not None:
                        if rparams.LAYER_CUTS[i-1] == "<":
                            lowbound = val
                        else:
                            highbound = val
                if i < len(rparams.LAYER_CUTS) - 2 and (rparams.LAYER_CUTS[i+1]
                                                        in ["<", ">"]):
                    try:
                        val = float(rparams.LAYER_CUTS[i+2])
                    except:
                        logger.warning("LAYER_CUTS: Error parsing right-hand "
                                        "boundary for "+ s)
                    if val is not None:
                        if rparams.LAYER_CUTS[i+1] == ">":
                            lowbound = val
                        else:
                            highbound = val
                if 'dc' in s:
                    cutoff *= (self.ucell[2,2] 
                               / np.linalg.norm(self.ucell[:,2]))
                for i in range(1,len(al)):
                    if ((abs(al[i].cartpos[2]-al[i-1].cartpos[2]) > cutoff) 
                            and al[i].pos[2] > lowbound 
                            and al[i].pos[2] < highbound
                            and al[i-1].pos[2] > lowbound 
                            and al[i-1].pos[2] < highbound):
                        ct.append(abs((al[i].pos[2]+al[i-1].pos[2])/2))
            elif s not in ["<", ">"]:
                try:
                    ct.append(float(s))
                except:
                    logger.warning("LAYER_CUTS: Could not parse value: " + s)
                    continue
        ct.sort()
        self.layers = []
        tmplist = self.atlist[:]
        self.sortByZ()
        laynum = 0
        b = True if rparams.N_BULK_LAYERS > 0 else False
        newlayer = tl.Layer(self,0,b)
        self.layers.append(newlayer)
        for atom in self.atlist:
            #only check for new layer if we're not in the top layer already
            if laynum < len(ct):
                if atom.pos[2] > ct[laynum]:
                    #if atom is higher than the next c cutoff, make a new layer
                    laynum += 1
                    b = True if rparams.N_BULK_LAYERS > laynum else False
                    newlayer = tl.Layer(self,laynum,b)
                    self.layers.append(newlayer)
                    check = True    #check for empty layer
                    while check:
                        if laynum >= len(ct):
                            check = False
                        elif atom.pos[2] <= ct[laynum]:
                            check = False
                        else:
                            laynum += 1
                            b = True if rparams.N_BULK_LAYERS>laynum else False
                            newlayer = tl.Layer(self,laynum,b)
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
            if layer.isBulk: self.layers[layer.num+1].isBulk = True
            self.layers.remove(layer)
            del layer
        self.layers.reverse()
        for i, layer in enumerate(self.layers):
            layer.getLayerPos()
            layer.num = i
        self.atlist = tmplist

    def createSublayers(self,eps=0.001):
        """Sorts the atoms in the slab into sublayers, sorted by element and Z
        coordinate."""
        self.sortByZ()
        subl = [] #will be a list of sublayers, using the Layer
                            #  class, where a sublayer is all atoms of the same
                            #  element at the same z position (+- eps)
        for el in self.elements:
            sublists = [[a for a in self.atlist if a.el == el]]
            # first, split at points where two atoms are more than eps apart
            i = 0
            while i < len(sublists):
                brk = False
                if len(sublists[i]) > 1:
                    tmplist = sublists[i][:]
                    for j in range(1,len(tmplist)):
                        if (abs(tmplist[j].cartpos[2] - tmplist[j-1].cartpos[2])
                              > eps):
                            sublists.append(tmplist[:j])
                            sublists.append(tmplist[j:])
                            sublists.pop(i)
                            brk = True
                            break
                    if not brk: i += 1
                else: i += 1
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
                if not brk: i += 1
            #now, create sublayers based on sublists:
            for l in sublists:
                newsl = tl.Layer(self,0,sublayer=True)
                subl.append(newsl)
                newsl.atlist = l
                newsl.cartbotz = l[0].cartpos[2]
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
        
        for (i,sl) in enumerate(self.sublayers): 
            sl.num = i

    def updateElements(self,rp):
        """Updates nelem based on the ELEMENT_MIX parameter, and warns in case 
        of a naming conflict."""
        if self.lastupdateelmix == rp.ELEMENT_MIX:
            return     #don't update if up to date
        # update nelem
        c = 0
        oldels = self.elements[:]
        for i, pel in enumerate(oldels):
            if not pel in rp.ELEMENT_MIX:
                c += 1
            else:
                c += len(rp.ELEMENT_MIX[pel])
                # check for overlapping names:
                for el in rp.ELEMENT_MIX[pel]:
                    if el in oldels:
                        logger.warning('element name '+el+' given in '
                                'ELEMENT_MIX is also an element name in '
                                'POSCAR. It is recommended you rename the '
                                'element in the POSCAR file.')
        self.nelem = c
        self.chemelem = []
        for el in self.elements:
            if el in rp.ELEMENT_MIX:
                self.chemelem.extend(rp.ELEMENT_MIX[el])
            else:
                self.chemelem.append(el)
        self.lastupdateelmix = rp.ELEMENT_MIX

    def updateElementCount(self):
        # update the number of atoms per element
        for i in range(0, len(self.elements)):
            self.nperelem[i] = len([at for at in self.atlist
                                   if at.el == self.elements[i]])
        # Remove elements that are not present any more
        i = 0
        while i < len(self.elements):
            if self.nperelem[i] == 0:  # the element is not present
                self.nperelem.pop(i)
                self.elements.pop(i)
                self.oriels.pop(i)
            else:
                i += 1

    def initSites(self,rparams):
        """Goes through the atom list and supplies them with appropriate
        Sitedef objects, based on the SITE_DEF parameters from the supplied
        Rparams."""
        atlist = self.atlist[:]     #copy to not have any permanent changes
        atlist.sort(key=lambda atom: atom.oriN)
        sl = []
        for el, sitedict in rparams.SITE_DEF.items():
            for sitename, sitelist in sitedict.items():
                newsite = tl.Sitetype(el, sitename)
                sl.append(newsite)
                for i in sitelist:
                    try:
                        if atlist[i-1].el != el:
                            logger.warning('SITE_DEF tries to assign '
                               'atom number '+str(i)+' as '+el+', but '
                               'POSCAR has it as '+atlist[i-1].el+'. Atom '
                               'will be skipped and left as default site '
                               'type!')
                            rparams.setHaltingLevel(1)
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
        self.sitelist = sl

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
                                               #create a list sorted by element
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
                for at in isoLists[i]:
                    sortedlist.append(at)
        self.atlist = sortedlist

    def sortOriginal(self):
        """Sorts atlist by original atom order from POSCAR"""
        self.atlist.sort(key=lambda atom: atom.oriN)

    def projectCToZ(self):
        """makes the c vector of the unit cell perpendicular to the surface,
        changing all atom coordinates to fit the new base"""
        if self.ucell[0,2] != 0.0 or self.ucell[1,2] != 0.0:
            self.getCartesianCoordinates()
            self.ucell[:,2] = np.array([0,0,self.ucell[2,2]])
            self.collapseCartesianCoordinates()     #implicitly also gets new
                           #  fractional coordinates based on the new unit cell

    def restoreOriState(self, keepDisp=False):
        """Resets the atom positions and site vibrational amplitudes to the 
        original state, and stores the deviations as offset instead."""
        for site in self.sitelist:
            siteats = [at for at in self.atlist if at.site == site]
            for el in site.occ:
                for at in siteats:
                    o = site.vibamp[el] - site.oriState.vibamp[el]
                    if el in at.offset_vib:
                        o += at.offset_vib[el]
                    at.offset_vib[el] = o
                    o = site.occ[el] - site.oriState.occ[el]
                    if el in at.offset_occ:
                        o += at.offset_occ[el]
                    at.offset_occ[el] = o
                    if abs(at.offset_vib[el]) < 1e-6:
                        del at.offset_vib[el]
                    if abs(at.offset_occ[el]) < 1e-6:
                        del at.offset_occ[el]
                site.vibamp[el] = site.oriState.vibamp[el]
                site.occ[el] = site.oriState.occ[el]
        uci = np.linalg.inv(self.ucell)
        self.getCartesianCoordinates()
        for at in self.atlist:
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
            at.disp_geo_offset = {"all": [np.array([0.0,0.0,0.0])]}
        self.collapseFractionalCoordinates()
        self.getCartesianCoordinates()
        self.updateLayerCoordinates()
        if keepDisp:
            return 0
        for at in self.atlist:
            at.deltasGenerated = []
            at.initDisp(force=True)
            at.constraints = {1: {}, 2: {}, 3: {}}
        return 0

    def rotate(self, axis, order):
        """Translates the atoms in the slab to have the axis in the origin,
        applies an order-fold rotation matrix, then translates back"""
        #these explicit definitions are likely useless, but sqrts might be
        #  marginally more accurate than sin/cos
        if order == 2:
            m = np.array([[-1,0],[0,-1]])
        elif order == 3:
            m = np.array([[-0.5,-np.sqrt(3)/2],[np.sqrt(3)/2,-0.5]])
        elif order == -3:
            m = np.array([[-0.5,np.sqrt(3)/2],[-np.sqrt(3)/2,-0.5]])
        elif order == 4:
            m = np.array([[0,1],[-1,0]])
        elif order == -4:
            m = np.array([[0,-1],[1,0]])
        elif order == 6:
            m = np.array([[0.5,np.sqrt(3)/2],[-np.sqrt(3)/2,0.5]])
        elif order == -6:
            m = np.array([[0.5,-np.sqrt(3)/2],[np.sqrt(3)/2,0.5]])
        else:
            ang = 2*np.pi/order
            m = np.array([[np.cos(ang),-np.sin(ang)],
                           [np.sin(ang),np.cos(ang)]])
        for at in self.atlist:
            at.cartpos[0:2] -= axis    # translate origin to candidate point
            at.cartpos[0:2] = np.dot(m, at.cartpos[0:2])    # rotation
            at.cartpos[0:2] += axis    # undo translation

    def mirror(self, symplane, glide=False):
        """Translates the atoms in the slab to have the symplane in the
        origin, applies a mirror or glide matrix, then translates back"""
        ang = angle(np.array([1,0]),symplane.dir)
        if symplane.dir[1] > 0: 
            ang *= -1
        rotm = np.array([[np.cos(ang),-np.sin(ang)],
                          [np.sin(ang),np.cos(ang)]])
        rotmirm = np.dot(np.linalg.inv(rotm),
                         np.dot(np.array([[1,0],[0,-1]]),rotm))
                            #rotates to have plane in x direction, mirrors on x
        if glide:
            abt = np.transpose(self.ucell[0:2,0:2])
            glidevec = (symplane.par[0]*abt[0]+symplane.par[1]*abt[1])/2
        for at in self.atlist:
            at.cartpos[0:2] -= symplane.pos     #translate to plane
            at.cartpos[0:2] = np.dot(rotmirm, at.cartpos[0:2])    #apply mirror
            at.cartpos[0:2] += symplane.pos     #translate back
            if glide:
                # if we're testing for glide plane, add the appropriate vector
                at.cartpos[0:2] += glidevec

    def getLowOccLayer(self):
        """Find and returns the lowest occupancy sublayer"""
        minlen = len(self.sublayers[0].atlist)
        lowocclayer = self.sublayers[0]
        for lay in self.sublayers:
            if len(lay.atlist) < minlen:
                lowocclayer = lay
                minlen = len(lay.atlist)
        return lowocclayer

    def isRotationSymmetric(self,axis,order,eps):
        """Evaluates whether the slab is equivalent to itself when rotated
        around the axis with the given rotational order"""
        m = rotMatrix(order)
        ab = self.ucell[0:2,0:2]
        abt = np.transpose(ab)
        releps = np.array([0.0,0.0])
        for j in range(0,2):
            releps[j] = eps / np.linalg.norm(abt[j])
        shiftv = axis.reshape(2,1)
        for sl in self.sublayers:
            coordlist = [at.cartpos[0:2] for at in sl.atlist]
            shiftm = np.tile(shiftv,len(coordlist))
                                     # matrix to shift all coordinates by axis
            oricm = np.array(coordlist)  # original cartesian coordinate matrix
            oripm = np.dot(np.linalg.inv(ab), oricm.transpose()) % 1.0
                            # collapse (relative) coordinates to base unit cell
            oricm = np.dot(ab, oripm).transpose()
                   # original cartesian coordinates collapsed to base unit cell
            tmpcoords = np.copy(oricm).transpose()
                                    # copy of coordinate matrix to be rotated
            tmpcoords -= shiftm
            tmpcoords = np.dot(m,tmpcoords)
            tmpcoords += shiftm
            tmpcoords = np.dot(ab, (np.dot(np.linalg.inv(ab),tmpcoords) % 1.0))
                                       # collapse coordinates to base unit cell
            # for every point in matrix, check whether is equal:
            for (i,p) in enumerate(oripm.transpose()):
                # get extended comparison list for edges/corners:
                addlist = []
                for j in range(0,2):
                    if abs(p[j]) < releps[j]:
                        addlist.append(oricm[i]+abt[j])
                    if abs(p[j]-1) < releps[j]:
                        addlist.append(oricm[i]-abt[j])
                if len(addlist) == 2:
                    # coner - add the diagonally opposed one
                    addlist.append(addlist[0]+addlist[1]-oricm[i])
                for v in addlist:
                    oricm = np.concatenate((oricm,v.reshape(1,2)))
            distances = sps.distance.cdist(tmpcoords.transpose(), oricm,
                                           'euclidean')
            for sublist in distances:
                if min(sublist) > eps:
                    return False
        return True

    def isTranslationSymmetric(self, tv, eps, z_periodic = True):
        """Evaluates whether the slab is equivalent to itself when translated
        along the given cartesian translation vector tv. 2- or 3-dimensional
        translation vectors are accepted."""
        if len(tv) == 2:  # two-dimensional displacement. append zero for z
            tv = np.append(tv, 0.)
        uc = copy.copy(self.ucell)
        uc[:,2] *= -1   # mirror c vector down
        uct = np.transpose(uc)
        releps = np.array([0.,0.,0.])
        for j in range(0,3):
            releps[j] = eps / np.linalg.norm(uct[j])
        shiftv = tv.reshape(3,1)
        # unlike in-plane operations, this one cannot be done sublayer-internal
        coordlist = [at.cartpos for at in self.atlist]
        shiftm = np.tile(shiftv,len(coordlist))
                                 # matrix to shift all coordinates
        oricm = np.array(coordlist)  # original cartesian coordinate matrix
        oricm[:,2] *= -1
        shiftm[2] *= -1
        oripm = np.dot(np.linalg.inv(uc), oricm.transpose()) % 1.0
                        # collapse (relative) coordinates to base unit cell
        oricm = np.dot(uc, oripm).transpose()
               # original cartesian coordinates collapsed to base unit cell
        tmpcoords = np.copy(oricm).transpose()
                                # copy of coordinate matrix to be manipulated
        tmpcoords += shiftm
        if not z_periodic:
            # discard atoms that moved out of cell in z
            tmpcoords = tmpcoords[:, tmpcoords[2] <= 0]
        tmpcoords = np.dot(uc, (np.dot(np.linalg.inv(uc),tmpcoords) % 1.0))
                        # collapse coordinates to base unit cell
        # for every point in matrix, check whether is equal:
        for (i,p) in enumerate(oripm.transpose()):
            # get extended comparison list for edges/corners:
            addlist = []
            for j in range(0,3):
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
                oricm = np.concatenate((oricm,v.reshape(1,3)))
        distances = sps.distance.cdist(tmpcoords.transpose(), oricm,
                                       'euclidean')
        # print(oricm)
        for sublist in distances:
            if min(sublist) > eps:
                return False
        return True

    def isBulkTransformSymmetric(self, matrix, sldisp, eps):
        """Evalues whether the slab is self-equivalent under a given symmetry 
        operation, and subsequent translation by a given number of sublayers"""
        uc = self.ucell
        uct = np.transpose(uc)
        releps = np.array([0.,0.,0.])
        for j in range(0,3):
            releps[j] = eps / np.linalg.norm(uct[j])
        # get translation vectors to check
        transVecs = []
        lowocclayer = self.getLowOccLayer()
        baseInd = self.sublayers.index(lowocclayer)
        ori = lowocclayer.atlist[0].cartpos
        for at in self.sublayers[(baseInd + sldisp) 
                                 % len(self.sublayers)].atlist:
            transVecs.append((at.cartpos - np.dot(matrix, ori)).reshape(3,1))
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
            for (i,p) in enumerate(oripm2.transpose()):
                # get extended comparison list for edges/corners:
                addlist = []
                for j in range(0,3):
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
                    oricm2 = np.concatenate((oricm2,v.reshape(1,3)))
            j = 0
            while j < len(transVecs):
                v = transVecs[j]
                shiftm = np.tile(v,len(coordlist))
                tmpcoords = transcoords + shiftm
                tmpcoords = np.dot(uc, (np.dot(np.linalg.inv(uc),tmpcoords) 
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
        m = np.array([[1,0,0],[0,1,0],[0,0,1]], dtype=float)
        m[:2, :2] = rotMatrix(order)
        return self.isBulkTransformSymmetric(m, sldisp, eps)

    
    def isBulkGlideSymmetric(self, symplane, sldisp, eps):
        """Evaluates whether the bulk has a glide plane along a given 
        direction, i.e. mirror at this direction, then some translation."""
        m = np.array([[1,0,0],[0,1,0],[0,0,1]], dtype=float)
        ang = angle(np.array([1,0]),symplane.dir)
        if symplane.dir[1] > 0: 
            ang *= -1
        rotm = np.array([[np.cos(ang),-np.sin(ang)],
                          [np.sin(ang),np.cos(ang)]])
        m[:2, :2] = np.dot(np.linalg.inv(rotm),
                         np.dot(np.array([[1,0],[0,-1]]),rotm))
        return self.isBulkTransformSymmetric(m, sldisp, eps)

    def isMirrorSymmetric(self, symplane, eps, glide=False):
        """Evaluates whether the slab is equivalent to itself when applying a
        mirror or glide operation at a given plane"""
        ang = angle(np.array([1,0]),symplane.dir)
        if symplane.dir[1] > 0: 
            ang *= -1
        rotm = np.array([[np.cos(ang),-np.sin(ang)],
                          [np.sin(ang),np.cos(ang)]])
        rotmirm = np.dot(np.linalg.inv(rotm),
                         np.dot(np.array([[1,0],[0,-1]]),rotm))
                            #rotates to have plane in x direction, mirrors on x
        ab = self.ucell[0:2,0:2]
        abt = np.transpose(ab)
        releps = np.array([0.0,0.0])
        for j in range(0,2):
            releps[j] = eps / np.linalg.norm(abt[j])
        shiftv = symplane.pos.reshape(2,1)
        if glide:
            glidev = ((symplane.par[0]*abt[0]+symplane.par[1]*abt[1])
                       / 2).reshape(2,1)
        for sl in self.sublayers:
            coordlist = [at.cartpos[0:2] for at in sl.atlist]
            shiftm = np.tile(shiftv,len(coordlist))
                                            # matrix to shift all coordinates
            if glide: glidem = np.tile(glidev,len(coordlist))
                                # matrix to shift all coordinates by glidev
            oricm = np.array(coordlist)  # original cartesian coordinate matrix
            oripm = np.dot(np.linalg.inv(ab), oricm.transpose()) % 1.0
                            # collapse (relative) coordinates to base unit cell
            oricm = np.dot(ab, oripm).transpose()
                   # original cartesian coordinates collapsed to base unit cell
            tmpcoords = np.copy(oricm).transpose()
                                    # copy of coordinate matrix to be rotated
            tmpcoords -= shiftm
            tmpcoords = np.dot(rotmirm,tmpcoords)
            tmpcoords += shiftm
            if glide: tmpcoords += glidem
            tmpcoords = np.dot(ab, (np.dot(np.linalg.inv(ab),tmpcoords) % 1.0))
                                    # collapse coordinates to base unit cell
            # for every point in matrix, check whether is equal:
            for (i,p) in enumerate(oripm.transpose()):
                # get extended comparison list for edges/corners:
                addlist = []
                for j in range(0,2):
                    if abs(p[j]) < releps[j]:
                        addlist.append(oricm[i]+abt[j])
                    if abs(p[j]-1) < releps[j]:
                        addlist.append(oricm[i]-abt[j])
                if len(addlist) == 2:
                    # coner - add the diagonally opposed one
                    addlist.append(addlist[0]+addlist[1]-oricm[i])
                for v in addlist:
                    oricm = np.concatenate((oricm,v.reshape(1,2)))
            distances = sps.distance.cdist(tmpcoords.transpose(), oricm,
                                           'euclidean')
            for sublist in distances:
                if min(sublist) > eps:
                    return False
        return True

    def isEquivalent(self,slab,eps=0.001):
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
        ab = self.ucell[0:2,0:2]
        for (i,sl) in enumerate(slab1.sublayers):
            if (len(sl.atlist) != len(slab2.sublayers[i].atlist)
                  or abs(sl.cartbotz-slab2.sublayers[i].cartbotz) > eps
                  or sl.atlist[0].el != slab2.sublayers[i].atlist[0].el):
                return False
            for at1 in sl.atlist:
                complist = [at1.cartpos[0:2]]
                # if we're close to an edge or corner, also check translations
                for j in range(0,2):
                    releps = eps / np.linalg.norm(ab[:,j])
                    if abs(at1.pos[j]) < releps:
                        complist.append(at1.cartpos[0:2]+ab[:,j])
                    if abs(at1.pos[j]-1) < releps:
                        complist.append(at1.cartpos[0:2]-ab[:,j])
                if len(complist) == 3:
                    # coner - add the diagonally opposed one
                    complist.append(complist[1]+complist[2]-complist[0])
                found = False
                for at2 in slab2.sublayers[i].atlist:
                    for p in complist:
                        if np.linalg.norm(p-at2.cartpos[0:2]) < eps:
                            found = True
                            break
                    if found: break
                if not found:
                    return False
        return True

    def revertUnitCell(self, restoreTo=[]):
        """If the unit cell in a and b was transformed earlier, restore the
        original form and coordinates. If a 'restoreTo' argument is passed,
        restore only back to the point defined by the argument."""
        if len(self.uCellMod) > 0:
            self.getCartesianCoordinates()
            oplist = self.uCellMod[len(restoreTo):]
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
            self.uCellMod = self.uCellMod[:len(restoreTo)]

    def getMinUnitCell(self, rp):
        """Checks whether there is a unit cell (a,b) with a smaller area than
        the current one. If so, returns True and the minimized unit cell, else
        returns False and the current unit cell."""
        eps = rp.SYMMETRY_EPS
        epsz = rp.SYMMETRY_EPS_Z
        abst = np.transpose(self.ucell[0:2,0:2])
        # create a testslab: C projected to Z
        ts = copy.deepcopy(self)
        ts.projectCToZ()
        ts.sortByZ()
        ts.createSublayers(epsz)
        # find the lowest occupancy sublayer
        lowocclayer = ts.getLowOccLayer()
        minlen = len(lowocclayer.atlist)
        if minlen < 2:
            return(False, abst) # nothing to check
        # create list of translation vectors to test
        plist = [at.cartpos[0:2] for at in lowocclayer.atlist]
        vlist = [(p1-p2) for (p1,p2) in itertools.combinations(plist,2)]
        # check & make list of real translation vectors
        tvecs = [v for v in vlist if ts.isTranslationSymmetric(v, eps)]
        if len(tvecs) == 0:
            return(False, abst)
        mincell = abst.copy()
        smaller = False
        for v in tvecs:
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
            mincell, _, _ = tl.leedbase.reduceUnitCell(mincell, eps)
            # cosmetic corrections
            if abs(mincell[0,0]) < eps and abs(mincell[1,1]) < eps:
                # if matrix is diagonal, swap a and b
                tmp = mincell[1].copy()
                mincell[1] = mincell[0].copy()
                mincell[0] = tmp
            if abs(mincell[1,0]) < eps and abs(mincell[0,1]) < eps:
                # if matrix is diagonal, make elements positive
                mincell = abs(mincell)
            # by convention, make the shorter vector the first one
            if np.linalg.norm(mincell[0]) > np.linalg.norm(mincell[1]) + eps:
                mincell = np.dot(np.array([[0,1],[-1,0]]), mincell)
            return(True, mincell)

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
            if not repeatC is None:
                break
        if repeatC is None:
            return None
        # optimize C vector to be close to Z, if possible
        cFracBase = np.dot(np.linalg.inv(ts.ucell[:2,:2]), repeatC[:2]) % 1.0
        newC = np.append(np.dot(ts.ucell[:2, :2], cFracBase), -repeatC[2])
        for (i,j) in [(0,-1),(-1,0),(-1,-1)]:
            v = np.dot(ts.ucell[:2, :2], cFracBase + np.array([i,j]))
            if (np.linalg.norm(np.append(v, -repeatC[2])) 
                                                    < np.linalg.norm(newC)):
                newC[:2] = v
        return newC

    def getCandidateLayerPeriod(self, eps = 1e-4):
        """For a bulk slab, find at what offsets the sublayers repeat,
        checking only element and number of atoms. Returns a list of integer
        offsets between layer indices that are potentially equivalent."""
        if len(self.sublayers) <= 1:
            return([])
        cl = []     # candidate layers
        h = self.ucell[2,2]  # cell height; periodicity cannot go beyond h/2
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
                if (self.sublayers[(j + cl[i]) % nl].atlist[0].el !=
                                            self.sublayers[j].atlist[0].el):
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

    def changeBulkCell(self, rp, newcell):
        """Takes a unit cell (a,b), calculates a new SUPERLATTICE parameter
        from it and creates a new bulk slab with that cell. If a different
        SUPERLATTICE was defined by the user, outputs an error and returns."""
        # calculate new SUPERLATTICE matrix
        abst = np.transpose(self.ucell[:2,:2])
        newSL = np.dot(abst, np.linalg.inv(newcell))
        if not np.all(abs(newSL - newSL.round()) < 1e-3):
            logger.error("Automatically detected bulk SUPERLATTICE is "
                "not integer-valued: "+str(newSL)+"\n"
                "# User-defined SUPERLATTICE will be used instead.")
            rp.setHaltingLevel(2)
            return
        else:
            newSL = newSL.round()
        if rp.superlattice_defined:
            logger.warning("Automatically detected minimum-area bulk "
                "unit cell differs from the cell defined by the "
                "SUPERLATTICE parameter. Consider changing the "
                "SUPERLATTICE parameter. Found matrix: \n"
                + str(newSL.astype(int)))
            rp.setHaltingLevel(2)
            return
        else:
            rp.SUPERLATTICE = newSL
            self.bulkslab = self.makeBulkSlab(rp)

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
        ts.ucell[:,2] *= 2
        ts.collapseCartesianCoordinates(updateOrigin=True)
        return ts

    def makeBulkSlab(self, rp):
        """Copies self to create a bulk slab, in which evething apart from the
        bulk layers is deleted. Returns that bulk slab."""
        # construct bulk slab
        bsl = copy.deepcopy(self)
        bsl.atlist = [at for at in bsl.atlist if at.layer.isBulk]
        bsl.layers = [l for l in bsl.layers if l.isBulk]
        bsl.getCartesianCoordinates()
        al = bsl.atlist[:]     #temporary copy
        al.sort(key=lambda atom: atom.pos[2])
        topat = al[-1]
        if type(rp.BULK_REPEAT) == np.ndarray:
            bulkc = rp.BULK_REPEAT
            if bulkc[2] < 0:
                bulkc = -bulkc
        else:
            cvec = bsl.ucell[:,2]
            if rp.BULK_REPEAT is None:
                # assume that interlayer vector from bottom non-bulk to top
                #  bulk layer is the same as between bulk units
                zdiff = (bsl.layers[-1].cartbotz
                         - self.layers[bsl.layers[0].num-1].cartbotz)
            elif type(rp.BULK_REPEAT) == float:
                zdiff = rp.BULK_REPEAT
            bulkc = cvec * zdiff / cvec[2]
        bsl.ucell[:,2] = bulkc
        # reduce dimensions in xy
        sl = np.array([[0,0,0],[0,0,0],[0,0,1]],dtype=float)
        sl[:2,:2] = np.transpose(rp.SUPERLATTICE)
        bsl.ucell = np.dot(bsl.ucell, np.linalg.inv(sl))
        if (rp.superlattice_defined and np.linalg.norm(bsl.ucell[:2,0]) > 
                                       np.linalg.norm(bsl.ucell[:2,1]) + 1e-4):
            logger.warning("The bulk unit cell defined by SUPERLATTICE does "
                "not follow standard convention: the first vector is larger "
                "than the second. Make sure beams are labelled correctly.")
            rp.setHaltingLevel(1)
        bsl.collapseCartesianCoordinates(updateOrigin=True)
        bsl.uCellMod = [] # if self.uCellMod is not empty, don't drag that into
                          #  the bulk slab.
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
                                               eps = rp.SYMMETRY_EPS):
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

    def getSurfaceAtoms(self, rp):
        """Checks which atoms are 'at the surface', returns them as a list."""
        abst = np.transpose(self.ucell[0:2,0:2])
        testats = copy.deepcopy(self.atlist)
        testats.sort(key=lambda atom: -atom.pos[2])
        covered = []
        surfats = []
        ptl = [el.lower() for el in tl.leedbase.periodic_table]
        for ta in testats:
            if ta.el in rp.ELEMENT_MIX:
                # radius as weighted average
                totalocc = 0.0
                r = 0.0
                for chemel in rp.ELEMENT_MIX[ta.el]:
                    if chemel.lower() in ptl:
                        if chemel in ta.site.occ:
                            r += (tl.leedbase.elementCovalentRadii[
                                   chemel.capitalize()]* ta.site.occ[chemel])
                            totalocc += ta.site.occ[chemel]
                    else:
                        logger.error("Error identifying surface atoms: Could "
                             "not identify "+chemel+" as a chemical element.")
                        rp.setHaltingLevel(2)
                        return []
                if totalocc == 0:
                    logger.error("Unexpected point encountered in "
                        "generateSearchInput: GS001")
                else:
                    r /= totalocc
            else:
                if ta.el.lower() in ptl:
                    r = tl.leedbase.elementCovalentRadii[ta.el.capitalize()]
                else:
                    logger.error("Error identifying surface atoms: Could not "
                         "identify "+ta.el+" as a chemical element.")
                    rp.setHaltingLevel(2)
                    return []
            r *= 1.1    # !!! test if this is enough
            points = [ta.cartpos[:2]]
            for i in range(-1,2):
                for j in range(-1,2):
                    if not (i == 0 and j == 0):
                        points.append(points[0] + i*abst[0] + j*abst[1])
            surfats.extend([a for a in self.atlist if (a.pos[2] >= ta.pos[2]
                                                  and a not in covered and
                                                  a not in surfats)])
            for at in [a for a in self.atlist if (a.pos[2] < ta.pos[2] and
                                             a not in covered)]:
                for p in points:
                    if np.linalg.norm(at.cartpos[:2] - p) < r:
                        covered.append(at)
                        break
            if len(covered)+len(surfats) >= len(testats):
                break   # that's all of them
        return surfats
    
    def getBulk3Dstr(self):
        """Returns a one-line string containing information about the bulk 
        screw axes and glide planes. Only to be used for bulk slabs. Format of 
        the string is 'r(2, 4), m([1,1], [ 1,-1])'. If neither screw axes nor 
        glide planes exist, returns string 'None'."""
        b3ds = ""
        if self.bulkScrews:
            b3ds += "r({})".format(", ".join([str(v) 
                                              for v in self.bulkScrews]))
        if self.bulkGlides:
            if b3ds:
                b3ds += ", "
            b3ds += "m({})".format(", ".join([np.array2string(gp.par, 
                                                              separator=",") 
                                              for gp in self.bulkGlides]))
        if not b3ds:
            return "None"
        return b3ds