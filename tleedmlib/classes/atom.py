# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Class storing position and other properties of individual atoms (to be used 
with Slab, Layer, etc)
"""

import logging
import numpy as np
import copy

# from tleedmlib import DEFAULT

logger = logging.getLogger("tleedm.atom")
      
class Atom:
    """To be used for Slab; each atom has an element and a position"""
    def __init__(self,el,pos,oriN,slab):
        self.el = el        #element type
        self.pos = pos      #atom position
        self.oriN = oriN    #original atom number in the POSCAR
        self.slab = slab     #which slab the atom belongs to
        self.layer = None     #which layer the atom belongs to. None if 
                              #  no layer was assigned yet.
        self.site = None      #the site type that the atom is in, supplied by 
                              #  the SITEDEF parameter, assigned by the 
                              #  updateSites routine of Slab
        self.cartpos = None   #position in cartesian coordinates, with the 
                                #  highest atom as z = 0, z going down
        self.posInLayer = None   #same as cartpos, but from the layer origin
        self.linklist = []   #defines to which other atoms the atom is linked 
                             #  by symmetry
        self.displist = []   #like linklist, but keeps track of which symmetry
                             #  was active when displacement was defined
        self.freedir = 1   #defines whether the atom can be moved or is locked 
                    # by symmetry. 0: no movements, 1: completely free, 
                    #  np.array([0 or 1, 0 or 1 or -1]): parallel movement to 
                    #  a, b, or diagonal
        self.symrefm = np.array([[1,0],[0,1]])  #matrix defines how a 
                    # translation of the atom at self.linklist[0] should be 
                    #  transformed to affect this atom.
        self.disp_vib = {"all": [0.0]} # list of displacements: vibration,
        self.disp_geo = {"all": [np.array([0.0,0.0,0.0])]}   #  geometry,
        self.disp_occ = {el: [1.0]}                          #  occupation
        self.disp_geo_offset = {"all": [np.array([0.0,0.0,0.0])]}
            # also formatted as list for convenience; should be one element.
        self.dispInitialized = False
          # disp_ variables get initialized after readVIBROCC by Atom.initDisp
        
        self.deltasGenerated = []   # list of filenames as strings: delta 
                                    #  files found or generated for this atom
        
        self.offset_geo = {}    # offset of position per element, compared to
                                #   self.cartpos
        self.offset_vib = {}    # offset of vibrational amplitude per element, 
                                #   compared to self.site.vib
        self.offset_occ = {}    # offset of occupation per element, compared
                                #   to self.site.occ
        self.constraints = {1: {}, 2: {}, 3: {}}
                                # parameter constraints for restrict.f per 
                                #   element, where 1,2,3 = geo,vib,occ. Can 
                                #   be integer-valued index in disp range or 
                                #   a tuple (atom, element)
        self.oriState = None    # deep copy of self before a search is applied
        self.duplicateOf = None  # if this atom is identical to another one 
                            # by translational symmetry (through 
                            # SYMMETRY_CELL_TRANSFORM or domain supercell 
                            # creation), this points to the atom in the 
                            # base cell.
        
    def __str__(self):
        return ("Atom({} {})".format(self.oriN, self.el))
    
    def storeOriState(self):
        """Stores the initial values from the input files for this atom."""
        if self.oriState is None:
            self.oriState = self.duplicate(addToAtlists=False)

    def copyOriState(self, at):
        """Copies positions and offsets from another atom, which may be from 
        another slab."""
        self.storeOriState()
        self.oriState.pos = copy.copy(at.pos)
        self.oriState.cartpos = copy.copy(at.cartpos)
        self.oriState.offset_geo = copy.copy(at.offset_geo)
        self.oriState.offset_vib = copy.copy(at.offset_vib)
        self.oriState.offset_occ = copy.copy(at.offset_occ)

    def initDisp(self, force=False):
        """"Initializes disp_vib, disp_geo and disp_occ, based on the atoms 
        site. Site needs to be assigned first."""
        if (not self.dispInitialized or force) and self.site is not None:
            self.dispInitialized = True
            self.disp_vib = {"all": [0.0]}
            self.disp_geo = {"all": [np.array([0.0,0.0,0.0])]}
            self.disp_occ = {}
            for k, v in self.site.occ.items():
                if v > 0 or k in self.site.mixedEls:
                    self.disp_occ[k] = [v]
        return 0
  
    def mergeDisp(self, el):
        """Merges the offsets from VIBROCC and DISPLACEMENTS into the 
        displacements lists from DISPLACEMENTS for the given element."""
        self.storeOriState()
        for (d, offsetlist) in [(self.disp_geo, self.offset_geo), 
                                (self.disp_vib, self.offset_vib),
                                (self.disp_occ, self.offset_occ),
                                (self.disp_geo, self.disp_geo_offset)]:
            if offsetlist == self.disp_geo_offset:
                if el not in offsetlist:
                    offset = offsetlist["all"][0]
                else:
                    offset = offsetlist[el][0]
            else:
                if el in offsetlist:
                    offset = offsetlist[el]
                else:
                    continue
            if el not in d:
                if d != self.disp_occ:
                    d[el] = copy.copy(d["all"])
                else:
                    logger.error("Atom "+str(self.oriN)+" has occupation "
                        "offset defined for element "+el+", but element "
                        "was not found in atom occupation list.")
            if el in d:
                d[el] = [v + offset for v in d[el]]
                if offsetlist != self.disp_geo_offset:
                    del offsetlist[el]
            if offsetlist == self.disp_geo_offset:
                self.disp_geo_offset = {"all": [np.array([0.0,0.0,0.0])]}
    
    def clearOffset(self, mode, targetel = "", primary = True, displist = []):
        """
        Reverts the atom's offsets for the given mode and element to the 
        original values (from POSCAR and VIBROCC)
        
        Parameters
        ----------
        mode : integer
            Defines what to displace. 1: geo, 2: vib, 3: occ
        targetel : string, optional
            If passed, assignment is made only for that element, otherwise for 
            all.
        primary : bool, optional
            Defines whether assignment should be passed along to linked atoms. 
            This will call assignDisp for these atoms, with primary=False.
        displist : list of Atom objects, optional
            Passed in secondary assignment to later link parameters (the 
            'linklist' defines how the 'displist' is defined, but can change 
            via the SYM_DELTA parameter).
            
        Returns
        -------
        None

        """
        if self.oriState is None or targetel.lower() == "vac":
            return
        if mode == 1:
            td = self.offset_geo
            od = self.oriState.offset_geo
        elif mode == 2:
            td = self.offset_vib
            od = self.oriState.offset_vib
        elif mode == 3:
            td = self.offset_occ
            od = self.oriState.offset_occ
        else:
            logger.warning("Atom.clearOffset: Unknown key for mode (Atom {})."
                            .format(self.oriN))
            return
        if targetel == "":
            els = list(td.keys())
        else:
            els = [targetel]
        for el in els:
            if el not in od:
                del td[el]
            else:
                td[el] = od[el]
        # assign to atoms in linklist:
        if primary:
            if len(self.displist) == 0:
                self.displist = [self]
                self.slab.displists.append(self.displist)
            for at in [at for at in self.linklist if at != self]:
                at.clearOffset(mode, targetel, primary=False, 
                              displist = self.displist)
            return
        if displist != self.displist and len(self.displist) > 0:
            logger.warning("Atom {} is being linked to different groups "
                    "in the DISPLACEMENTS file. This will interfere with  "
                    "correct parameter linking in the search! Check "
                    "SYM_DELTA settings.".format(self.oriN))
        if self not in displist:
            displist.append(self)
        self.displist = displist
        return

    def assignDisp(self, mode, disprange, targetel = "", primary = True, 
                   displist = []):
        """
        Assigns a list of displacements to this atom, for all or a given 
        element.
        
        Parameters
        ----------
        mode : integer
            Defines what to displace. 1: geo, 2: vib, 3: occ, 4: geo offset
        disprange : 
            The list of displacements. For geometrical offsets, pass a list 
            with only one element.
        targetel : string, optional
            If passed, assignment is made only for that element, otherwise for 
            all.
        primary : bool, optional
            Defines whether the assigned displacement should be passed along 
            to linked atoms. This will call assignDisp for these atoms, with 
            primary=False.
        displist : list of Atom objects, optional
            Passed in secondary assignment to later link parameters (the 
            'linklist' defines how the 'displist' is defined, but can change 
            via the SYM_DELTA parameter).
            
        Returns
        -------
        None

        """
        if targetel.lower() == "vac":
            # don't write stuff for vacancies - they don't get displaced, and
            #   occupation can be determined implicitly
            return
        eps = 1e-5  #tolerance for comparing displacements
        dr = copy.copy(disprange)  # to make sure multiple atoms do not get 
                                   #  the same array object 
        if mode == 1:
            td = self.disp_geo
        elif mode == 2:
            td = self.disp_vib
        elif mode == 3:
            td = self.disp_occ
        elif mode == 4:
            td = self.disp_geo_offset
        else:
            logger.warning("Atom.assignDisp: Unknown key for mode (Atom {})."
                            .format(self.oriN))
            return
        if targetel == "":
            els = list(td.keys())
        else:
            els = [targetel]
        for el in els:
            if mode == 4:  # special case: offset
                if not el in td:
                    td[el] = dr
                else:
                    match = True
                    if (abs(dr[0][2]) > eps and 
                            abs(td[el][0][2] - dr[0][2]) > eps):
                        if abs(td[el][0][2]) < eps:
                            td[el][0][2] = dr[0][2]
                        else:
                            match = False
                    if (np.linalg.norm(dr[0][:2]) > eps and 
                            np.linalg.norm(td[el][0][:2] - dr[0][:2]) > eps):
                        if all([(abs(f) < eps) for f in td[el][0][:2]]):
                            td[el][0][:2] = dr[0][:2]
                        else:
                            match = False
                if not match:
                    logger.warning("Atom.assignDisp: Trying to assign "
                        "offset, but atom "+str(self.oriN)+" already has an "
                        "offset defined. Skipping second assignment.")
                continue
            if el not in td:    # did not assign for this element yet -> OK
                td[el] = dr
                continue
            # check whether every value in td[el] is also in disprange
            match = True
            if el in td:
                if len(td[el]) != len(dr):
                    match = False
                else:
                    for v in td[el]:
                        # check whether all values from td[el] are in disprange
                        found = False
                        for v2 in dr:
                            if np.linalg.norm(v-v2) < eps:
                                found = True
                        if not found:
                            match = False
                            break
                if not match and len(td[el]) > 1:
                    logger.warning("Atom.assignDisp: Trying to assign "
                        "displacement list, but atom "+str(self.oriN)+
                        " already has displacements assigned. Skipping second "
                        "assignment.")
                else:
                    td[el] = dr
        # assign to atoms in linklist:
        if primary:
            if len(self.displist) == 0:
                self.displist = [self]
                self.slab.displists.append(self.displist)
            for at in [at for at in self.linklist if at != self]:
                if mode == 1 or mode == 4:
                    tm = np.identity(3)
                    tm[:2,:2] = np.dot(at.symrefm, np.linalg.inv(self.symrefm))
                    newdr = [np.dot(tm, v) for v in dr]
                else:
                    newdr = dr[:]
                at.assignDisp(mode, newdr, targetel, primary=False, 
                              displist = self.displist)
            return
        if displist != self.displist and len(self.displist) > 0:
            logger.warning("Atom {} is being linked to different groups "
                    "in the DISPLACEMENTS file. This will interfere with  "
                    "correct parameter linking in the search! Check "
                    "SYM_DELTA settings.".format(self.oriN))
        if self not in displist:
            displist.append(self)
        self.displist = displist
        return
    
    def assignConstraint(self, mode, targetel = "", value=None, linkAtEl=None,
                         index = None):
        """
        Assign a displacement constraint to this atom. Can be assigned for all 
        elements of only one. Constraint is either a fixed value, or another 
        (Atom, element) pair to link to.

        Parameters
        ----------
        mode : integer
            Defines what parameter to constrain. 1: geo, 2: vib, 3: occ
        targetel : string, optional
            If undefined, constrain for all elements. Otherwise, only constrain 
            for the given element.
        value : float, optional
            The value to which the given parameter should be constrained. If 
            the value is not in the disprange, assignment will be skipped. 
            Leave at default (None) if 'linkAtEl' or 'index' argument is 
            passed instead.
        linkAtEl : tuple (Atom, float), optional
            The atom and element to which the parameter should be linked. If 
            that atom's displist has a different length, assigment will be 
            skipped. Leave at default (None) if 'value' or 'index' argument 
            is passed instead.
        index:
            The index to which the given parameter should be constrained. 
            Leave at default (None) if 'linkAtEl' or 'value' argument is 
            passed instead.
        Returns
        -------
        None
        
        """
        eps = 1e-5
        pars = len([p for p in [value, linkAtEl, index] if p is not None])
        if pars > 1:
            logger.warning("Atom.assignConstraint: Trying to constrain both "
                            "to a fixed value and to another atom for atom {}"
                            .format(self.oriN))
            return
        if pars == 0:
            logger.warning("Atom.assignConstraint: No value, index or target "
                            "atom passed for atom {}".format(self.oriN))
            return
        if mode == 1:
            td = self.disp_geo
        elif mode == 2:
            td = self.disp_vib
        elif mode == 3:
            td = self.disp_occ
        else:  # offset is not allowed here
            logger.warning("Atom.assignConstraint: Unknown key for mode "
                            "(Atom {}).".format(self.oriN))
            return
        if index is not None or value is not None:
            if targetel == "":
                els = list(td.keys())
            else:
                if targetel in td:
                    els = [targetel]
                elif "all" in td:
                    els = ["all"]
                else:
                    logger.warning("Cannot assign constraint for atom {}: "
                            "Element {} does not have displacements assigned."
                            .format(self.oriN, targetel))
                    return
            for el in els:
                if index:
                    if index > len(td[el]) or index < 1:
                        logger.warning("Cannot assign constraint for atom {}"
                                ", element {}: index {} is out of bounds"
                                .format(self.oriN, el, index))
                    else:
                        self.constraints[mode][el] = index - 1
                    continue
                # else value
                if mode == 1:
                    dirvec = td[el][-1] - td[el][0]  # dir of disp range
                    dirvec = dirvec / np.linalg.norm(dirvec)
                    match = [(np.linalg.norm(v - value*dirvec) < eps) 
                              for v in td[el]]
                else:
                    match = [(abs(v - value) < eps) for v in td[el]]
                if any(match):
                    self.constraints[mode][el] = match.index(True)
                else:
                    logger.warning("Cannot assign constraint for atom {}, "
                            "element {}: value {} is not in displacement list"
                            .format(self.oriN, el, value))
        else: # linkAtEl
            (at, el2) = linkAtEl
            if mode == 1:
                td2 = at.disp_geo
            elif mode == 2:
                td2 = at.disp_vib
            elif mode == 3:
                td2 = at.disp_occ
            if el2 in td2:
                l = len(td2[el2])
            else:
                l = len(list(td2.values())[-1])
            if targetel == "":
                els = list(td.keys())
            else:
                els = [targetel]
            for el in els:
                if len(td[el]) != l:
                    logger.warning("Cannot constrain atom {} to atom {}: "
                                    "displacement list lengths differ."
                                    .format(self.oriN, at.oriN))
                    return
            for el in els:
                self.constraints[mode][el] = linkAtEl
        return

    def duplicate(self, addToAtlists = True):
        """Creates a new instance of this atom, using deepcopy for attributes 
        like position and elements, but without making copies of the slab or 
        layer, instead adding the new atom to the existing objects."""
        newat = Atom(self.el, np.copy(self.pos),len(self.slab.atlist),
                     self.slab)
        if addToAtlists:
            self.slab.atlist.append(newat)
            if self.layer is not None: 
                self.layer.atlist.append(newat)
                newat.layer = self.layer
            self.slab.nperelem[self.slab.elements.index(self.el)] += 1
        newat.duplicateOf = self
        newat.site = self.site
        newat.dispInitialized = True
        newat.disp_vib = self.disp_vib
        newat.disp_geo = self.disp_geo
        newat.disp_occ = self.disp_occ
        newat.cartpos = np.copy(self.cartpos)
        return newat

    def isSameXY(self,pos,eps=0.001):
        """Checks whether the atom is at the given x/y coordinates 
        (+- epsilon), taking shifts by one unit vector into account if the 
        atom is close to an edge or corner."""
        abt = np.transpose(self.slab.ucell[0:2,0:2])
        complist = [self.cartpos[0:2]]
        # if we're close to an edge or corner, also check translations:
        for j in range(0,2):
            releps = eps / np.linalg.norm(abt[j])
            if abs(self.pos[j]) < releps:
                complist.append(self.cartpos[0:2]+abt[j])
            if abs(self.pos[j]-1) < releps:
                complist.append(self.cartpos[0:2]-abt[j])
        if len(complist) == 3:
            complist.append(complist[1]+complist[2]-complist[0])    
                                    # coner - add the diagonally opposed one
        for p in complist:
            if np.linalg.norm(p-pos) < eps: 
                return True
        return False
