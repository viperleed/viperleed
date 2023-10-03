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

from viperleed.tleedmlib.base import add_edges_and_corners

logger = logging.getLogger("tleedm.atom")


class Atom:
    """To be used for Slab; each atom has an element and a position.

    Attributes
    ----------
    el : str
        Element type
    pos : numpy.array
        Atom position as fractional coordinate
    oriN : int
        Original atom number in the POSCAR
    slab : Slab
        The Slab that the atom belongs to.
    layer : Layer
        Layer object that the atom belongs to
    site : Sitetype
        Site type of the atom, containing vibrational amplitude and occupation.
        Supplied by the SITEDEF parameter, assigned by the updateSites routine
        of Slab
    cartpos : numpy.array
        Position in cartesian coordinates, with the highest atom as z = 0,
        positive z pointing into the surface
    linklist : list of Atom
        Defines to which other atoms the atom is linked
    displist : list of Atom
        Like linklist, but keeps track of which symmetry was active when
        displacement was defined
    freedir : int
        Defines whether the atom can be moved or is locked by symmetry.
        0: no movements, 1: completely free,
        np.array([0|1, 0|1|-1]): parallel movement to a, b, or diagonal
    symrefm : numpy.array
        Defines how a translation of the atom at self.linklist[0] should be
        transformed to affect this atom.
    disp_vib, disp_geo, disp_occ : dict
        Keys are elements, values are lists of vibration/geometry/occupation
        offsets
    disp_geo_offset : dict
        Keys are elements; values are also formatted as lists for convenience,
        but should be one element long.
    disp_center_index : dict of dict
        Which index in the displacement range corresponds to "no change"
    dispInitialized : bool
        disp_ variables get initialized after readVIBROCC by Atom.initDisp
    deltasGenerated : list of str
        Filenames of delta files generated or found for this atom
    offset_geo, offset_vib, offset_occ : dict
        Offsets from self.cartpos, self.site.vib, self.site.occ per element
    constraints : dict
        Parameter constraints for restrict.f per element, where
        1,2,3 = geo,vib,occ. Can be integer-valued index in disp range or
        a tuple (atom, element)
    oriState : Atom
        Deep copy of self before a search is applied
    duplicateOf : Atom
        If this atom is identical to another one by translational symmetry
        (through SYMMETRY_CELL_TRANSFORM or domain supercell creation),
        this points to the atom in the base cell.
    """

    def __init__(self, el, pos, oriN, slab):
        self.el = el
        self.pos = pos
        self.oriN = oriN
        self.slab = slab
        self.layer = None
        self.site = None
        self.cartpos = None
        self.linklist = []
        self.displist = []
        self.freedir = 1
        self.symrefm = np.identity(2)
        self.disp_vib = {"all": [0.]}
        self.disp_geo = {"all": [np.zeros(3)]}
        self.disp_occ = {el: [1.]}
        self.disp_labels = {
            "geo": "",
            "vib": "",
            "occ": "",
        }
        self.disp_lin_steps = {
            "geo": [],
            "vib": [],
            "occ": [],
        }
        self.disp_geo_offset = {"all": [np.zeros(3)]}
        self.disp_center_index = {"vib": {"all": 0},
                                  "geo": {"all": 0},
                                  "occ": {el: 0}}
        self.dispInitialized = False
        self.deltasGenerated = []
        self.offset_geo = {}
        self.offset_vib = {}
        self.offset_occ = {}
        self.constraints = {1: {}, 2: {}, 3: {}}
        self.oriState = None
        self.duplicateOf = None

    def __str__(self):
        return ("Atom({} {})".format(self.oriN, self.el))

    def storeOriState(self):
        """Stores the initial values from the input files for this atom."""
        if self.oriState is None:
            self.oriState = self.duplicate(addToAtlists=False)

    def copyOriState(self, at):
        """
        Copies positions and offsets from another atom, which may be from
        another slab.
        """
        self.storeOriState()
        self.oriState.pos = copy.copy(at.pos)
        self.oriState.cartpos = copy.copy(at.cartpos)
        self.oriState.offset_geo = copy.copy(at.offset_geo)
        self.oriState.offset_vib = copy.copy(at.offset_vib)
        self.oriState.offset_occ = copy.copy(at.offset_occ)

    def initDisp(self, force=False):
        """
        Initializes disp_vib, disp_geo and disp_occ, based on the atoms site.
        Site needs to be assigned first.
        """
        if (not self.dispInitialized or force) and self.site is not None:
            self.dispInitialized = True
            self.disp_vib = {"all": [0.]}
            self.disp_geo = {"all": [np.zeros(3)]}
            self.disp_occ = {}
            self.disp_center_index = {"vib": {"all": 0},
                                      "geo": {"all": 0},
                                      "occ": {}}
            for k, v in self.site.occ.items():
                if v > 0 or k in self.site.mixedEls:
                    self.disp_occ[k] = [v]
                    self.disp_center_index["occ"][k] = 0
        return None

    def mergeDisp(self, el):
        """
        Merges the offsets from VIBROCC and DISPLACEMENTS into the
        displacements lists from DISPLACEMENTS for the given element.

        For vibrational and occupational offsets a consistency check is
        performed. The offset lists will be emptied.
        """
        self.storeOriState()
        # geometric offsets from DISPLACEMENTS
        geo_d_offset = self.disp_geo_offset.get(el,
                                          self.disp_geo_offset['all'])[0]
        if el not in self.disp_geo:
            self.disp_geo[el] = copy.copy(list(self.disp_geo['all']))
        self.disp_geo[el] = list(self.disp_geo[el] + geo_d_offset)
        self.disp_geo_offset = {"all": [np.zeros(3)]}

        # geometric offsets from VIBROCC
        if el in self.offset_geo:
            geo_offset = self.offset_geo[el]
            self.disp_geo[el] = [geo_step + geo_offset for geo_step in self.disp_geo[el]]
            del self.offset_geo[el]

        # vibrational offsets from VIBROCC
        if el not in self.disp_vib:
            self.disp_vib[el] = copy.copy(self.disp_vib['all'])
        if el in self.offset_vib:
            vib_offset = self.offset_vib[el]
            final_vib_steps = [vib_step + vib_offset for vib_step in self.disp_vib[el]]
            if any(np.array(final_vib_steps) + self.site.vibamp[el] < 0):
                logger.error(f"Vibrational offset for {self} defined in "
                             "VIBROCC would result in negative vibrational "
                             "amplitude. Offset will be ignored.")
            else:
                self.disp_vib[el] = final_vib_steps
            del self.offset_vib[el]

        # vibrational offsets from VIBROCC
        if el in self.offset_occ:
            occ_offset = self.offset_occ[el]
            final_occ_steps = [occ_step + occ_offset for occ_step in self.disp_occ[el]]
            if any(np.array(final_occ_steps) < 0) or any(np.array(final_occ_steps) > 1):
                logger.error(f"Occupational offset for {self} defined in "
                            "VIBROCC would result in unphysical concentration"
                            "(occupation <0 or >1). Offset will be ignored.")
            else:
                self.disp_occ[el] = final_occ_steps
            del self.offset_occ[el]


    def clearOffset(self, mode, targetel="", primary=True, displist=[]):
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
                               displist=self.displist)
            return
        if displist != self.displist and len(self.displist) > 0:
            logger.warning(
                "{} is being linked to different groups in the DISPLACEMENTS "
                "file. This will interfere with correct parameter linking in "
                "the search! Check SYM_DELTA settings.".format(self))
        if self not in displist:
            displist.append(self)
        self.displist = displist
        return

    def assignDisp(self, mode, disprange, targetel="", primary=True,
                   displist=[], disp_label="", disp_lin_steps = []):
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
        eps = 1e-5  # tolerance for comparing displacements
        dr = copy.copy(disprange)
        # to make sure multiple atoms do not get the same list object
        if mode == 1:
            td = self.disp_geo
            self.disp_labels["geo"] = disp_label
            self.disp_lin_steps["geo"] = disp_lin_steps
        elif mode == 2:
            td = self.disp_vib
            self.disp_labels["vib"] = "N/A"  # direction not applicable for vib
            self.disp_lin_steps["vib"] = disp_lin_steps
        elif mode == 3:
            td = self.disp_occ
            self.disp_labels["occ"] = "N/A"  # direction not applicabel for occ
            self.disp_lin_steps["occ"] = disp_lin_steps
        elif mode == 4:
            td = self.disp_geo_offset
        else:
            logger.warning("Atom.assignDisp: Unknown key for mode ({})."
                           .format(self))
            return
        if targetel == "":
            els = list(td.keys())
        else:
            els = [targetel]
        for el in els:
            if mode == 4:  # special case: offset
                if el not in td:
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
                    logger.warning(
                        "Atom.assignDisp: Trying to assign offset, but atom "
                        "{} already has an offset defined. Skipping second "
                        "assignment.".format(self))
                continue
            if (el not in td 
                    or (len(td[el]) == 1 and np.linalg.norm(td[el]) < 1e-5)
                    or (mode == 3 and len(td[el]) == 1)):
                # did not assign for this element yet -> OK, store
                td[el] = dr
                # also store center
                if mode != 3:
                    n = [np.linalg.norm(v) for v in dr]
                else:
                    n = [abs(v - self.site.occ[el]) for v in dr]
                smode = {1: "geo", 2: "vib", 3: "occ"}
                self.disp_center_index[smode[mode]][el] = n.index(min(n))
                continue
            # is warning required: every value in td[el] also in disprange?
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
                                break
                        if not found:
                            match = False
                            break
                if not match and len(td[el]) > 1:
                    logger.warning(
                        "Atom.assignDisp: Trying to assign displacement list, "
                        "but atom {} already has displacements assigned. "
                        "Skipping second assignment.".format(self))
                    return    # also don't assign to other atoms
                # in case of linking, base assignments on current dr
                dr = td[el][:]
        # assign to atoms in linklist:
        if primary:
            if len(self.displist) == 0:
                self.displist = [self]
                self.slab.displists.append(self.displist)
            for at in [at for at in self.linklist if at != self]:
                if mode == 1 or mode == 4:
                    tm = np.identity(3)
                    tm[:2, :2] = np.dot(at.symrefm,
                                        np.linalg.inv(self.symrefm))
                    newdr = [np.dot(tm, v) for v in dr]
                else:
                    newdr = dr[:]
                at.assignDisp(mode, newdr, targetel, primary=False,
                              displist=self.displist,
                              disp_label=disp_label,
                              disp_lin_steps=disp_lin_steps)
            return
        if displist != self.displist and len(self.displist) > 0:
            logger.warning(
                "{} is being linked to different groups in the DISPLACEMENTS "
                "file. This will interfere with correct parameter linking in "
                "the search! Check SYM_DELTA settings.".format(self))
        if self not in displist:
            displist.append(self)
        self.displist = displist
        return

    def assignConstraint(self, mode, targetel="", value=None, linkAtEl=None,
                         index=None):
        """
        Assign a displacement constraint to this atom. Can be assigned for all
        elements or only one. Constraint is either a fixed value, or another
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
                           "to a fixed value and to another atom for {}"
                           .format(self))
            return
        if pars == 0:
            logger.warning("Atom.assignConstraint: No value, index or target "
                           "atom passed for {}".format(self))
            return
        if mode == 1:
            td = self.disp_geo
        elif mode == 2:
            td = self.disp_vib
        elif mode == 3:
            td = self.disp_occ
        else:  # offset is not allowed here
            logger.warning("Atom.assignConstraint: Unknown key for mode "
                           "({}).".format(self))
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
                    logger.warning(
                        "Cannot assign constraint for {}: Element {} does not "
                        "have displacements assigned.".format(self, targetel))
                    return
            for el in els:
                if index:
                    if index > len(td[el]) or index < 1:
                        logger.warning(
                            "Cannot assign constraint for {}, element {}: "
                            "index {} is out of bounds".format(self, el,
                                                               index))
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
                    logger.warning(
                        "Cannot assign constraint for {}, element {}: value "
                        "{} is not in displacement list".format(self, el,
                                                                value))
        else:  # linkAtEl
            (at, el2) = linkAtEl
            if mode == 1:
                td2 = at.disp_geo
            elif mode == 2:
                td2 = at.disp_vib
            elif mode == 3:
                td2 = at.disp_occ
            if el2 in td2:
                listlen = len(td2[el2])
            else:
                listlen = len(list(td2.values())[-1])
            if targetel == "":
                els = list(td.keys())
            else:
                els = [targetel]
            for el in els:
                if len(td[el]) != listlen:
                    logger.warning(
                        "Cannot constrain {} to {}: displacement list lengths "
                        "differ.".format(self, at))
                    return
            for el in els:
                self.constraints[mode][el] = linkAtEl
        return

    def duplicate(self, addToAtlists=True):
        """
        Creates a new instance of this atom, using deepcopy for attributes
        like position and elements, but without making copies of the slab or
        layer, instead adding the new atom to the existing objects.

        Parameters
        ----------
        addToAtlists : bool, optional
            Set to False to not add the atom to atom lists in the existing
            Slab and Layer objects.

        Returns
        -------
        newat : Atom
            The duplicate atom that was created.

        """
        newat = Atom(self.el, np.copy(self.pos), self.slab.n_atoms, self.slab)
        if addToAtlists:
            self.slab.atlist.append(newat)
            if self.layer is not None:
                self.layer.atlist.append(newat)
                newat.layer = self.layer
            self.slab.n_per_elem[self.el] += 1
        newat.duplicateOf = self
        newat.site = self.site
        newat.dispInitialized = True
        newat.disp_vib = self.disp_vib
        newat.disp_geo = self.disp_geo
        newat.disp_occ = self.disp_occ
        newat.cartpos = np.copy(self.cartpos)
        return newat

    def isSameXY(self, pos, eps=1e-3):
        """
        Checks whether the atom is at the given x/y coordinates
        (+- epsilon), taking shifts by one unit vector into account if the
        atom is close to an edge or corner.

        Parameters
        ----------
        pos : numpy.array
            Cartesian xy coordinates to check against the position of this
            atom.
        eps : float, optional
            The precision to which positions are expected to match. The default
            is 1e-3.

        Returns
        -------
        bool
            True if positions match, else False.

        """
        abt = self.slab.ab_cell.T
        releps = eps / np.linalg.norm(abt, axis=1)
        complist, _ = add_edges_and_corners([self.cartpos[:2]], (self.pos,),
                                            releps, abt)
        return any(np.linalg.norm(pos - complist, axis=1) < eps)
