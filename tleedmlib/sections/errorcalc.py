# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 10:28:52 2021

@author: Florian Kraushofer

TensErLEED Manager section Error calculation
"""

import logging
import copy
import numpy as np
import os

import viperleed.tleedmlib as tl
import viperleed.tleedmlib.files.ioerrorcalc as tl_io

logger = logging.getLogger("tleedm.error")


class R_Error():
    """
    Data structure for storing errors after they have been calculated.
    Stores the displacements that led to the R-factors even if they are
    modified later. Displacements are stored for one atom only, the other are
    linked.
    """

    def __init__(self, atoms, mode, rfacs):
        self.atoms = atoms  # atoms that have been varied together
        self.mode = mode    # vib, geo, or occ
        self.rfacs = rfacs  # the r-factors from the variations
        self.displacements = []   # displacements of atoms[0]
        self.main_element = ""    # element occupation displayed in output
        d = {}
        at = atoms[0]   # store displacements for this one; main element
        if mode == "occ":
            if atoms[0].el in atoms[0].disp_occ:
                self.main_element = atoms[0].el  # main site element
            else:
                # element with highest occupation in refcalc
                self.main_element = max(atoms[0].site.occ,
                                        key=lambda k: atoms[0].site.occ[k])
            self.displacements = copy.deepcopy(atoms[0].disp_occ[
                                                        self.main_element])
        else:
            if mode == "geo":
                d = atoms[0].disp_geo
            elif mode == "vib":
                d = atoms[0].disp_vib
        if len(d) > 0:
            if at.el in d:
                k = at.el
            else:
                k = "all"
            self.displacements = copy.deepcopy(d[k])


def errorcalc(sl, rp):
    """Runs 1D error calculations based on DISPLACEMENTS."""

    if rp.domainParams:
        # !!! Should be straightforward. Just take the beams from other domains
        # as constant (from refcalc), vary for the one domain as always
        logger.error("Error calculation not implemented for multiple domains.")
        return

    if rp.best_v0r is None:
        logger.info("Error calculation called without a stored inner "
                    "potential shift. Running R-factor calculation "
                    "from refcalc-fd.out to determine inner potential shift.")
        tl.sections.rfactor(sl, rp, index=11)
        logger.info("Finished R-factor pre-run, now starting error "
                    "calculation.\n")

    errors = []
    if rp.search_index > 0 and rp.search_index >= len(rp.disp_blocks):
        # running after one or more searches; read the last DISPLACEMENTS block
        rp.search_index -= 1
    seg_info = {"geo": "geometrical", "vib": "vibrational",
                "occ": "occupation"}
    for mode in "geo", "vib", "occ":
        sl.restoreOriState()  # reset positions, store any changes as offsets
        # read DISPLACEMENTS block - ONLY geo OR vib
        deltas_required = tl.files.displacements.readDISPLACEMENTS_block(
            rp, sl, rp.disp_blocks[rp.search_index],
            only_mode=mode)
        rp.disp_block_read = True  # to prevent deltas segment from re-reading
        if not deltas_required:
            continue
        logger.info("\nStarting error calculations for "
                    + seg_info[mode] + " displacements.")
        # run delta calculations
        tl.sections.deltas(sl, rp)
        sl.deltas_initialized = True
        if rp.STOP:     # since this may stop the deltas, also check here
            return
        rp.generateSearchPars(sl)
        # group atoms under variation by linking
        atom_groups = [[at] for at in rp.search_atlist
                       if [sp for sp in rp.searchpars
                           if sp.atom == at and sp.steps > 1]]
        i = 0
        while i < len(atom_groups):
            found = None
            for at in atom_groups[i]:
                sps = [sp for sp in rp.searchpars if
                       sp.atom == at and sp.mode == mode and sp.el != "vac"
                       and ((isinstance(sp.linkedTo, tl.SearchPar)
                             and sp.linkedTo.atom not in atom_groups[i]) or
                            (isinstance(sp.restrictTo, tl.SearchPar)
                             and sp.restrictTo.atom not in atom_groups[i]))]
                if not sps:
                    continue
                if isinstance(sps[0].linkedTo, tl.SearchPar):
                    found = [ag for ag in atom_groups
                             if sps[0].linkedTo.atom in ag][0]
                else:
                    found = [ag for ag in atom_groups
                             if sps[0].restrictTo.atom in ag][0]
                if found:
                    break
            if found:  # merge the two groups
                found.extend(atom_groups[i])
                atom_groups.pop(i)
            else:
                i += 1
        # run superpos and rfactor:
        for ag in atom_groups:
            # rp.search_atlist = ag  # limit variations to only that group
            only_vary = [sp for sp in rp.searchpars if sp.atom in ag]
            logger.info("\nNow calculating " + seg_info[mode] + " errors for "
                        "atom group: " + ", ".join(str(at) for at in ag))
            logger.info("Running superpos...")
            tl.sections.superpos(sl, rp, for_error=True, only_vary=only_vary)
            if rp.halt >= rp.HALTING:
                return
            if os.path.isfile("ROUTSHORT"):
                os.remove("ROUTSHORT")
            logger.info("Starting R-factor calculation...")
            rfaclist = tl.sections.rfactor(sl, rp, index=12, for_error=True,
                                           only_vary=only_vary)
            logger.info("Finished with " + seg_info[mode] + " errors for "
                        "atom group: " + ", ".join(str(at) for at in ag))
            errors.append(R_Error(ag, mode, rfaclist))
            if rp.halt >= rp.HALTING:
                return
    if len(errors) == 0:
        logger.info("Error calculation: Returning with no output.")
        return
    tl_io.write_errors_csv(errors)
    tl_io.write_errors_pdf(errors, v0i=rp.V0_IMAG, energy_range=rp.total_energy_range())
    return
