# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 10:28:52 2021

@author: Florian Kraushofer

TensErLEED Manager section Error calculation
"""

import logging
import copy
import numpy as np

import viperleed.tleedmlib as tl
import viperleed.tleedmlib.files.ioerrorcalc as io

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
        d = {}
        at = atoms[0]   # store displacements for this one; main element
        if mode == "occ":
            pass  # !!! TODO - not sure what that looks like...
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
    exclude_mode = {"geo": "vib", "vib": "geo"}
    seg_info = {"geo": "geometrical", "vib": "vibrational"}
    for mode in "geo", "vib":
        sl.restoreOriState()  # reset positions, store any changes as offsets
        # read DISPLACEMENTS block - ONLY geo OR vib
        deltas_required = tl.files.displacements.readDISPLACEMENTS_block(
            rp, sl, rp.disp_blocks[rp.search_index],
            exclude_mode=exclude_mode[mode])
        rp.disp_block_read = True  # to prevent deltas segment from re-reading
        if not deltas_required:
            continue   # !!! how to deal with occupation?
        logger.info("\nStarting error calculations for "
                    + seg_info[mode] + " displacements.")
        # run delta calculations
        tl.sections.deltas(sl, rp)
        sl.deltas_initialized = True
        if rp.STOP:     # since this may stop the deltas, also check here
            return
        rp.generateSearchPars(sl)
        # group atoms under variation by linking
        atom_groups = [[at] for at in rp.search_atlist]
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
            rp.search_atlist = ag  # limit variations to only that group
            logger.info("Now calculating " + seg_info[mode] + " errors for "
                        "atom group: " + ", ".join(str(at) for at in ag))
            logger.info("Running superpos...")
            tl.sections.superpos(sl, rp, for_error=True)
            logger.info("Starting R-factor calculation...")
            rfaclist = tl.sections.rfactor(sl, rp, index=12, for_error=True)
            logger.info("Finished with " + seg_info[mode] + " errors for "
                        "atom group: " + ", ".join(str(at) for at in ag)
                        + "\n")
            errors.append(R_Error(ag, mode, rfaclist))
    varR = np.sqrt(8*np.abs(rp.V0_IMAG) / rp.total_energy_range())
    io.write_errors_csv(errors)
    # !!! PLOT
    return
