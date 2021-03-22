# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 13:47:08 2020

@author: Florian Kraushofer

Functions for writing files relevant to the superpos calculation
"""

import logging
import numpy as np

from viperleed import fortranformat as ff

logger = logging.getLogger("tleedm.files.iosuperpos")


def writeSuperposInput(sl, rp, config, param_name="PARAM",
                       contrin_name="superpos-CONTRIN",
                       for_error=False):
    """Writes a PARAM and a CONTRIN file as input for the Superpos
    calculation. Returns the contents of CONTRIN, as this will also be needed
    as an input stream for the superpos calculation."""
    mnfiles = 0
    occupations = []    # occupation per delta file
    deltanames = []     # names of delta files
    surfaceChecks = []  # is that atom at the surface or not?
    indices = []        # vib and geo indices from search, cast to 1D array
    n_var = 1
    n_conc = 1
    if for_error:
        n_var = max([sp.steps for sp in rp.searchpars
                     if sp.mode in ("geo", "vib")
                     and sp.atom in rp.search_atlist])
        n_conc = max([sp.steps for sp in rp.searchpars if sp.mode == "occ"
                      and sp.atom in rp.search_atlist])
    # determine which atoms are "at the surface"
    surfats = sl.getSurfaceAtoms(rp)
    # make lists to print
    for at in rp.search_atlist:
        sps = [sp for sp in rp.searchpars if sp.atom == at]
        if not for_error:
            occpar = [sp for sp in sps if sp.mode == "occ"][0]  # exactly one
            occind = rp.searchpars.index(occpar)
        totalocc = np.array([0.] * n_conc)
        for el in at.disp_occ:
            if at in surfats:
                surfaceChecks.append(1)
            else:
                surfaceChecks.append(0)
            if not for_error:
                o = [at.disp_occ[el][config[occind]-1]]
            else:
                try:
                    o = [at.disp_occ[el][j] for j in range(n_conc)]
                except IndexError:
                    logger.error("Superpos for error: Inconsistent number of "
                                 "variation steps. Found maximum {} steps,"
                                 " but only {} for {}, element {}."
                                 .format(n_conc, len(at.disp_occ[el]), at, el))
                    rp.setHaltingLevel(2)
                    return ""
            occupations.append(o)
            totalocc += np.array(o)
            pl = [sp for sp in sps if sp.el == el]
            if len(pl) == 0:
                logger.error("No search parameters found for atom {}."
                             "Aborting...".format(at.oriN))
                rp.setHaltingLevel(2)
                return ""
            deltanames.append(pl[0].deltaname)
            mnfiles += 1
            if not for_error:
                if el in at.disp_geo:
                    k = el
                else:
                    k = "all"
                ngeo = len(at.disp_geo[k])
                vibpar = [sp for sp in pl if sp.mode == "vib"]
                if len(vibpar) > 0:
                    vibind = config[rp.searchpars.index(vibpar[0])] - 1
                else:
                    vibind = 0
                geopar = [sp for sp in pl if sp.mode == "geo"]
                if len(geopar) > 0:
                    geoind = config[rp.searchpars.index(geopar[0])] - 1
                else:
                    geoind = 0
                # cast vib and geo indices from 2D to 1D:
                indices.append([(ngeo * vibind) + geoind + 1])
            else:
                if pl[0].steps == 1:
                    indices.append([1] * n_var)
                elif pl[0].steps == n_var:
                    indices.append(list(range(1, n_var + 1)))
                else:
                    logger.error("Superpos for error: Inconsistent number of "
                                 "variation steps. Found maximum {} steps,"
                                 " but only {} for {}."
                                 .format(n_var, pl[0].steps, at))
                    rp.setHaltingLevel(2)
                    return ""
        vp = [sp for sp in sps if sp.el == "vac"]
        if len(vp) != 0 and np.any(totalocc < 1.0):
            occupations.append(list(1 - totalocc))
            deltanames.append(vp[0].deltaname)
            mnfiles += 1
            indices.append([1] * n_var)
            if at in surfats:
                surfaceChecks.append(1)
            else:
                surfaceChecks.append(0)

    # write PARAM
    param = """C  DIMENSIONS MAY BE CHANGED IN PARAMETER-STATEMENT
C
C  MNFILES: (maximum) number of different files containing amplitude changes
C           to be used as input files for the I(E) calculation
C  MNCONCS: (maximum) number of different combinations of occupation probabilities
C           for the lattice sites in the surface - i.e., weight factors multiplied
C           to the delta amplitudes from each file
C  MNT0:    (maximum) number of beams for which I(E) spectra are calculated
C  MNCSTEP: (maximum) number of parameter combinations (geometry x vibrations)
C           per energy tabulated in an individual delta amplitude file
C  MNATOMS: currently inactive - set to 1
C
"""
    param += "      PARAMETER (MNFILES={})\n".format(mnfiles)
    param += "      PARAMETER (MNCONCS={})\n".format(n_conc)
    param += ("      PARAMETER (MNT0={}, MNCSTEP={}, MNATOMS=1)\n"
              .format(len(rp.ivbeams), rp.mncstep))
    try:
        with open(param_name, "w") as wf:
            wf.write(param)
    except Exception:
        logger.error("Failed to write to PARAM file for superpos")
        raise

    # collect contrin output
    i3 = ff.FortranRecordWriter("I3")
    i3x100 = ff.FortranRecordWriter("100I3")
    f74x10 = ff.FortranRecordWriter('10F7.4')
    contrin = (i3.write([mnfiles]) + "  0               no. of files, VarAll: "
               "all possible parameter combinations? (1/0)\n")
    contrin += i3.write([n_conc]) + ("                  "
                                     "number of concentration steps\n")
    occ_trans = map(list, zip(*occupations))  # transpose occupations
    for occ in occ_trans:
        contrin += f74x10.write(occ) + "\n"
    for i in range(0, len(deltanames)):
        contrin += deltanames[i].ljust(15) + i3.write([n_var])
        contrin += (i3.write([surfaceChecks[i]]) + "  1  FILENAME(A15 !!!),"
                    "VARIATIONS,SURFACE?,FORMATTED?\n")
        contrin += i3x100.write(indices[i]) + "\n"
    try:
        with open(contrin_name, "w") as wf:
            wf.write(contrin)
    except Exception:
        logger.error("Failed to write " + contrin_name + " file. Execution "
                     "will continue...")

    return contrin
