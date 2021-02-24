# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 13:47:08 2020

@author: Florian Kraushofer

Functions for writing files relevant to the superpos calculation
"""

import logging
from viperleed import fortranformat as ff

logger = logging.getLogger("tleedm.files.iosuperpos")

def writeSuperposInput(sl, rp, config, param_name = "PARAM",
                                       contrin_name = "superpos-CONTRIN"):
    """Writes a PARAM and a CONTRIN file as input for the Superpos 
    calculation. Returns the contents of CONTRIN, as this will also be needed 
    as an input stream for the superpos calculation."""
    mnfiles = 0
    occupations = []    # occupation per delta file
    deltanames = []     # names of delta files
    surfaceChecks = []  # is that atom at the surface or not?
    indices = []        # vib and geo indices from search, cast to 1D array
    # determine which atoms are "at the surface"
    surfats = sl.getSurfaceAtoms(rp)
    # make lists to print
    for at in rp.search_atlist:
        sps = [sp for sp in rp.searchpars if sp.atom == at]
        occpar = [sp for sp in sps if sp.mode == "occ"][0] # exactly one
        i = rp.searchpars.index(occpar)
        totalocc = 0.0
        for el in at.disp_occ:
            if at in surfats:
                surfaceChecks.append(1)
            else:
                surfaceChecks.append(0)
            o = at.disp_occ[el][config[i]-1]
            occupations.append(o)
            totalocc += o
            pl = [sp for sp in sps if sp.el == el]
            if len(pl) == 0:
                logger.error("No search parameters found for atom {}."
                              "Aborting...".format(at.oriN))
                rp.setHaltingLevel(2)
                return ""
            deltanames.append(pl[0].deltaname)
            mnfiles += 1
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
            indices.append((ngeo * vibind) + geoind + 1)
        vp = [sp for sp in sps if sp.el == "vac"]
        if len(vp) != 0 and totalocc < 1.0:
            occupations.append(1 - totalocc)
            deltanames.append(vp[0].deltaname)
            mnfiles += 1
            indices.append(1)
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
    param += "      PARAMETER (MNCONCS=1)\n"
    param += ("      PARAMETER (MNT0={}, MNCSTEP={}, MNATOMS=1)\n"
               .format(len(rp.ivbeams), rp.mncstep))
    try:
        with open(param_name, "w") as wf:
            wf.write(param)
    except:
        logger.error("Failed to write to PARAM file for superpos")
        raise

    # collect contrin output
    i3 = ff.FortranRecordWriter("I3")
    f74x10 = ff.FortranRecordWriter('10F7.4')
    contrin = (i3.write([mnfiles]) + "  0               no. of files, VarAll: "
              "all possible parameter combinations? (1/0)\n")
    contrin += "  1                  number of concentration steps\n"
    contrin += f74x10.write(occupations) + "\n"
    for i in range(0, len(deltanames)):
        contrin += deltanames[i].ljust(15) + "  1"
        contrin += (i3.write([surfaceChecks[i]]) + "  1  FILENAME(A15 !!!),"
                   "VARIATIONS,SURFACE?,FORMATTED?\n")
        contrin += i3.write([indices[i]]) + "\n"
    try:
        with open(contrin_name, "w") as wf:
            wf.write(contrin)
    except:
        logger.error("Failed to write " + contrin_name + " file. Execution "
                      "will continue...")

    return contrin