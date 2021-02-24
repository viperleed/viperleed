# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 17:20:43 2020

@author: Florian Kraushofer

Functions for reading and writing the PHASESHIFTS file
"""

import logging
import numpy as np

from viperleed import fortranformat as ff

logger = logging.getLogger("tleedm.files.phaseshifts")

def readPHASESHIFTS(sl, rp, readfile='_PHASESHIFTS', check=True, 
                    ignoreEnRange=False):
    """Reads from a _PHASESHIFTS file, returns the data as a list of tuples
    (E, enps), where enps is a list of lists, containing one list of values
    (for different L) for each element. Therefore, len(phaseshifts) is the
    number of energies found, len(phaseshifts[0][1]) should match the number
    of elements, and len(phaseshifts[0][1][0]) is the number of different
    values of L in the phaseshift file. The "check" parameters controls
    whether the phaseshifts that were found should be checked against the
    parameters / slab. If it is set to False, then passing "None" as sl and rp
    will work."""
    rf74x10 = ff.FortranRecordReader('10F7.4')
    ri3 = ff.FortranRecordReader('I3')

    try:
        rf = open(readfile, 'r')
    except FileNotFoundError:
        logger.error("_PHASESHIFTS file not found.")
        raise

    filelines = []
    for line in rf:
        filelines.append(line[:-1])
    rf.close()

    try:
        nel = ri3.read(filelines[0])[0]
    except:
        logger.error("Exception while trying to read _PHASESHIFTS: could not "
                      "find number of blocks in first line.")
        raise
    phaseshifts = []

    firstline = filelines[0]
    readline = 1
    linesperblock = 0
    while readline < len(filelines):
        if linesperblock:
            en = rf74x10.read(filelines[readline])[0]
            enps = []
            for i in range(0,nel):
                elps = []
                for j in range(0,linesperblock):
                    llist = rf74x10.read(filelines[readline
                                                   +(i*linesperblock)+j+1])
                    llist = [f for f in llist if f is not None]
                    elps.extend(llist)
                enps.append(elps)
            phaseshifts.append((en,enps))
            readline += linesperblock*nel+1
        else:
            #first check how many lines until the next energy:
            lineit = 1
            llist = rf74x10.read(filelines[readline+lineit])
            llist = [f for f in llist if f is not None]
            longestline = len(llist)
            shortestline = longestline
            lastlen = longestline
            cont = True
            while cont:
                lineit += 1
                llist = rf74x10.read(filelines[readline+lineit])
                llist = [f for f in llist if f is not None]
                if len(llist) == 1:
                    if lastlen == 1 or (shortestline > 1
                                        and shortestline < longestline):
                        cont = False  #found next energy
                    else:
                        shortestline = 1
                        #blocklines += 1
                elif len(llist) != longestline:
                    shortestline = len(llist)
                lastlen = len(llist)
            linesperblock = int((lineit-1)/nel)
            if not linesperblock or (((lineit-1)/nel) - linesperblock != 0.0):
                logger.warning("Error while trying to read _PHASESHIFTS: "
                    "Could not parse file: The number of blocks may not match "
                    "the number given in the first line. A new _PHASESHIFTS "
                    "file will be generated.")
                rp.setHaltingLevel(1)
                return ("", [], True, True)
            # don't increase readline -> read the same block again afterwards

    if not check:
        newpsGen, newpsWrite = False, False
    else:
        # check whether the phaseshifts that were found fit the data:
        newpsGen, newpsWrite = True, True
                    # recommend that new values should be generated / written
        psblocks = 0
        for el in sl.elements:
            if el in rp.ELEMENT_MIX:
                n = len(rp.ELEMENT_MIX[el])
            else:
                n = 1
            psblocks += n*len([s for s in sl.sitelist if s.el == el])
        # check for MUFTIN parameters:
        muftin = True
        llist = firstline.split()
        if len(llist) >= 6:
            for i in range(1,5):
                try:
                    float(llist[i])
                except ValueError:
                    muftin = False
        else:
            muftin = False
        if rp.V0_REAL == "default" and not muftin:
            logger.warning("Could not convert first line of "
                "_PHASESHIFTS file to MUFTIN parameters. A new "
                "_PHASESHIFTS file will be generated.")
            rp.setHaltingLevel(1)
        elif len(phaseshifts[0][1]) == psblocks:
            logger.debug("Found "+str(psblocks)+" blocks in _PHASESHIFTS "
                          "file, which is consistent with PARAMETERS.")
            newpsGen, newpsWrite = False, False
        elif len(phaseshifts[0][1]) == sl.nelem:
            logger.warning("Found fewer blocks than expected in the "
                "_PHASESHIFTS file. However, the number of blocks matches "
                "the number of chemical elements. A new _PHASESHIFTS file "
                "will be generated, assuming that each block in the old "
                "file should be used for all atoms of one element.")
            rp.setHaltingLevel(1)
            oldps = phaseshifts[:]
            phaseshifts = []
            for (en,oldenps) in oldps:
                enps = []
                j = 0   #block index in old enps
                for el in sl.elements:
                    if el in rp.ELEMENT_MIX:
                        m = len(rp.ELEMENT_MIX[el])
                    else:
                        m = 1
                    n = len([s for s in sl.sitelist if s.el == el])
                    for i in range(0,m):    # repeat for chemical elements
                        for k in range(0, n):   # repeat for sites
                            enps.append(oldenps[j])
                        j += 1  # count up the block in old enps
                phaseshifts.append((en,enps))
            newpsGen = False
            firstline = str(len(phaseshifts[0][1])).rjust(3) + firstline[3:]
        else:
            logger.warning("_PHASESHIFTS file was read but is "
                "inconsistent with PARAMETERS. A new _PHASESHIFTS file "
                "will be generated.")
            rp.setHaltingLevel(1)

    if check and not ignoreEnRange:
        # check whether energy range is large enough:
        checkfail = False
        er = np.arange(rp.THEO_ENERGIES[0], rp.THEO_ENERGIES[1]+1e-4, 
                       rp.THEO_ENERGIES[2])
        psmin = round(phaseshifts[0][0]*27.2116, 2)
        psmax = round(phaseshifts[-1][0]*27.2116, 2)
        if rp.V0_REAL == "default":
            llist = firstline.split()
            c = []
            try:
                for i in range(0,4):
                    c.append(float(llist[i+1]))
            except:
                checkfail = True
            else:
                er_inner = [e + (rp.FILAMENT_WF - max(c[0], 
                                        c[1] + (c[2]/np.sqrt(e + c[3] 
                                                          + rp.FILAMENT_WF))))
                            for e in er] # energies at which scattering occurs
        else:
            try:
                v0r = float(rp.V0_REAL)
            except:
                checkfail = True
            else:
                er_inner = [e + v0r for e in er]
        if not checkfail:
            if (psmin > min(er_inner) or psmax < max(er_inner)):
                if (psmin > min(er_inner) and psmin <= 20. 
                                          and psmax >= max(er_inner)):
                    # can lead to re-calculation of phaseshifts every run if 
                    #  V0r as calculated by EEASiSSS differs from 'real' V0r.
                    #  Don't automatically correct.
                    logger.warning("Lowest value in _PHASESHIFTS file ({:.1f} "
                        "eV) is larger than the lowest predicted scattering "
                        "energy ({:.1f} eV). If this causes problems in the "
                        "reference calculation, try deleting the _PHASESHIFTS "
                        "file to generate a new one, or increase the starting "
                        "energy in the THEO_ENERGIES parameter."
                        .format(psmin, min(er_inner)))
                else:
                    logger.warning("The energy range found in the _PHASESHIFTS"
                        " file is smaller than the energy range requested for "
                        "theoretical beams. A new _PHASESHIFTS file will be "
                        "generated.")
                    newpsGen, newpsWrite = True, True
        else:
            logger.warning("Could not check energy range in _PHASESHIFTS "
                "file. If energy range is insufficient, try deleting the "
                "_PHASESHIFTS file to generate a new one.")
    return (firstline, phaseshifts, newpsGen, newpsWrite)

def writePHASESHIFTS(firstline, phaseshifts, filename='_PHASESHIFTS'):
    """Takes phaseshift data and writes it to a _PHASESHIFTS file."""
    output = firstline
    if output[-1] != "\n":
        output += "\n"
    f74x10 = ff.FortranRecordWriter('10F7.4')
    f74 = ff.FortranRecordWriter('F7.4')
    for (en,enps) in phaseshifts:
        output += f74.write([en])+"\n"
        for block in enps:
            output += f74x10.write(block)+"\n"
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
        logger.debug("Wrote to "+filename+" successfully.")
    except:
        logger.error("Exception while writing _PHASESHIFTS file: ", 
                      exc_info=True)
        raise
    return