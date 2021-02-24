# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 13:30:12 2020

@author: Florian Kraushofer

Functions for reading and writing files relevant to the delta calculation
"""

import logging
import numpy as np
import fortranformat as ff
import os
import shutil

from viperleed.tleedmlib.files.beams import writeAUXBEAMS

logger = logging.getLogger("tleedm.files.iodeltas")


def checkDelta(filename, at, el, rp):
    """Checks whether a given delta file corresponds to the requested
    displacements of a given atom. Returns True or False."""
    eps = 1e-4
    fgeo = []                               # found geo disp
    fvib = []                               # found vib disp
    at.mergeDisp(el)
    if el == "vac":
        dgeo = [np.array([0.0, 0.0, 0.0])]  # requested geo disp
        dvib = [0.0]                        # requested vib disp
    else:
        if el in at.disp_geo:
            dgeo = at.disp_geo[el]
        else:
            dgeo = at.disp_geo["all"]
        if el in at.disp_vib:
            dvib = at.disp_vib[el]
        else:
            dvib = at.disp_vib["all"]
    try:
        with open(filename, "r") as rf:
            lines = rf.readlines()
    except FileNotFoundError:
        logger.error("Error reading file "+filename)
        raise
    try:
        nbeams = int(lines[1][0:3])  # number of beams
        nvar = int(lines[1][6:9])    # number of variations (geo*vib)
    except Exception:
        logger.error("Error parsing file " + filename)
        raise
    if nbeams != len(rp.ivbeams):
        return False
    if nvar != len(dgeo)*len(dvib):
        return False
    beams = []      # read beams from delta file
    atline = 2
    for i in range(atline, len(lines)):
        try:
            fl = [float(s) for s in lines[i].split()]
        except Exception:
            logger.error("Error parsing file "+filename)
            raise
        for j in range(0, int(len(fl)/2)):
            beams.append((fl[2*j], fl[(2*j) + 1]))
        if len(beams) == nbeams:
            atline = i+1
            break
    for (i, hk) in enumerate(beams):   # check beams
        if not rp.ivbeams[i].isEqual_hk(hk, eps=eps):
            return False
    atline += 1   # skip the line after beams
    # geo displacements start here
    rf74x10 = ff.FortranRecordReader("10F7.4")
    parselist = []
    repeats = False
    endgeo = False
    entrycount = 0
    for i in range(atline, len(lines)):
        try:
            fl = [f for f in rf74x10.read(lines[i]) if f is not None]
        except Exception:
            logger.error("Error parsing file "+filename)
            raise
        if len(fl) < 10:
            atline = i+1
            endgeo = True  # short line -> end of geo block
        parselist = parselist + fl
        while len(parselist) >= 3:
            v = parselist[:3]
            new = np.array([v[1], v[2], v[0]])
            parselist = parselist[3:]
            entrycount += 1
            if not repeats:
                append = True
                if fgeo and np.linalg.norm(new - fgeo[0]) < eps:
                    repeats = True
                    append = False
                if append:
                    fgeo.append(new)
                if len(fgeo) == nvar:
                    atline = i+1
                    endgeo = True
                    break
            elif not endgeo:
                if (entrycount-1) % len(fgeo) == 0:  # should repeat here
                    if np.linalg.norm(new - fgeo[0]) > eps:
                        atline = i
                        endgeo = True
                        break
        if endgeo:
            break
    ngeo = len(fgeo)
    # check geometry
    if ngeo != len(dgeo):
        return False
    for (i, gd) in enumerate(fgeo):
        if np.linalg.norm(gd - dgeo[i]) > eps:
            return False
    # vib displacement starts here
    nvib = nvar / ngeo
    if int(nvib) - nvib > 1e-4:
        logger.error("Error reading file "+filename+": number of geometry "
                     "variations found does not match header.")
        return False
    nvib = int(nvib)
    for i in range(atline, len(lines)):
        try:
            fl = [f for f in rf74x10.read(lines[i]) if f is not None]
        except Exception:
            logger.error("Error parsing file "+filename)
            raise
        parselist = parselist + fl
        while len(parselist) >= ngeo:
            fvib.append(parselist[0])
            if any([f != 0. for f in parselist[1:ngeo]]):
                logger.warning("File "+filename+": Found unexpected entries "
                               "in list of vibrational displacements.")
            parselist = parselist[ngeo:]
        if len(fvib) >= nvib:
            break
    # check vibrations:
    if el.lower() == "vac":
        voff = 0.
    else:
        voff = at.site.vibamp[el]
    for (i, f) in enumerate(fvib):
        bv = round(round(dvib[i] + voff, 4) / 0.529, 4)
        # in bohr radii; rounding twice to account for 1. writing to
        #  delta-input, 2. reading from DELTA file. Precision 0.529 taken from
        #  TensErLEED GLOBAL
        if abs(f - bv) >= 1e-4:
            return False
    return True


def generateDeltaInput(atom, targetel, sl, rp, deltaBasic="", auxbeams="",
                       phaseshifts=""):
    """
    Generates a PARAM file and delta input for one element of one atom.

    Parameters
    ----------
    atom : Atom
        Atom object for which input should be generated.
    targetel : str
        The element of that atom for which input should be generated. This may
        not be the main atom element.
    sl : Slab
        The Slab object containing atom information.
    rp : Rparams
        The run parameters object.
    deltaBasic : str, optional
        Part of delta input that is the same for all atoms. If not passed,
        will call generateDeltaBasic directly.
    auxbeams : str, optional
        The contents of the AUXBEAMS file. If not passed, will attempt to find
        the AUXBEAMS file and read it.
    phaseshifts : str, optional
        The contents of the _PHASESHIFTS file. If not passed, will attempt to
        find the _PHASESHIFTS file and read it.

    Returns
    -------
    (str, str, str).
        The delta input, a shortened version of that input for logging, and
        the contents of the required PARAM file.
    """

    if deltaBasic == "":
        deltaBasic = generateDeltaBasic()
    if auxbeams == "":
        # if AUXBEAMS is not in work folder, check AUX folder
        if not os.path.isfile(os.path.join(".", "AUXBEAMS")):
            if os.path.isfile(os.path.join(".", "AUX", "AUXBEAMS")):
                try:
                    shutil.copy2(os.path.join(".", "AUX", "AUXBEAMS"),
                                 "AUXBEAMS")
                except Exception:
                    logger.warning("Failed to copy AUXBEAMS from AUX folder")
            else:
                logger.warning("generateDeltaInput: AUXBEAMS not found")
        try:
            with open("AUXBEAMS", "r") as rf:
                auxbeams = rf.read()
            if auxbeams[-1] != "\n":
                auxbeams += "\n"
        except Exception:
            logger.error("generateDeltaInput: Could not read AUXBEAMS")
            raise
    if phaseshifts == "":
        try:
            with open("_PHASESHIFTS", "r") as rf:
                phaseshifts = rf.read()
            if phaseshifts[-1] != "\n":
                phaseshifts += "\n"
        except Exception:
            logger.error("generateDeltaInput: Could not read _PHASESHIFTS")
            raise
    MLMAX = [19, 126, 498, 1463, 3549, 7534, 14484, 25821, 43351, 69322,
             106470, 158067, 227969, 320664, 441320]
    try:
        beamlist, _, _ = writeAUXBEAMS(ivbeams=rp.ivbeams,
                                       beamlist=rp.beamlist, write=False)
    except Exception:
        logger.error("writeDeltaInput: Exception while getting data from "
                     "writeAUXBEAMS")
        raise

    # merge offsets with displacement lists
    atom.mergeDisp(targetel)

    # generate delta.in
    din = ("""   1                         FORMOUT - 1: formatted output
-------------------------------------------------------------------
--- chemical nature of displaced atom                           ---
-------------------------------------------------------------------
""")
    if targetel.lower() == "vac":
        iel = 0
    else:
        # find number of target element
        i = 0
        for el in sl.elements:
            # this reproduces the order of blocks contained in _PHASESHIFTS:
            if el in rp.ELEMENT_MIX:
                chemelList = rp.ELEMENT_MIX[el]
            else:
                chemelList = [el]
            siteList = [s for s in sl.sitelist if s.el == el]
            for cel in chemelList:
                for s in siteList:
                    i += 1
                    if s.isEquivalent(atom.site) and cel == targetel:
                        iel = i
    i4 = ff.FortranRecordWriter("I4")
    ol = i4.write([iel])
    din += ol.ljust(29) + "IEL  - element in PHASESHIFT list"
    din += """
-------------------------------------------------------------------
--- undisplaced position of atomic site in question             ---
-------------------------------------------------------------------
"""
    f74x3 = ff.FortranRecordWriter('3F7.4')
    ol = f74x3.write([0.0, 0.0, 0.0])
    din += ol.ljust(29) + "CUNDISP - displacement offset"
    din += """
-------------------------------------------------------------------
--- displaced positions of atomic site in question              ---
-------------------------------------------------------------------
"""
    if targetel == "vac":
        geolist = [(0., 0., 0.)]
    elif targetel in atom.disp_geo:
        geolist = atom.disp_geo[targetel]
    else:
        geolist = atom.disp_geo["all"]
    geosteps = len(geolist)
    ol = i4.write([geosteps])
    din += ol.ljust(29) + "NCSTEP - number of displaced positions\n"
    for disp in geolist:
        ol = f74x3.write([disp[2], disp[0], disp[1]])
        din += ol.ljust(29)+"CDISP(z,x,y) - z pointing towards bulk\n"
    din += (
        """-------------------------------------------------------------------
--- vibrational displacements of atomic site in question        ---
-------------------------------------------------------------------
""")
    if targetel == "vac":
        viblist = [0.]
    elif targetel in atom.disp_vib:
        viblist = atom.disp_vib[targetel]
    else:
        viblist = atom.disp_vib["all"]
    vibsteps = len(viblist)
    ol = i4.write([vibsteps])
    din += ol.ljust(29) + "NDEB - number of vib. amplitudes to be considered\n"
    f74 = ff.FortranRecordWriter("F7.4")
    if targetel.lower() != "vac":
        # "default" vibamp + offset, not just offset
        vibamps = [v + atom.site.vibamp[targetel] for v in viblist]
        if any([v <= 0 for v in vibamps]):
            logger.warning(
                "Vibrational amplitudes for {} contain values <= 0 "
                "(smallest: {:.4f}). Shifting displacement list to avoid "
                "non-positive numbers.".format(atom, min(vibamps)))
            corr = min([v for v in vibamps if v > 0]) - min(vibamps)
            vibamps = [v + corr for v in vibamps]
            for i in range(len(viblist)):
                # can't be done by list comprehension because it should modify
                #   the list that 'viblist' is pointing to, not make a copy
                viblist[i] += corr
    else:
        vibamps = [0.]
    for vibamp in vibamps:
        ol = f74.write([vibamp])
        din += ol.ljust(29)+"DRPER_A\n"

    din_main = deltaBasic + auxbeams + phaseshifts + din
    din_short = deltaBasic + "[AUXBEAMS]\n" + "[_PHASESHIFTS]\n" + din

    # write PARAM
    param = ("""C  Parameter statements for delta amplitude calculation, v1.2
C  parameters must be consistent with preceding reference calculation!

C  MLMAX: maximum angular momentum to be considered in calculation
C  MNLMB: number of Clebsh-Gordon coefficients needed in tmatrix() subroutine -
C         set according to current LMAX
"""
             "C         MLMAX:  1  2   3    4    5    6    7     8     9     "
             "10    11     12     13     14     15\n"
             "C         MNLMB: 19 126 498 1463 3549 7534 14484 25821 43351 "
             "69322 106470 158067 227969 320664 441320\n" """
C  MNPSI: number of phase shift values tabulated in phase shift file
C  MNEL : number of elements for which phase shifts are tabulated
C  MNT0 : number of beams for which delta amplitude calculation is required
C  MNATOMS: currently must be set to 1. In principle number of different atomic
C      positions in a superlattice wrt the reference periodicity when computing
C      TLEED beams for a superlattice not present in the reference structure
C  MNDEB: number of thermal variation steps to be performed (outer var. loop)
C  MNCSTEP: number of geometric variation steps to be performed """
             + "(inner var. loop)\n\n")
    param += "      PARAMETER( MLMAX = {} )\n".format(rp.LMAX)
    param += "      PARAMETER( MNLMB = {} )\n".format(MLMAX[rp.LMAX-1])
    param += ("      PARAMETER( MNPSI = {}, MNEL = {} )\n"
              .format(len(rp.phaseshifts), (len(rp.phaseshifts[0][1]))))
    param += "      PARAMETER( MNT0 = {} )\n".format(len(beamlist))
    param += "      PARAMETER( MNATOMS = 1 )\n"
    param += "      PARAMETER( MNDEB = {} )\n".format(vibsteps)
    param += "      PARAMETER( MNCSTEP = {} )\n".format(geosteps)
    return din_main, din_short, param


def generateDeltaBasic(sl, rp):
    """Generates the part of the input for delta-amplitudes that is the same
    for all atoms, and returns it as a string."""
    output = ""
    output += rp.systemName+" "+rp.timestamp+"\n"
    f72x2 = ff.FortranRecordWriter('2F7.2')
    ol = f72x2.write([rp.THEO_ENERGIES[0], rp.THEO_ENERGIES[1]+0.01])
    output += ol.ljust(24) + "EI,EF\n"
    f74x2 = ff.FortranRecordWriter('2F7.4')
    ucsurf = sl.ucell[:2, :2].T
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    ucbulk = sl.bulkslab.ucell[:2, :2].T
    ol = f74x2.write(ucbulk[0])
    output += ol.ljust(24) + 'ARA1\n'
    ol = f74x2.write(ucbulk[1])
    output += ol.ljust(24) + 'ARA2\n'
    ol = f74x2.write(ucsurf[0])
    output += ol.ljust(24) + 'ARB1\n'
    ol = f74x2.write(ucsurf[1])
    output += ol.ljust(24) + 'ARB2\n'
    ol = f72x2.write([rp.THETA, rp.PHI])
    output += ol.ljust(24) + "THETA PHI\n"
    return output
