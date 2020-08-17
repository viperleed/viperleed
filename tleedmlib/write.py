# -*- coding: utf-8 -*-
"""
Created on Jun 13 2019

@author: Florian Kraushofer

Contains functions for writing the output files of the TensErLEED interface
"""

import logging
import numpy as np
import fortranformat as ff
import os
import copy
import random
import shutil

import tleedmlib as tl
from guilib.base import project_to_first_domain

logger = logging.getLogger("tleedm.write")

try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')  # !!! check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
except:
    plotting = False
else:
    plotting = True

def modifyPARAMETERS(rp, modpar, new="", comment=""):
    """Looks for 'modpar' in the PARAMETERS file, comments that line out, and 
    replaces it by the string specified by 'new'"""
    oriname = "PARAMETERS_ori_"+rp.timestamp
    if not oriname in rp.manifest:
        try:
            shutil.copy2("PARAMETERS", oriname)
        except:
            logger.error("modifyPARAMETERS: Could not copy PARAMETERS file "
                "to PARAMETERS_ori. Proceeding, original file will be lost.")
        rp.manifest.append(oriname)
    if not "PARAMETERS" in rp.manifest:
        rp.manifest.append("PARAMETERS")
    output = ""
    headerPrinted = False

    try:
        with open("PARAMETERS", "r") as rf:
            plines = rf.readlines()
    except:
        logger.error("Error reading PARAMETERS file.")
        raise
    found = False
    for line in plines:
        if "! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #" in line:
            headerPrinted = True
        valid = False
        param = ""
        if "=" in line:    #ignore all lines that don't have an "=" sign at all
            param = line.split('=')[0]        #parameter is defined left of "="
            if param: 
                valid = True
                plist = param.split()  
                if plist: 
                    param = plist[0]
                if param[0] == '!': 
                    valid = False
        if valid and param == modpar:
            found = True
            if new:
                if comment == "":
                    comment = "line automatically changed to:"
                output += "!"+line[:-1] + " ! " + comment + "\n"
                output += new + "\n"
            else:
                if comment == "":
                    comment = "line commented out automatically"
                output += ("!"+line.rstrip()).ljust(35)+ " ! " + comment + "\n"
        else:
            output += line
    if not found:
        if not headerPrinted:
            output += """

! ######################################################
! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #
! ######################################################
"""
            headerPrinted = True
        output += "\n" + modpar + " = " + new + " ! " + comment
    try:
        with open("PARAMETERS", "w") as wf:
            wf.write(output)
    except:
        logger.error("modifyPARAMETERS: Failed to write PARAMETERS file.")
        raise
    return



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
    return 0

def writeSuperposPARAM(rp, filename = "PARAM"):
    """Writes a PARAM file for the Superpos calculation."""
    mnfiles = 0
    for at in rp.search_atlist:
        mnfiles += len(at.deltasGenerated)
        
    output = """C  DIMENSIONS MAY BE CHANGED IN PARAMETER-STATEMENT
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
    output += "      PARAMETER (MNFILES={})\n".format(mnfiles)
    output += "      PARAMETER (MNCONCS=1)\n"
    output += ("      PARAMETER (MNT0={}, MNCSTEP={}, MNATOMS=1)\n"
               .format(len(rp.ivbeams), rp.mncstep))
    try:
        with open(filename, "w") as wf:
            wf.write(output)
    except:
        logger.error("Failed to write to PARAM file for superpos")
        raise
    return 0

def writeCONTRIN(sl, rp, config, filename = "superpos-CONTRIN", write=True):
    """Writes a CONTRIN file as input for the Superpos calculation."""
    mnfiles = 0
    occupations = []    # occupation per delta file
    deltanames = []     # names of delta files
    surfaceChecks = []  # is that atom at the surface or not?
    indices = []        # vib and geo indices from search, cast to 1D array
    # determine which atoms are "at the surface"
    surfats = sl.getSurfaceAtoms(rp)
    # make lists to print
    for at in rp.search_atlist:
        mnfiles += len(at.deltasGenerated)
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
            ind = (ngeo * vibind) + geoind + 1
            indices.append(ind)
        vp = [sp for sp in sps if sp.el == "vac"]
        if len(vp) != 0 and totalocc < 1.0:
            occupations.append(1 - totalocc)
            deltanames.append(vp[0].deltaname)
            indices.append(1)
            if at in surfats:
                surfaceChecks.append(1)
            else:
                surfaceChecks.append(0)
            
    # collect output
    i3 = ff.FortranRecordWriter("I3")
    f74x10 = ff.FortranRecordWriter('10F7.4')
    output = (i3.write([mnfiles]) + "  0               no. of files, VarAll: "
              "all possible parameter combinations? (1/0)\n")
    output += "  1                  number of concentration steps\n"
    output += f74x10.write(occupations) + "\n"
    for i in range(0, len(deltanames)):
        output += deltanames[i].ljust(15) + "  1"
        output += (i3.write([surfaceChecks[i]]) + "  1  FILENAME(A15 !!!),"
                   "VARIATIONS,SURFACE?,FORMATTED?\n")
        output += i3.write([indices[i]]) + "\n"
    try:
        with open(filename, "w") as wf:
            wf.write(output)
    except:
        logger.error("Failed to write " + filename + " file. Execution will "
                      "continue...")
    return output
    

def writeRfInfo(sl, rp, filename="rf.info"):
    """Generates r-factor parameters for the search, combines them with the
    experimental beams in AUXEXPBEAMS format to make the entire input for the 
    search, returns that as a string. If 'write' is True, will also write the 
    result to 'filename'."""
    expEnergies = []
    for b in rp.expbeams:
        expEnergies.extend([k for k in b.intens if k not in expEnergies])
    expEnergies.sort()
    minen = min(min(expEnergies), rp.THEO_ENERGIES[0])
    maxen = max(max(expEnergies), rp.THEO_ENERGIES[1])
    # extend energy range if they are close together
    if abs(min(expEnergies) - rp.THEO_ENERGIES[0]) < abs(rp.IV_SHIFT_RANGE[0]):
        minen = (max(min(expEnergies), rp.THEO_ENERGIES[0]) 
                 - rp.IV_SHIFT_RANGE[0])
    if abs(max(expEnergies) - rp.THEO_ENERGIES[1]) < abs(rp.IV_SHIFT_RANGE[1]):
        maxen = (min(max(expEnergies), rp.THEO_ENERGIES[1]) 
                 + rp.IV_SHIFT_RANGE[1]) + 0.01
    step = min(expEnergies[1]-expEnergies[0], rp.THEO_ENERGIES[2])
    if rp.IV_SHIFT_RANGE[2] > 0:
        vincr = rp.IV_SHIFT_RANGE[2]
        step = min(step, vincr)
    else:
        vincr = step
    # find correspondence experimental to theoretical beams:
    beamcorr = tl.leedbase.getBeamCorrespondence(sl, rp)
    # integer & fractional beams
    iorf = []
    for beam in rp.expbeams:
        if beam.hkfrac[0].denominator != 1 or beam.hkfrac[1].denominator != 1:
            iorf.append(2)
        else:
            iorf.append(1)

    f72 = ff.FortranRecordWriter("F7.2")
    i3 = ff.FortranRecordWriter("I3")
    i3x25 = ff.FortranRecordWriter("25I3")
    f31x25 = ff.FortranRecordWriter("25F3.1")
    output = (f72.write([minen]).ljust(16) + "EMIN\n")
    output += (f72.write([maxen]).ljust(16) + "EMAX\n")
            # !!! BULLSHIT RESULTS WHEN EMAX > MAX ENERGY IN DELTA FILES
    output += (f72.write([step]).ljust(16) + "EINCR\n")
    output += "  0             IPR - determines amount of output to stdout  \n"
    output += (f72.write([rp.V0_IMAG]).ljust(16) + "VI\n")
            # !!! bulk or surface V0_IMAG?
    output += "   0.00         V0RR\n"   
    output += (f72.write([rp.IV_SHIFT_RANGE[0]]).ljust(16) + "V01\n")
    output += (f72.write([rp.IV_SHIFT_RANGE[1]]).ljust(16) + "V02\n")
    output += (f72.write([vincr]).ljust(16) + "VINCR\n")
    output += (i3.write([rp.R_FACTOR_SMOOTH]) + "             ISMOTH\n")
    output += "  0             EOT - 0: exp format, 1: van Hove format\n"
    output += ((i3.write([len(rp.ivbeams)]) + i3.write([len(rp.expbeams)]))
               .ljust(16) + "NTH NEX\n")
    # numbers of theoretical beams, as they correspond to experimental beams
    output += (i3x25.write([n+1 for n in beamcorr]) + "\n")
    output += " DATA MITTEL (integer & half order beams) :\n"
    output += i3x25.write(iorf) + "\n"
    output += " exp - th relationship IBP, beam weights WB\n"
    output += i3x25.write([n+1 for n in range(0,len(rp.expbeams))]) + "\n"
    output += f31x25.write([1]*len(rp.expbeams))+"\n"
    auxexpbeams = writeAUXEXPBEAMS(rp.expbeams, header = rp.systemName,
                                   write = True, numbers = True)
    output += auxexpbeams
    
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to "+filename+" successfully")
    return output

def writeWEXPEL(sl, rp, theobeams, filename="WEXPEL"):
    """Writes input file WEXPEL for R-factor calculation."""
    theoEnergies = []
    for b in theobeams:
        theoEnergies.extend([k for k in b.intens if k not in theoEnergies])
    theoEnergies.sort()
    expEnergies = []
    for b in rp.expbeams:
        expEnergies.extend([k for k in b.intens if k not in expEnergies])
    expEnergies.sort()
    minen = max(min(expEnergies), rp.THEO_ENERGIES[0])
    maxen = min(max(expEnergies), rp.THEO_ENERGIES[1])
    # extend energy range if they are close together
    if abs(min(expEnergies) - rp.THEO_ENERGIES[0]) < abs(rp.IV_SHIFT_RANGE[0]):
        minen = (max(min(expEnergies), rp.THEO_ENERGIES[0]) 
                 - rp.IV_SHIFT_RANGE[0])
    if abs(max(expEnergies) - rp.THEO_ENERGIES[1]) < abs(rp.IV_SHIFT_RANGE[1]):
        maxen = (min(max(expEnergies), rp.THEO_ENERGIES[1]) 
                 + rp.IV_SHIFT_RANGE[1]) + 0.01
    step = min(expEnergies[1]-expEnergies[0], theoEnergies[1]-theoEnergies[0])
    if rp.IV_SHIFT_RANGE[2] > 0:
        vincr = rp.IV_SHIFT_RANGE[2]
        step = min(step, vincr)
    else:
        vincr = step
    # find correspondence experimental to theoretical beams:
    beamcorr = tl.leedbase.getBeamCorrespondence(sl, rp)
    # integer & fractional beams
    iorf = []
    for (i, beam) in enumerate(rp.ivbeams):
        if beamcorr[i] == -1:
            iorf.append(0)
        elif beam.hk[0] % 1.0 != 0.0 or beam.hk[1] % 1.0 != 0.0:
            iorf.append(2)
        else:
            iorf.append(1)

    f72 = ff.FortranRecordWriter("F7.2")
    i3x25 = ff.FortranRecordWriter("25I3")
    i3 = ff.FortranRecordWriter("I3")
    output = " &NL1\n"
    output += (" EMIN=" + f72.write([minen]).rjust(9) + ",\n")
    output += (" EMAX=" + f72.write([maxen]).rjust(9) + ",\n")
    output += (" EINCR=" + f72.write([step]).rjust(8) + ",\n")
    output += " LIMFIL=      1,\n" # number of consecutive input files
    output += " IPR=         0,\n" # output formatting
    output += (" VI=" + f72.write([rp.V0_IMAG]).rjust(11) + ",\n")
    output += " V0RR=      0.0,\n"
    output += (" V01=" + f72.write([rp.IV_SHIFT_RANGE[0]]).rjust(10) + ",\n") 
    output += (" V02=" + f72.write([rp.IV_SHIFT_RANGE[1]]).rjust(10) + ",\n")
    output += (" VINCR=" + f72.write([vincr]).rjust(8) + ",\n")
    output += " ISMOTH=" + i3.write([rp.R_FACTOR_SMOOTH]).rjust(7) + ",\n"
    output += " EOT=         0,\n"
    output += " PLOT=        1,\n"
    output += " GAP=         0,\n"
    output += " &END\n"
    output += (i3x25.write([n+1 for n in beamcorr]) + "\n")
    for i in range(0,2):   # redundant since indices are already taken care of
        output += i3x25.write([n+1 for n in range(0,len(rp.expbeams))]) + "\n"
    output += i3x25.write(iorf) + "\n"
    output += "\n&NL2\n"
    output += " NSSK=    0,\n"
    if rp.R_FACTOR_TYPE == 1:
        output += " WR=      0.,0.,1.,\n" # Pendry
    elif rp.R_FACTOR_TYPE == 2:
        output += " WR=      1.,0.,0.,\n" # R2
    else:
        output += " WR=      0.,1.,0.,\n" # Zanazzi-Jona
    output += """ &END
 &NL3
 NORM=           1,
 INTMAX=    999.99,
 PLSIZE=   1.0,1.0,
 XTICS=         50,
 &END
 """
    auxexpbeams = writeAUXEXPBEAMS(rp.expbeams, header = rp.systemName,
                                   write = True, numbers = False)
    output += auxexpbeams
    output += "\n"  # information about gaps in the experimental spectra would
                    #   go here
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to R-factor input file "+filename+" successfully")

def writeRfactPARAM(rp, theobeams):
    """Generates the PARAM file for the rfactor calculation."""
    theoEnergies = []
    for b in theobeams:
        theoEnergies.extend([k for k in b.intens if k not in theoEnergies])
    theoEnergies.sort()
    expEnergies = []
    for b in rp.expbeams:
        expEnergies.extend([k for k in b.intens if k not in expEnergies])
    expEnergies.sort()
    minen = min(min(expEnergies), min(theoEnergies))
    maxen = max(max(expEnergies), max(theoEnergies))
    if rp.IV_SHIFT_RANGE[2] > 0:
        step = rp.IV_SHIFT_RANGE[2]
    else:
        step = min(expEnergies[1]-expEnergies[0], rp.THEO_ENERGIES[2])
    ngrid = int(np.ceil(((maxen-minen)/step)*1.1))
    output = """
C  MNBED  : number of beams in experimental spectra before averaging
C  MNBTD  : number of beams in theoretical spectra before averaging

      PARAMETER (MNBED = {}, MNBTD = {})""".format(len(rp.expbeams), 
                                                   len(theobeams))
    output += """

C  MNET   : number of energies in theoretical beam at time of reading in
C  MNGP   : greater equal number of grid points in energy working grid (ie after
C           interpolation)
C  MNS    : number of geometries including those, that are skipped

      PARAMETER (MNET = {}, MNGP = {})""".format(len(theoEnergies), ngrid)
    output += """
      PARAMETER (MNS = 1)

C  MNGAP  : number of gaps in the experimental spectra (NOTE: if there are no
C           gaps in the spectra set MNGAP to 1 to avoid zero-sized arrays)

      PARAMETER (MNGAP = 1)
"""         # !!! check - any case in which MSN or MNGAP are != 1?
    # write PARAM
    try:
        with open("PARAM", "w") as wf:
            wf.write(output)
    except:
        logger.error("Failed at writing PARAM file for R-factor calculation.")
        raise

def generateSearchInput(sl, rp, steuOnly=False, cull=False):
    """Generates a PARAM and a search.steu file for the search. If steuOnly is 
    set True, PARAM and restrict.f will not be written. If cull is set True, 
    part of the population (given by rp.SEARCH_CULL) will be removed and 
    replaced by copies of random surviving members of the population."""
    # first generate list of SearchPar objects and figure out search parameters
    r = rp.generateSearchPars(sl, rp)
    if r != 0:
        logger.debug(r)
        logger.error("Error getting search parameters.")
        return 1
    
    # if search population is undefined, calculate a default:
    if rp.SEARCH_POPULATION == 0:
        spop = min(48, 15 + rp.indyPars)
        rp.SEARCH_POPULATION = int(np.ceil(spop / rp.N_CORES)) * rp.N_CORES
    
    # calculate some more things for later
    expEnergies = []
    for b in rp.expbeams:
        expEnergies.extend([k for k in b.intens if k not in expEnergies])
    theoEnergies = int((rp.THEO_ENERGIES[1]-rp.THEO_ENERGIES[0])
                       / rp.THEO_ENERGIES[2]) + 1
    if theoEnergies >= len(expEnergies):
        logger.warning("Theoretical beams have more data points than "
                        "experimental beams")
        rp.setHaltingLevel(1)
        # TODO: will this crash? If so, return here
    
    # determine which atoms are "at the surface"
    surfats = sl.getSurfaceAtoms(rp)

    # merge offsets with displacement lists
    for at in rp.search_atlist:
        if at.oriState is None:
            for el in at.offset_occ:
                if el not in at.disp_occ:
                    logger.error("Atom "+str(at.oriN)+" has occupation "
                            "offset defined for element "+el+", which does "
                            "not appear in atom displacement list. Value "
                            "will be skipped, this may cause errors.")
                else:
                    at.disp_occ[el] = [v + at.offset_occ[el] 
                                       for v in at.disp_occ[el]]
                del at.disp_occ[el]
            for (d, o) in [(at.disp_geo, at.offset_geo),
                           (at.disp_vib, at.offset_vib),
                           (at.disp_geo, at.disp_geo_offset)]:
                dl = []
                for el in o:
                    if el not in d:
                        d[el] = copy.copy(d["all"])
                    d[el] = [v + o[el] for v in d[el]]
                    dl.append(el)
                for el in dl:
                    if o != at.disp_geo_offset:
                        del o[el]
            at.disp_geo_offset = {"all": [np.array([0.0,0.0,0.0])]}

    # PARAM
    output = """C Here are some global parameterization to be performed

C MNBED IS NUMBER OF EXPERIMENTAL BEAMS
C MNBTD IS NUMBER OF THEORETICAL BEAMS
C MNBMD IS MAX(MNBED,MNBTD)
      PARAMETER(MNBED = {})\n""".format(len(rp.expbeams))
    output += "      PARAMETER(MNBTD = {})\n".format(len(rp.ivbeams))
    output += "      PARAMETER(MNBMD = {})".format(max(len(rp.expbeams),
                                                         len(rp.ivbeams)))
    output += """
C MNDATA IS MAX. NUMBER OF DATA POINTS IN EXPERIMENTAL BEAMS
      PARAMETER(MNDATA = {})""".format(int(len(expEnergies)*1.1))
              # array size parameter only - add some buffer -> *1.1
    output += """
C MNDATT IS NUMBER OF THEORETICAL DATA POINTS IN EACH BEAM
      PARAMETER(MNDATT = {})""".format(int(theoEnergies*1.1))
    output += """
C MPS IS POPULATION SIZE (number of independent trial structures)
      PARAMETER(MPS = {})""".format(rp.SEARCH_POPULATION)
    output += """
C MNDOM is number of domains to be incoherently averaged
      parameter (MNDOM = 1)"""  # !!! HARDCODED FOR NOW - ADD LATER !!!
    output += """
C MNPLACES IS NUMBER OF DIFFERENT ATOMIC SITES IN CURRENT VARIATION
      PARAMETER(MNPLACES = {})""".format(len(rp.search_atlist))
    output += """
C MNFILES IS MAXIMUM NUMBER OF TLEED-FILES PER SITE
      PARAMETER(MNFILES = {})""".format(rp.search_maxfiles)
    output += """
C MNCONCS IS MAX NUMBER OF CONCENTRATION STEPS PER SITE
      PARAMETER(MNCONCS = {})\n""".format(rp.search_maxconc)
    output += ("C MNPRMK IS NUMBER OF PARAMETERS - INCL ONE CONC PARAMETER "
               "PER SITE - incl. 1 per domain!\n")
    output += "      PARAMETER(MNPRMK = {})".format(len(rp.searchpars)+1)
                                # !!! the "+1" comes from hardcoding domains!
    output += """
C MNCSTEP IS MAX NUMBER OF VARIATIONS (geo. times therm.) IN 1 FILE
C MNATOMS IS RELICT FROM OLDER VERSIONS 
      PARAMETER (MNCSTEP = {}, MNATOMS = 1)""".format(rp.mncstep)
    # write PARAM
    if not steuOnly:
        try:
            with open("PARAM", "w") as wf:
                wf.write(output)
        except:
            logger.error("Failed at writing PARAM file for search.")
            raise
    
    # now search.steu
    output = ""
    info = ("#".rjust(4) + "label".rjust(7) + "atom".rjust(7) + "elem".rjust(7)
            + "mode".rjust(7) + "steps".rjust(7) + "constr".rjust(7) + "\n")
                # while producing search.steu, collect information on 
                #  parameters for searchpars.info
    nsteps = []     # keep track of number of steps per parameter for init
    parcount = 0
    i3 = ff.FortranRecordWriter("I3")
    i1 = ff.FortranRecordWriter("I1")
    output += (i3.write([rp.indyPars]).ljust(16)
               + "number of independent parameters\n")
    f74 = ff.FortranRecordWriter('F7.4')
    output += (f74.write([rp.GAUSSIAN_WIDTH]).ljust(16)
               + "gaussian width control parameter RMUT\n")
    output += ("  0             initialisation for random number generator - "
               "0: use system time, 1,...: use init\n")
    output += (i3.write([rp.R_FACTOR_TYPE]).ljust(16)
               + "1: use RPe -- 2: use R2\n")
    output += (i3.write([rp.SEARCH_BEAMS]).ljust(16)
               + "Optimization of which beam group do you want? "
               "(0=Aver,1=Int,2=Half)\n")
    output += " 1000           output intervall\n"
    output += (str(rp.SEARCH_MAX_GEN).ljust(16) + "desired number of "
               "generations to be performed\n")
    output += """100             area fraction step width (%)
SD.TL           name of search document file (max. 10 characters)
  1             Number of domains under consideration
"""             # !!! some more hardcoded domain parameters here
    output += ("======= Information about Domain 1: =========================="
               "==================\n")
    output += (i3.write([len(rp.search_atlist)]).ljust(16)
               + "Number of atomic sites in variation: Domain 1\n")
    displistcount = len(sl.displists)+1
    for (i, at) in enumerate(rp.search_atlist):
        printvac = False
        output += ("------- Information about site {}: -----------------------"
                   "-----------------------\n".format(i+1))
        surf = 1 if at in surfats else 0
        output += (i3.write([surf]).ljust(16) + "Surface (0/1)\n")
        if at.displist in sl.displists:
            dlind = sl.displists.index(at.displist) + 1
        else:
            dlind = displistcount
            displistcount += 1
        output += (i3.write([dlind]).ljust(16) + "Atom number\n")
        output += (i3.write([len(at.deltasGenerated)]).ljust(16) 
                   + "No. of different files for Atom no. {}\n".format(i+1))
        for (j, deltafile) in enumerate(at.deltasGenerated):
            output += "****Information about file {}:\n".format(j+1)
            output += (deltafile.ljust(16) + "Name of file {} (max. 15 "
                       "characters)\n".format(j+1))
            output += "  1             Formatted(0/1)\n"
            el = deltafile.split("_")[2]
            if el.lower() == "vac":
                geo = 1
                vib = 0
                types = 1
                printvac = True
            else:
                geo0 = True
                vib0 = True
                for (mode, d) in [(1, at.disp_geo), (2, at.disp_vib)]:
                    if el in d:
                        dl = d[el]
                    else:
                        dl = d["all"]
                    if mode == 1:
                        geo = len(dl)
                        if geo == 1 and np.linalg.norm(dl[0]) != 0.:
                            geo0 = False
                    else:
                        vib = len(dl)
                        if vib == 1 and dl[0] != 0.:
                            vib0 = False
                if vib == 1 and vib0:
                    vib = 0
                elif geo == 1 and geo0:
                    geo = 0
                if geo > 0 and vib > 0:
                    types = 2
                else:
                    types = 1
            output += (i3.write([types]).ljust(16) + "Types of parameters in "
                       "file {}\n".format(j+1))
            label = i1.write([i+1])+i1.write([j+1])
            constr = {"vib": "-", "geo": "-"}
            for mode in ["vib", "geo"]:
                spl = [s for s in rp.searchpars if s.el == el 
                      and s.atom == at and s.mode == mode]
                if spl and spl[0].restrictTo is not None:
                    sp = spl[0]
                    if type(sp.restrictTo) == int:
                        constr[mode] = str(sp.restrictTo)
                    else:
                        constr[mode] = "#"+str(rp.searchpars.index(
                                                            sp.restrictTo)+1)
                elif spl and spl[0].linkedTo is not None:
                    constr[mode] = "#"+str(rp.searchpars.index(
                                                            spl[0].linkedTo)+1)
            if vib > 0:
                output += (i3.write([vib]).ljust(16) + "vibrational steps\n")
                parcount += 1

                info += (str(parcount).rjust(4) + ("P"+label).rjust(7) 
                         + str(at.oriN).rjust(7) + el.rjust(7) + "vib".rjust(7)
                         + str(vib).rjust(7) + constr["vib"].rjust(7) + "\n")
                nsteps.append(vib)
            if geo > 0:
                output += (i3.write([geo]).ljust(16) + "geometrical steps\n")
                parcount += 1
                info += (str(parcount).rjust(4) + ("P"+label).rjust(7) 
                         + str(at.oriN).rjust(7) + el.rjust(7) + "geo".rjust(7)
                         + str(geo).rjust(7) + constr["geo"].rjust(7) + "\n")
                nsteps.append(geo)
        output += "****concentration steps for site no. {}\n".format(i+1)
        occsteps = len(next(iter(at.disp_occ.values())))
        output += (i3.write([occsteps]).ljust(16) + "no. of concentration "
                   "steps - sum must equal 1 !\n")
        for j in range(0, occsteps):
            ol = ""
            totalocc = 0
            for el in at.disp_occ:
                ol += f74.write([at.disp_occ[el][j]])
                totalocc += at.disp_occ[el][j]
            if printvac:
                ol += f74.write([1 - totalocc])
            output += (ol.ljust(16) 
                       + "   concentration step no. {}\n".format(j+1))
        label = i1.write([i+1])+i1.write([i+1])
        parcount += 1
        constr = "-"
        spl = [s for s in rp.searchpars if s.atom == at and s.mode == "occ"]
        if spl and spl[0].restrictTo is not None:
            sp = spl[0]
            if type(sp.restrictTo) == int:
                constr = str(sp.restrictTo)
            else:
                constr = "#"+str(rp.searchpars.index(sp.restrictTo)+1)
        elif spl and spl[0].linkedTo is not None:
            constr = "#"+str(rp.searchpars.index(spl[0].linkedTo)+1)
        info += (str(parcount).rjust(4) + ("C"+label).rjust(7) 
                 + str(at.oriN).rjust(7) + "-".rjust(7) + "occ".rjust(7)
                 + str(occsteps).rjust(7) + constr.rjust(7) + "\n")
        nsteps.append(occsteps)
    output += ("Information about start configuration: (Parameters for each "
               "place, conc for each place)\n")
    
    if rp.SEARCH_START == "control":
        controlpath = ""
        if os.path.isfile("control.chem"):
            controlpath = "control.chem"
        elif os.path.isfile(os.path.join("AUX","control.chem")):
            controlpath = os.path.join("AUX","control.chem")
            logger.warning("No control.chem file found in working folder, "
                            "using AUX/control.chem")
            rp.setHaltingLevel(1)
        else:
            logger.warning("No control.chem file found. Defaulting to random "
                          "starting configuration.")
            rp.SEARCH_START = "random"
            rp.setHaltingLevel(2)
    if rp.SEARCH_START == "control":
        # try reading control
        try:
            with open(controlpath, "r") as rf:
                controllines = rf.readlines()
        except:
            logger.error("Error reading control.chem file. Defaulting to "
                          "random starting configuration.")
            rp.SEARCH_START = "random"
            rp.setHaltingLevel(2)
        if len(controllines) < rp.SEARCH_POPULATION + 2:
            if not rp.controlChemBackup:
                logger.error("Information in control.chem is incomplete. "
                              "Defaulting to random starting configuration.")
                rp.SEARCH_START = "random"
                rp.setHaltingLevel(2)
            else:
                logger.debug("Information in control.chem is incomplete. "
                              "Loading from stored data...")
                controllines = [s+"\n" 
                                for s in rp.controlChemBackup.split("\n")]
    if rp.SEARCH_START == "random":
        output += ("  0             Certain start position (1) or random "
                   "configuration (0)") 
    else:
        output += ("  1             Certain start position (1) or random "
                   "configuration (0)\n")
        if rp.SEARCH_START == "control":
            if cull and rp.SEARCH_CULL > 0:
                if rp.SEARCH_CULL < 1:
                    ncull = int(round(rp.SEARCH_POPULATION * rp.SEARCH_CULL))
                else:
                    if rp.SEARCH_CULL < rp.SEARCH_POPULATION:
                        ncull = rp.SEARCH_CULL
                    else:
                        logger.warning("SEARCH_CULL parameter too large: "
                            "would cull entire population. Culling will be "
                            "skipped.")
                        ncull = 0
                nsurvive = rp.SEARCH_POPULATION - ncull
                clines = controllines[2:]
                csurvive = []
                if rp.SEARCH_CULL_TYPE == "genetic":  # prepare readable clines
                    try:
                        csurvive = [tl.base.readIntLine(s, width=3) 
                                    for s in clines[:nsurvive]]
                    except:
                        logger.warning("SEARCH_CULL: Failed to read old "
                                "configuration from control.chem, cannot "
                                "run genetic algorithm. Defaulting to "
                                "cloning.")
                        rp.setHaltingLevel(1)
                        csurvive = []
                for (i, line) in enumerate(clines):
                    if i < nsurvive:
                        output += line
                    elif (rp.SEARCH_CULL_TYPE == "random" or 
                          (rp.SEARCH_CULL_TYPE == "genetic" and csurvive)):
                        if rp.SEARCH_CULL_TYPE == "random":
                            nc = rp.getRandomConfig()
                        else:  # "genetic"
                            nc = rp.getOffspringConfig(csurvive)
                        ol = ""
                        for v in nc:
                            ol += i3.write([v])
                        output += ol + "\n"
                    else:    #  "cloning"
                        output += clines[random.randrange(0, nsurvive)]
            else:  # don't cull
                for (i, line) in enumerate(controllines):
                    if i > 1: # first 2 lines are empty
                        output += line
        elif rp.SEARCH_START == "centered":
            for i in range(0, rp.SEARCH_POPULATION):
                ol = ""
                for n in nsteps:
                    ol += i3.write([(n+1)/2])
                output += ol + "  1\n"
                        # !!! to check: why does this sometimes have value 2?
                        # !!! DOMAINS: Last value is domain configuration
        else:    # rp.SEARCH_START == "crandom"
            pop = [rp.getCenteredConfig()]
            for i in range(1, rp.SEARCH_POPULATION):
                pop.append(rp.getRandomConfig())
            for p in pop:
                ol = ""
                for v in p:
                    ol += i3.write([v])
                output += ol + "\n"
    # write search.steu
    try:
        with open("search.steu", "w") as wf:
            wf.write(output)
    except:
        logger.error("Failed to write search.steu file for search.")
        raise
    # write searchpars.info
    try:
        with open("searchpars.info", "w") as wf:
            wf.write(info)
    except:
        logger.error("Failed to write searchpars.info file. This file is "
                    "only for user information and will not affect operation.")
    
    if steuOnly:
        logger.debug("Wrote search.steu input for search.")
        return 0
    
    # now restrict.f
    output = ("C  This subroutine serves to restrict individual parameters to "
"""certain values inside
C  the search. These values are stored in the PARIND(IPARAM,IPOP) array.
C  Parameters are counted as listed in control file search.steu, including
C  a concentration parameter after each atomic site, and the domain weight """
"""parameters
C  for each domain at the end of the PARIND array.
C
C  perform restrictions in the following fashion: e.g.
C
C     PARIND(1,IPOP) = 5
C
C  or
C
C     PARIND(5,IPOP) = PARIND(1,IPOP)
C
C  etc.

      subroutine restrict(NPRMK,NPS,PARIND,IPOP)

      INTEGER NPRMK,NPS
      INTEGER PARIND
      DIMENSION PARIND(NPRMK,NPS)

C  begin restrictions

""")
    for (i, sp) in enumerate(rp.searchpars):
        if sp.restrictTo is None:
            continue
        if type(sp.restrictTo) == int:
            output += "      PARIND({},IPOP) = {}\n".format(i+1, sp.restrictTo)
        else:
            output += "      PARIND({},IPOP) = PARIND({},IPOP)\n".format(i+1,
                                    rp.searchpars.index(sp.restrictTo)+1)
    output += ("""
C  end restrictions

      RETURN

      END           
""")
    # write restrict.f
    try:
        with open("restrict.f", "w") as wf:
            wf.write(output)
    except:
        logger.error("Failed to write restrict.f file for search.")
        raise
    
    logger.debug("Wrote input files for search: PARAM, search.steu, "
                 "restrict.f")
    return 0

def generateDeltaInput(atom, targetel, sl, rp, deltaBasic = "", auxbeams = "",
                       phaseshifts = ""):
    """Generates a PARAM file and input for one element of one atom. 
    If deltaBasic is not passed, will call generateDeltaBasic directly. 
    Returns the delta input as string, a shortened version of that input for 
    logging as string, and the contents of the required PARAM file as 
    string."""
    if deltaBasic == "":
        deltaBasic = generateDeltaBasic()
    if auxbeams == "":
        # if AUXBEAMS is not in work folder, check AUX folder 
        if not os.path.isfile(os.path.join(".","AUXBEAMS")):
            if os.path.isfile(os.path.join(".","AUX","AUXBEAMS")):
                try:
                    shutil.copy2(os.path.join(".","AUX","AUXBEAMS"), 
                                 "AUXBEAMS")
                except:
                    logger.warning("Failed to copy AUXBEAMS from AUX folder")
            else:
                logger.warning("generateDeltaInput: AUXBEAMS not found")
        try:
            with open("AUXBEAMS", "r") as rf:
                auxbeams = rf.read()
            if auxbeams[-1] != "\n":
                auxbeams += "\n"
        except:
            logger.error("generateDeltaInput: Could not read AUXBEAMS")
            raise
    if phaseshifts == "":
        try:
            with open("_PHASESHIFTS", "r") as rf:
                phaseshifts = rf.read()
            if phaseshifts[-1] != "\n":
                phaseshifts += "\n"
        except:
            logger.error("generateDeltaInput: Could not read _PHASESHIFTS")
            raise
    MLMAX = [19, 126, 498, 1463, 3549, 7534, 14484, 25821, 43351, 69322, 
             106470, 158067, 227969, 320664, 441320]
    try:
        beamlist, _, _ = tl.writeAUXBEAMS(ivbeams=rp.ivbeams, 
                                        beamlist=rp.beamlist, write=False)
    except:
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
        geolist = [(0.0,0.0,0.0)]
    elif targetel in atom.disp_geo:
        geolist = atom.disp_geo[targetel]
    else:
        geolist = atom.disp_geo["all"]
    geosteps = len(geolist)
    ol = i4.write([geosteps])
    din += ol.ljust(29) + "NCSTEP - number of displaced positions\n"
    for disp in geolist:
        ol = f74x3.write([disp[2],disp[0],disp[1]])
        din += ol.ljust(29)+"CDISP(z,x,y) - z pointing towards bulk\n"
    din += (
"""-------------------------------------------------------------------
--- vibrational displacements of atomic site in question        ---
-------------------------------------------------------------------
""")
    if targetel == "vac":
        viblist = [0.0]
    elif targetel in atom.disp_vib:
        viblist = atom.disp_vib[targetel]
    else:
        viblist = atom.disp_vib["all"]
    vibsteps = len(viblist)
    ol = i4.write([vibsteps])
    din += ol.ljust(29) + "NDEB - number of vib. amplitudes to be considered\n"
    f74 = ff.FortranRecordWriter("F7.4")
    for disp in viblist:
        # "default" vibamp + offset, not just offset
        if targetel.lower() == "vac":
            vibamp = 0.0
        else:
            vibamp = atom.site.vibamp[targetel] + disp
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
              "69322 106470 158067 227969 320664 441320\n"
    """C  MNPSI: number of phase shift values tabulated in phase shift file
C  MNEL : number of elements for which phase shifts are tabulated
C  MNT0 : number of beams for which delta amplitude calculation is required
C  MNATOMS: currently must be set to 1. In principle number of different atomic
C      positions in a superlattice wrt the reference periodicity when computing
C      TLEED beams for a superlattice not present in the reference structure
C  MNDEB: number of thermal variation steps to be performed (outer var. loop)
C  MNCSTEP: number of geometric variation steps to be performed """
            +"(inner var. loop)\n\n")
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
    ucsurf = np.transpose(sl.ucell[:2,:2])
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    ucbulk = np.transpose(sl.bulkslab.ucell[:2,:2])
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

def writePARAM(sl, rp):
    """Writes PARAM file for the refcalc"""
    try:
        beamlist,beamblocks,beamN=tl.writeAUXBEAMS(ivbeams=rp.ivbeams, 
                                        beamlist=rp.beamlist, write=False)
    except:
        logger.error("generatePARAM: Exception while getting data from "
                      "writeAUXBEAMS")
        raise
        
    # define Clebsh-Gordon coefficient tables:
    mnlmo = [1, 70, 264, 759, 1820, 3836, 7344, 13053, 21868, 34914, 53560, 
             79443, 114492, 160952, 221408] 
    mnlm = [1, 76, 284, 809, 1925, 4032, 7680, 13593, 22693, 36124, 55276, 
            81809, 117677, 165152, 226848]
    
    # start generating output
    output = ('C  Dimension statements for Tensor LEED reference calculation, '
                '\nC  version v1.2\n\n')
    output += 'C  1. lattice symmetry\n\n'
    m = rp.SUPERLATTICE.copy()
    if m[1,1] != 0:     # m[1] not parallel to a_bulk
        if m[0,1] != 0: # m[0] not parallel to a_bulk
            # find basis in which m[0] is parallel to a_bulk
            f = tl.base.lcm(abs(int(m[0,1])),abs(int(m[1,1])))
            m[0] *= f/m[0,1]
            m[1] *= f/m[1,1]
            m[0] -= m[1]*np.sign(m[0,1])*np.sign(m[1,1])
        nl1 = abs(int(round(m[0,0])))
        nl2 = abs(int(round(abs(np.linalg.det(rp.SUPERLATTICE))/nl1)))
    else:               # m[1] already parallel to a_bulk
        nl2 = abs(int(round(m[1,0])))
        nl1 = abs(int(round(abs(np.linalg.det(rp.SUPERLATTICE))/nl2)))
    ideg = 2  # any 2D point grid is at least 2fold symmetric
    # if sl.planegroup in ['p2','pmm','pmg','pgg','cmm','rcmm']:
    #     ideg = 2
    if sl.planegroup in ['p3','p3m1','p31m']:
        ideg = 3
    elif sl.planegroup in ['p4','p4m','p4g']:
        ideg = 4
    elif sl.planegroup in ['p6','p6m']:
        ideg = 3    # should be 6, but according to TensErLEED fortran 
                    #   comments, 3 works better
    output += ('      PARAMETER (MIDEG='+str(ideg)+',MNL1='+str(nl1)
                +',MNL2='+str(nl2)+')\n')
    output += '      PARAMETER (MNL = MNL1*MNL2)\n'
    output += '\nC  2. General calculational quantities\n\n'
    output += '      PARAMETER (MKNBS = '+str(beamblocks)+')\n'
    output += '      PARAMETER (MKNT =  '+str(beamN)+')\n'
    output += ('      PARAMETER (MNPUN = '+str(len(beamlist))+', MNT0 = '
                +str(len(beamlist))+')\n')
    output += ('      PARAMETER (MNPSI = '+str(len(rp.phaseshifts))+', MNEL = '
                +str(len(rp.phaseshifts[0][1]))+')\n')
    output += '      PARAMETER (MLMAX = '+str(rp.LMAX)+')\n'
    output += ('      PARAMETER (MNLMO = '+str(mnlmo[rp.LMAX-1])+', MNLM = '
                +str(mnlm[rp.LMAX-1])+')\n')
    output += '\nC  3. Parameters for (3D) geometry within (2D) unit mesh\n\n'
    output += '      PARAMETER (MNSITE  = '+str(len(sl.sitelist))+')\n'
    output += '      PARAMETER (MNLTYPE = '+str(len(sl.layers))+')\n'
    mnbrav = 0
    mnsub = 0
    mnstack = 0
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    for layer in [l for l in sl.layers if not l.isBulk]:
        mnstack += 1
        if len(layer.atlist) == 1:
            mnbrav += 1
        if len(layer.atlist) > mnsub:
            mnsub = len(layer.atlist)
    for i, layer in enumerate([l for l in sl.layers if l.isBulk]):
        if len(sl.bulkslab.layers[i].atlist) == 1:
            mnbrav += 1
        if len(sl.bulkslab.layers[i].atlist) > mnsub:
            mnsub = len(layer.atlist)
    output += '      PARAMETER (MNBRAV  = '+str(mnbrav)+')\n'
    output += '      PARAMETER (MNSUB   = '+str(mnsub)+')\n'
    output += '      PARAMETER (MNSTACK = '+str(mnstack)+')\n'
    output += ('\nC  4. some derived quantities that must be treated '
               'explicitly (dummy dimensions for\n')
    output += 'C     special cases necessary\n\n'
    output += '      PARAMETER (MLMAX1=MLMAX+1)\n'
    output += '      PARAMETER (MLMMAX = MLMAX1*MLMAX1)\n\n'
    output += ('      PARAMETER (MNBRAV2 = '+('MNBRAV' if mnbrav > 0 
                                              else '1')+')\n\n')
    output += ('      PARAMETER (MNCOMP= '+('MNLTYPE-MNBRAV' if mnsub > 1 
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLMT  = '+('MNSUB*MLMMAX' if mnsub > 1 
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MNSUB2= '+('MNSUB * (MNSUB-1)/2' if mnsub > 1 
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLMG  = '+('MNSUB2*MLMMAX*2' if mnsub > 1 
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLMN  = '+('MNSUB * MLMMAX' if mnsub > 1 
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLM2N = '+('2*MLMN' if mnsub > 1 
                                            else '1 ')+')\n')
    output += ('      PARAMETER (MLMNI = '+('MNSUB*MLMMAX' if mnsub > 1 
                                            else '1 ')+')\n')
    try:
        with open('PARAM', 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write PARAM file")
        raise
    logger.debug("Wrote to PARAM successfully.")
    return 0

def writeAUXLATGEO(sl, rp):
    """Writes AUXLATGEO, which is part of the input FIN for the refcalc."""
    output = ''
    output += rp.systemName+' '+rp.timestamp+'\n'
                        #this originally contained V0i, but only for user info
    f72x3 = ff.FortranRecordWriter('3F7.2')
    ens = [rp.THEO_ENERGIES[0], rp.THEO_ENERGIES[1]+0.01, rp.THEO_ENERGIES[2]]
    ol = f72x3.write(ens).ljust(24)
    output += ol + 'EI,EF,DE\n'
    f74x2 = ff.FortranRecordWriter('2F7.4')
    ucsurf = np.transpose(sl.ucell[:2,:2])
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    ucbulk = np.transpose(sl.bulkslab.ucell[:2,:2])
    ol = f74x2.write(ucbulk[0]).ljust(24)
    output += ol + 'ARA1\n'
    ol = f74x2.write(ucbulk[1]).ljust(24)
    output += ol + 'ARA2\n'
    output += ' 0.0    0.0             SS1\n'
    output += ' 0.0    0.0             SS2\n'
    output += ' 0.0    0.0             SS3\n'
    output += ' 0.0    0.0             SS4\n'
    ol = f74x2.write(ucsurf[0]).ljust(24)
    output += ol + 'ARB1\n'
    ol = f74x2.write(ucsurf[1]).ljust(24)
    output += ol + 'ARB2\n'
    output += ' 0.0    0.0             SO1\n'
    output += ' 0.0    0.0             SO2\n'
    output += ' 0.0    0.0             SO3\n'
    ol = f74x2.write([0.5, rp.V0_Z_ONSET])
    ol = ol.ljust(24)
    output += ol + 'FR ASE\n'
    
    try:
        with open('AUXLATGEO', 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write AUXLATGEO file")
        raise
    logger.debug("Wrote to AUXLATGEO successfully.")
    return 0

def writeAUXNONSTRUCT(sl, rp):
    """Writes AUXNONSTRUCT, which is part of the input FIN for the refcalc."""
    try:
        beamnums,_,_=tl.writeAUXBEAMS(ivbeams=rp.ivbeams, beamlist=rp.beamlist, 
                                      write=False)
    except:
        logger.error("generatePARAM: Exception while getting data from "
                      "writeAUXBEAMS")
        raise
    
    # start generating output
    output = ''
    f74 = ff.FortranRecordWriter('F7.4')
    ol = f74.write([rp.ATTENUATION_EPS])
    output += ol+'           >>>>> ! <<<<<              TST\n'
    i4x15 = ff.FortranRecordWriter('15I4')
    ol = i4x15.write(beamnums)
    output += ol+'\n'
    f72f61 = ff.FortranRecordWriter('(F7.2, F6.1)')
    ol = f72f61.write([rp.THETA, rp.PHI]).ljust(45)
    output += ol+'THETA FI\n'
    ol = f74.write([rp.BULKDOUBLING_EPS]).ljust(45)
    output += ol+'EPS\n'
    i3 = ff.FortranRecordWriter('I3')
    ol = i3.write([rp.BULKDOUBLING_MAX]).ljust(45)
    output += ol+'LITER\n'
    ol = i3.write([rp.LMAX]).ljust(45)
    output += ol+'LMAX\n'
    
    try:
        with open('AUXNONSTRUCT', 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write AUXNONSTRUCT file")
        raise
    logger.debug("Wrote to AUXNONSTRUCT successfully.")
    return 0

def writeAUXBEAMS(ivbeams=[], beamlist=[], beamsfile='IVBEAMS', 
                  readfile='_BEAMLIST', writefile='AUXBEAMS', write=True):
    """"Reads from a _BEAMLIST file (full list of beams for calculation), 
    finds the beams listed in IVBEAMS (if not passed as a list 'beams', will 
    attempt to call readIVBEAMS directly) and copies the corresponding lines 
    to AUXBEAMS. Returns a list of the corresponding beam numbers in _BEAMLIST.
    """
    if not readfile == "_BEAMLIST":
        blstr = "Beam list (filename "+readfile+")"
    else:
        blstr = "_BEAMLIST"
        
    if ivbeams == []:         #if 'ivbeams' is empty, try to fill it
        ivbeams = tl.readIVBEAMS(beamsfile)
        
    output = '   1               IFORM\n'      # !!! WHAT IS THIS VALUE? WHERE TO GET IT FROM?
    
    # read BEAMLIST
    if beamlist == []:
        logger.warning("writeAUXBEAMS routine: no beamlist passed, "
            "attempting to read _BEAMLIST directly.")
        try:
            beamlist = tl.readBEAMLIST(readfile)
        except:
            logger.error("Error getting beamlist not found.")
            raise
    
    err = 1e-4           #since beams are saved as floats, give some error 
                         #  tolerance when comparing e.g. 1/3
    #read BEAMLIST - very little error handling here, I'm assuming 
    #  the beamlists are safe; could be added later
    foundbeams = []
    numlist = []
    blocks = 0
    totalbeams = 0
    for line in beamlist:
        llist = line.split()
        if len(llist) > 1:
            totalbeams += 1
            for b in ivbeams:
                if b.isEqual_hk((float(llist[0]), float(llist[1])), eps=err):
                    output += line
                    foundbeams.append(b)
                    numlist.append(int(line.split('.')[-1]))
        else:
            blocks += 1
    for beam in [b for b in ivbeams if b not in foundbeams]:
        logger.warning('IVBEAMS contains beam '+ beam.label +', which was '
                        'not found in '+blstr)
    #write output
    if write:
        if not writefile == 'AUXBEAMS':
            wfstr = 'AUXBEAMS file (filename '+writefile+')'
        else:
            wfstr = 'AUXBEAMS'
        try:
            with open(writefile, 'w') as wf:
                wf.write(output)
        except:
            logger.error("Failed to write "+wfstr+" file")
            raise
        logger.debug("Wrote to "+wfstr+" successfully.")
    return numlist, blocks, totalbeams

def writeAUXGEO(sl, rp):
    """Writes AUXGEO, which is part of the input FIN for the refcalc."""
    output = ''
    output += ('---------------------------------------------------------'
               '----------\n')
    output += ('--- define chem. and vib. properties for different atomic '
               'sites ---\n')
    output += ('---------------------------------------------------------'
               '----------\n')
    i3 = ff.FortranRecordWriter('I3')
    ol = i3.write([len(sl.sitelist)])
    ol = ol.ljust(26)
    output += ol + 'NSITE: number of different site types\n'
    f74x2 = ff.FortranRecordWriter('2F7.4')
    for i, site in enumerate(sl.sitelist):
        output += '-   site type '+str(i+1)+' ---\n'
        
        for el in sl.elements:
            # this reproduces the order of blocks contained in _PHASESHIFTS:
            if el in rp.ELEMENT_MIX:
                chemelList = rp.ELEMENT_MIX[el]
            else:
                chemelList = [el]
            siteList = [s for s in sl.sitelist if s.el == el]
            for cel in chemelList:
                for s in siteList:
                    if s.isEquivalent(site):
                        occ, vib = site.occ[cel], site.vibamp[cel]
                        comment = ('Occ & VibAmp for '+cel+' in '+site.label
                                   +' site')
                    else:
                        occ, vib = 0.0, 0.0
                        comment = ''
                    try:
                        ol = f74x2.write([occ, vib])
                    except:
                        logger.error("Exception while trying to write \
                            occupation / vibrational amplitude for site "
                            +site.label, exc_info=True)
                    ol = ol.ljust(26)
                    output += ol + comment + '\n'
        
    output += ('-----------------------------------------------------'
               '--------------\n')
    output += ('--- define different layer types                     '
               '           ---\n')
    output += ('-----------------------------------------------------'
               '--------------\n')
    ol = i3.write([len(sl.layers)])
    ol = ol.ljust(26)
    output += ol + 'NLTYPE: number of different layer types\n'
    f74x3 = ff.FortranRecordWriter('3F7.4')
    blayers = [l for l in sl.layers if l.isBulk]
    nblayers = [l for l in sl.layers if not l.isBulk]
    layerOffsets = [np.array([0.,0.,0.]) for i in range(0, len(sl.layers)+1)]
    if sl.bulkslab is None:
        sl.bulkslab = sl.makeBulkSlab(rp)
    for i, layer in enumerate(sl.layers):
        output += '-   layer type '+str(i+1)+' ---\n'
        if layer.isBulk:
            output += ('  2                       LAY = 2: layer type no. '
                       +str(i+1)+' has bulk lateral periodicity\n')
        else:
            output += ('  1                       LAY = 1: layer type no. '
                       +str(i+1)+' has overlayer lateral periodicity\n')
        if layer.isBulk:
            bl = sl.bulkslab.layers[blayers.index(layer)]
            bulknums = [at.oriN for at in bl.atlist]
            bulkUnique = [at for at in layer.atlist if at.oriN in bulknums]
            natoms = len(bulkUnique)
            #sanity check: ratio of unit cell areas (given simply by 
            #  SUPERLATTICE) should match ratio of written vs skipped atoms:
            arearatio = 1/abs(np.linalg.det(rp.SUPERLATTICE))
            atomratio = len(bulkUnique)/len(layer.atlist)
            if abs(arearatio-atomratio) > 1e-3:
                logger.warning('Ratio of bulk atoms inside/outside the bulk '
                    'unit cell differs from bulk/slab unit cell size ratio. '
                    'This means that the actual periodicity of the POSCAR '
                    'does not match the periodicity given in the SUPERLATTICE '
                    'parameter. Check SUPERLATTICE parameter and bulk '
                    'symmetry!')
                rp.setHaltingLevel(2)
        else:
            natoms = len(layer.atlist)
        ol = i3.write([natoms])
        ol = ol.ljust(26)
        output += ol+'number of Bravais sublayers in layer '+str(i+1)+'\n'
        if layer.isBulk:
            writelist = bulkUnique
        else:
            writelist = layer.atlist
        writelist.sort(key=lambda atom: -atom.pos[2])
        for atom in writelist:
            if not rp.LAYER_STACK_VERTICAL:
                writepos = atom.posInLayer
            else:
                writepos = atom.cartpos - np.array([nblayers[-1].cartori[0], 
                                                    nblayers[-1].cartori[1], 
                                                    atom.layer.cartori[2]])
            ol = i3.write([sl.sitelist.index(atom.site)+1])
            if natoms != 1:
                ol += f74x3.write([writepos[2], writepos[0], writepos[1]])
            else:
                # Bravais layers need to have coordinate (0.0, 0.0, 0.0)
                #  -> store actual position for later, it will go into the 
                #  interlayer vector
                ol += f74x3.write([0.0, 0.0, 0.0])
                layerOffsets[layer.num] += writepos
                layerOffsets[layer.num+1] -= writepos
            ol = ol.ljust(26)
            output += ol+'Atom N='+str(atom.oriN)+' ('+atom.el+')\n'
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('--- define bulk stacking sequence                             '
               '  ---\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('  0                       TSLAB = 0: compute bulk using layer '
               'doubling\n') #change if this is supposed to be a parameter...
    
    #########################################################################
    
    if rp.BULK_REPEAT is None:
        # assume that interlayer vector from bottom non-bulk to top bulk layer
        #   is the same as between bulk units
        # save BULK_REPEAT value for later runs, in case atom above moves
        if rp.N_BULK_LAYERS == 2:
            rp.BULK_REPEAT = (blayers[1].cartbotz 
                              - sl.layers[blayers[0].num-1].cartbotz)
        else:
            rp.BULK_REPEAT = (blayers[0].cartbotz 
                              - sl.layers[blayers[0].num-1].cartbotz)
        tl.modifyPARAMETERS(rp, "BULK_REPEAT", "{:.4f}".format(rp.BULK_REPEAT),
                            comment = "Keeps bulk spacing constant during "
                                      "search")
        logger.warning("The BULK_REPEAT parameter was undefined, which may "
                "lead to unintended changes in the bulk unit cell during "
                "optimization if the lowest non-bulk atom moves.\n"
                "# The BULK_REPEAT value determined from the POSCAR was "
                "written to the PARAMETERS file.")
        
    if type(rp.BULK_REPEAT) == np.ndarray:
        bulkc = rp.BULK_REPEAT
        if bulkc[2] < 0:
            bulkc = -bulkc
        if rp.N_BULK_LAYERS == 2:
            zdiff = bulkc[2] - (blayers[1].cartbotz - blayers[0].cartbotz)
            asaZ = bulkc[2] - (blayers[1].cartbotz - blayers[0].cartori[2])
        else:
            zdiff = bulkc[2]
            asaZ = bulkc[2] - (blayers[0].cartbotz - blayers[0].cartori[2])
        bulkc = bulkc * zdiff/bulkc[2]
        bulkc[2] = asaZ
        bvectors_ASA = bulkc
    else:
        if rp.N_BULK_LAYERS == 2:
            zdiff = rp.BULK_REPEAT - (blayers[1].cartbotz 
                                      - blayers[0].cartbotz)
            asaZ = rp.BULK_REPEAT - blayers[1].cartbotz + blayers[0].cartori[2]
        else:
            zdiff = rp.BULK_REPEAT
            asaZ = rp.BULK_REPEAT - blayers[0].cartbotz + blayers[0].cartori[2]
        # zdiff = rp.ASAZ + blayers[0].cartbotz - blayers[0].cartori[2]
        bvectors_ASA = [-sl.ucell[0][2] * zdiff/sl.ucell[2][2], 
                        -sl.ucell[1][2] * zdiff/sl.ucell[2][2], 
                        asaZ]

    # determine ASBULK - interlayer vector between bulk layers
    if rp.N_BULK_LAYERS == 2:
        # add layerOffsets for Bravais layers:
        bvectors_ASA += layerOffsets[blayers[0].num+1]
        # calculate ASBULK:
        bvectors_ASBULK = blayers[1].cartori - blayers[0].cartori
        bvectors_ASBULK[2] = blayers[1].cartori[2] - blayers[0].cartbotz
        bl2num = blayers[1].num
        # add layerOffsets for Bravais layers:
        bvectors_ASBULK -= layerOffsets[blayers[0].num+1]
        # bvectors_ASBULK += (layerOffsets[blayers[1].num+1] 
        #                     + layerOffsets[blayers[0].num])
    else:
        bl2num = blayers[0].num
        bvectors_ASBULK = bvectors_ASA
    
    ol = f74x3.write([bvectors_ASA[2], bvectors_ASA[0], bvectors_ASA[1]])
    ol = ol.ljust(26)
    output += ol + 'ASA interlayer vector between different bulk units\n'
    ol = i3.write([blayers[0].num+1])
    ol = ol.ljust(26)
    output += (ol + 'top layer of bulk unit: layer type '+str(blayers[0].num+1)
                +'\n')
    ol = i3.write([bl2num+1])
    ol = ol.ljust(26)
    output += ol + 'bottom layer of bulk unit: layer type '+str(bl2num+1)+'\n'
    ol = f74x3.write([bvectors_ASBULK[2], bvectors_ASBULK[0], 
                      bvectors_ASBULK[1]])
    ol = ol.ljust(26)
    output += (ol + 'ASBULK between the two bulk unit layers (may differ from '
                     'ASA)\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('--- define layer stacking sequence and Tensor LEED output     '
               '  ---\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    nonbulk = len(sl.layers)-rp.N_BULK_LAYERS
    if len(rp.TENSOR_OUTPUT) < nonbulk:    #check TENSOR_OUTPUT parameter
        if rp.TENSOR_OUTPUT:
            logger.warning('Parameters TENSOR_OUTPUT is defined, but contains '
                'less values than there are non-bulk layers. Missing values '
                'will be set to 1.')
            rp.setHaltingLevel(1)
        for i in range(0,nonbulk-len(rp.TENSOR_OUTPUT)):
            rp.TENSOR_OUTPUT.append(1)
    if len(rp.TENSOR_OUTPUT) > nonbulk:
        logger.warning('Parameters TENSOR_OUTPUT is defined, but contains '
            'more values than there are non-bulk layers. Excess values will '
            'be ignored.')
        rp.setHaltingLevel(1)
    ol = i3.write([len(sl.layers)-rp.N_BULK_LAYERS])
    ol = ol.ljust(26)
    output += ol + 'NSTACK: number of layers stacked onto bulk\n'
    for layer in list(reversed(nblayers)):
        n = layer.num+1
        if layer == nblayers[-1] or not rp.LAYER_STACK_VERTICAL:
            v = sl.layers[n].cartori-layer.cartori
        else:
            v = np.array([0.0,0.0,0.0])
        v[2] = sl.layers[n].cartori[2]-layer.cartbotz
        v = v + layerOffsets[n]   # add layerOffsets for Bravais layers
        ol = i3.write([n]) + f74x3.write([v[2],v[0],v[1]])
        ol = ol.ljust(26)
        output += (ol + 'layer '+str(n)+': layer type '+str(n)+', interlayer '
                   'vector below\n')     #every layer is also a layer type
        ol = i3.write([rp.TENSOR_OUTPUT[layer.num]])
        ol = ol.ljust(26)
        output += (ol + '0/1: Tensor output is required for this layer '
                   '(TENSOR_OUTPUT)\n')
        i = 1
        for atom in layer.atlist:
            # ol = 'T_'+atom.el+str(atom.oriN)
            ol = 'T_'+str(atom.oriN)
            ol = ol.ljust(26)
            output += (ol + 'Tensor file name, current layer, sublayer '+str(i)
                        +'\n')
            i += 1
    output += ('--------------------------------------------------------------'
               '-----\n')
    output += ('--- end geometrical input                                     '
               '  ---\n')
    output += ('--------------------------------------------------------------'
               '-----\n')
    
    try:
        with open('AUXGEO', 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write AUXGEO file")
        raise
    logger.debug("Wrote to AUXGEO successfully.")
    return 0

def writeMuftin(sl, rp):
    """Writes a muftin.f file, which will be compiled for the refcalc."""
    output = """
C  Subroutine muftin contains explicit energy dependence of inner 
C  potentials. The functional form should be left to the user entirely,
C  thus the subroutine is included in the input explicitly.

C  All quantities should be set in eV. Conversion to Hartree (atomic units)
C  is handled by ref-calc.f itself. Both the real and imaginary part should
C  be given as positive values (the program then handles the sign correctly).

      subroutine muftin(EEV,VO,VV,VPI,VPIS,VPIO)

C  global variable

C  EEV :  electron energy in the vacuum region
C  VO  :  difference between real part of inner potential in the bulk and in
C         topmost layer. Usually irrelevant -> set to zero.
C  VV  :  real part of the inner potential in the bulk.
C  VPI :  imaginary part of the inner potential.
C  VPIS:  in case different values for the imaginary part of the inner
C  VPIO:  potential for bulk and topmost layer were desired, VPIS would
C         correspond to the bulk value, VPIO would be the respective
C         value for the topmost layer. Usually, the effect is irrelevant,
C         i.e. both are set equal to VPI.

      real EEV,VO,VV,VPI,VPIS,VPIO

C  local variable

C  workfn: Work function of the LEED filament material. Theoretical predictions
C         for the inner potential usually have their energy zero fixed at the
C         Fermi level; the experimental energy scale (-> EEV) nominally
C         does the same, except for the fact that electrons need to overcome
C         the work function of the cathode first. Note that the this value is
C         thus formally determined by a LEED fit, yet don't trust its accuracy.

      real workfn

C  set work function of cathode
C  work function should be positive (added to exp. energy EEV)

      workfn = """
    output += str(round(rp.FILAMENT_WF, 4))+"\n"
    output += """
C  set real part of inner potential

"""
    oline = "      VV = "+rp.V0_REAL
    output += tl.base.fortranContLine(oline) + "\n"
    output += """
      write(6,*) workfn, EEV
      write(6,*) VV

c  set difference between bulk and overlayer inner potential

      VO = 0.

c  set imaginary part of inner potential - energy independent value used here

"""
    oline = "      VPI = "+str(rp.V0_IMAG)
    output += tl.base.fortranContLine(oline) + "\n"
    output += """
C  set substrate / overlayer imaginary part of inner potential

"""
    oline = "      VPIS = "+str(rp.V0_IMAG)
    output += tl.base.fortranContLine(oline) + "\n"
    oline = "      VPIO = "+str(rp.V0_IMAG)
    output += tl.base.fortranContLine(oline) + "\n"
    output += """
      return
      end"""
    try:
        with open("muftin.f", "w") as wf:
            wf.write(output)
        logger.debug("Wrote to muftin.f successfully.")
    except:
        logger.error("Exception while writing muftin.f file: ", 
                      exc_info=True)
        raise
    return 0

def writeIVBEAMS(sl, rp, filename="IVBEAMS"):
    """Writes an IVBEAMS file based on rp.exbeams. Returns 
    those beams in IVBEAMS form."""
    d = tl.leedbase.getLEEDdict(sl, rp)
    if d is None:
        logger.error("Failed to write IVBEAMS")
        return []
    makebeams = project_to_first_domain(d, [b.hkfrac for b in rp.expbeams])
    output = "This IVBEAMS file was automatically generated from EXPBEAMS\n"
    ivbeams = [tl.Beam(hk) for hk in makebeams]
    for b in ivbeams:
        output += "{: 10.6f} {: 10.6f}\n".format(b.hk[0], b.hk[1])
    try:
        with open(filename, "w") as wf:
            wf.write(output)
        logger.debug("Wrote IVBEAMS file successfully.")
    except:
        logger.error("Exception while writing IVBEAMS file: ", 
                      exc_info=True)
        raise
    return ivbeams

def writeCONTCAR(sl, filename='CONTCAR', reorder=False, comments='none',
                 silent=False):
    """Takes a Slab object and writes a CONTCAR based on the atlist. If a 
    POSCAR is in the folder, it will take header data from it, otherwise is 
    will fill in a generic header. Defining a filename other than CONTCAR is 
    optional. If the reorder tag is set as True, atoms will be sorted by 
    element and by z coordinate; otherwise, the original order will be used. 
    The 'comments' tag defines whether additional information for each atom 
    should be printed; 'none' writes a normal POSCAR, 'all' all annotations, 
    'nodir' all annotations except directions relating to a or b 
    (meant to be used with POSCAR_oricell), 'bulk' is like 'none' but 
    writes the space group (meant to be used with POSCAR_bulk)."""
    if reorder:
        sl.sortByZ()
        sl.sortByEl()
    else:
        sl.sortByEl()   #this is the minimum that has to happen, otherwise the 
                        #  output will be buggy
    output = '' 	# store output here temporarily before writing it
    # see if a POSCAR is there
    epos = True
    try:
        rf = open('POSCAR', 'r')
    except:
        epos = False
    #header
    if epos:
        #copy the first line:
        output += rf.readline()
        rf.close()
    else:
        #fill dummy header:
        output += 'unknown\n'
    output += '   1.0\n'    #scaling will be 1.0, since we apply the scaling 
                            #  to the unit cell
    #unit cell
    m = np.transpose(sl.ucell)
    for vec in m:
        for val in vec:
            s = '%.16f'%val
            output += '{:>22}'.format(s)
        output += '\n'
    #atom types and numbers
    for el in sl.oriels:
        output += '{:>5}'.format(el)
    output += '\n'
    for n in sl.nperelem:
        output += '{:>5}'.format(n)
    output += '\n'
    #header line 'Direct'
    ol = 'Direct'
    if comments != 'none':
        ol = ol.ljust(16)
        ol += 'Plane group = '
        if comments == 'nodir': ol += '*'
        if (sl.planegroup in ["pm", "pg", "cm", "rcm", "pmg"] 
               and comments != 'nodir'):
            ol += sl.planegroup + str(sl.orisymplane.par)
        else:
            ol += sl.planegroup
        if sl.planegroup != sl.foundplanegroup.split('[')[0]:
            ol += ' (found '
            if '[' in sl.foundplanegroup and comments == 'nodir':
                ol += sl.foundplanegroup.split('[')[0]
            else:
                ol += sl.foundplanegroup
            ol += ')'
        ol = ol.ljust(60)
        if comments == 'all':
            ol += ('N'.rjust(5)+'SiteLabel'.rjust(12)
                   +'Layer'.rjust(7)+'Linking'.rjust(9)+'FreeDir'.rjust(12))
        elif comments == 'nodir':
            ol += ('N'.rjust(5) + 'SiteLabel'.rjust(12) 
                   + 'Layer'.rjust(7) + 'Linking'.rjust(9) 
                   + 'FreeDir'.rjust(12) + ': see POSCAR')
    output += ol+'\n'
    # NperElCount = 0
    # lastEl = sl.atlist[0].el
    for i, at in enumerate(sl.atlist):
        # if lastEl != at.el:
        #     lastEl = at.el
        #     NperElCount = 1
        # else:
        #     NperElCount += 1
        ol = ''
        for coord in at.pos:
            s = '%.16f'%coord
            ol += '{:>20}'.format(s)
        if comments != 'none' and comments != 'bulk':
            ol += str(i+1).rjust(5)                         #N
            # ol += (at.el+str(NperElCount)).rjust(9)         #NperEl
            ol += at.site.label.rjust(12)                   #SiteLabel
            ol += str(at.layer.num+1).rjust(7)              #Layer
            if len(at.linklist) <= 1:                       #Linking
                ol += 'none'.rjust(9)
            else:
                ol += str((sl.linklists.index(at.linklist))+1).rjust(9)
            if comments != 'nodir':                         #FreeDir
                if at.layer.isBulk:
                    ol += 'bulk'.rjust(12)
                elif type(at.freedir) == int:
                    if at.freedir == 0:
                        ol += 'locked'.rjust(12)
                    else:
                        ol += 'free'.rjust(12)
                else:
                    ol += str(at.freedir).rjust(12)
        output += ol+'\n'
    #write output
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write "+filename)
        raise
    if not silent:
        logger.info("Wrote to "+filename+" successfully")
    return 0

def writeOUTBEAMS(beams, filename="THEOBEAMS.csv", sep="; "):
    """Takes a list of Beam objects and writes them to a comma-separated 
    file."""
    nan = "NaN" #what to put when no value
    output = "E".rjust(7)+sep
    w = max(11, beams[0].lwidth*2 + 3)
    energies = []
    for b in beams:
        output += b.label.rjust(w)+sep
        energies.extend([k for k in b.intens if k not in energies])
    output = output[:-len(sep)]
    output += "\n"
    energies.sort()
    for en in energies:
#        output += str(round(en,2)).rjust(7)+sep
        output += '{:7.2f}'.format(en) + sep
        for b in beams:
            if en in b.intens:
                output += '{:0.5E}'.format(b.intens[en]).rjust(w) + sep
            else:
                output += nan.rjust(w) + sep
        output = output[:-len(sep)]
        output += "\n"
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to "+filename+" successfully")
    return 0

def writeVIBROCC_OUT(sl, rp, filename="VIBROCC_OUT", silent=False):
    """Writes a new VIBROCC file with the optimized parameters obtained after 
    the search. The new file will not follow the ordering of the input VIBROCC 
    file."""
    output = "= Vibrational Amplitudes\n"
    for site in sl.sitelist:
        ol = site.label + " = "
        targetels = [el for el in site.vibamp if site.occ[el] != 0]
        if len(targetels) == 1:
            ol += "{:.4g}".format(site.vibamp[targetels[0]])
        else:
            for i, el in enumerate(targetels):
                ol += el + " {:.4g}".format(site.vibamp[el])
                if i < len(targetels)-1:
                    ol += ", "
        output += ol + "\n"
        
    output += "\n= Occupations\n"
    for site in sl.sitelist:
        write = True
        ol = site.label + " = "
        targetels = [el for el in site.occ if site.occ[el] != 0]
        if len(targetels) == 1:
            ol += "{:.4g}".format(site.occ[targetels[0]])
            if site.occ[targetels[0]] == 1 and site.el == targetels[0]:
                write = False
        else:
            for i, el in enumerate(targetels):
                ol += el + " {:.4g}".format(site.occ[el])
                if i < len(targetels)-1:
                    ol += ", "
        if write:
            output += ol + "\n"
    
    output += "\n= Search offsets\n"
    # figure out which atoms to write, in which order
    offsetList = []    # metric per atom
    for at in [at for at in sl.atlist if not at.layer.isBulk]:
        to = 0
        # !!! TODO: Think about weights for the three
        for el in at.offset_occ:
            to += abs(at.offset_occ[el])
        for el in at.offset_vib:
            to += abs(at.offset_vib[el])
        for el in at.offset_geo:
            to += np.linalg.norm(at.offset_geo[el])
        if to >= 1e-4: # !!! TODO: enough?
            offsetList.append((to, at))
    offsetList.sort(key=lambda tup: -tup[0])
    for (_, at) in offsetList:
        # POSITION OFFSET
        write = False
        ol = "POS {} = ".format(at.oriN)
        for el in at.offset_geo:
            if np.linalg.norm(at.offset_geo[el]) >= 1e-4:
                write = True
                ol += el + " {:g} {:g} {:g}, ".format(at.offset_geo[el][0],
                                                      at.offset_geo[el][1],
                                                      -at.offset_geo[el][2])
        if write:
            output += ol[:-2]+"\n"
        # VIBRATION OFFSET
        write = False
        ol = "VIB {} = ".format(at.oriN)
        for el in at.offset_vib:
            if abs(at.offset_vib[el]) >= 1e-3:
                write = True
                ol += el + " {:g}, ".format(at.offset_vib[el])
        if write:
            output += ol[:-2]+"\n"
        # OCCUPATION OFFSET
        write = False
        ol = "OCC {} = ".format(at.oriN)
        for el in at.offset_occ:
            if abs(at.offset_occ[el]) >= 1e-3:
                write = True
                ol += el + " {:g}, ".format(at.offset_occ[el])
        if write:
            output += ol[:-2]+"\n"
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write "+filename)
        raise
    if not silent:
        logger.info("Wrote to "+filename+" successfully")
    return 0

def writeAUXEXPBEAMS(beams, filename="AUXEXPBEAMS", header="Unknown system", 
                     write=True, numbers=False):
    """Takes a list of Beam objects and writes them in the format required by 
    TensErLEED for experimental beam data. Returns the whole output as a 
    string. 'numbers' defines whether a sequence of beam numbers in line 2 
    is expected."""
    output = header+"\n"
    if numbers:
        i3x25 = ff.FortranRecordWriter("25I3")
        output += i3x25.write([n+1 for n in range(0, len(beams))]) + "\n"
    output += " (12F6.2)\n"
    f62x12 = ff.FortranRecordWriter('12F6.2')
    i4 = ff.FortranRecordWriter('I4')
    for beam in beams:
        # renormalize
        scaling = 999.99 / max(beam.intens.values())
        for k in beam.intens:
            beam.intens[k] *= scaling
        # write
        output += "*"+beam.label.replace("|"," ") + "*\n"
        ol = i4.write([len(beam.intens)]).ljust(8)
        output += ol + '{:10.4E}\n'.format(scaling)
        outlist = [val for tup in beam.intens.items() for val in tup]
                             # zip & flatten energies and values into one list
        output += f62x12.write(outlist)
        if output[-1] != "\n":
            output += "\n"
    #write output
    if write:
        try:
            with open(filename, 'w') as wf:
                wf.write(output)
        except:
            logger.warning("Failed to write "+filename)
        logger.debug("Wrote to "+filename+" successfully")
    return output

def writePatternInfo(sl, rp, filename="PatternInfo.tlm"):
    """Writes a PatternInfo file that can be used by the TLEEDMAP GUI utility 
    to display the expected LEED pattern and show beam labelling."""
    output = "eMax = {:.2f}\n".format(rp.THEO_ENERGIES[1])
    mstring = "[[{}, {}], [{}, {}]]".format(sl.ucell[0,0], sl.ucell[1,0],
                                            sl.ucell[0,1], sl.ucell[1,1])
    output += "surfBasis = "+mstring+"\n"
    mstring = ("[[{:.0f}, {:.0f}], [{:.0f}, {:.0f}]]"
               .format(rp.SUPERLATTICE[0,0], rp.SUPERLATTICE[0,1],
                       rp.SUPERLATTICE[1,0], rp.SUPERLATTICE[1,1]))
    output += "superlattice = "+mstring+"\n"
    if sl.planegroup in ["pm", "pg", "cm", "rcm", "pmg"]:
        pgstring = sl.planegroup+str(sl.orisymplane.par)
    else:
        pgstring = sl.planegroup
    output += "surfGroup = "+pgstring+"\n"
    if sl.bulkslab is None:
        logger.error("PatternInfo.tlm: bulk slab has not been initialized.")
        return 1
    output += "bulkGroup = "+sl.bulkslab.foundplanegroup+"\n"
    output += "bulk3Dsym = "+sl.bulkslab.getBulk3Dstr()
    #write output
    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to "+filename+" successfully")
    return 0

def writeSearchOutput(sl, rp, parinds=None, silent=False):
    """Modifies data in sl and rp to reflect the search result given by 
    parinds, then writes POSCAR_OUT and VIBROCC_OUT. If no configuration 
    'parinds' is passed, will use rp.searchResultConfig instead"""
    if parinds is None:
        if rp.searchResultConfig is None:
            logger.error("Failed to write search output: No configuration "
                          "passed.")
            return 1
        else:
            parinds = rp.searchResultConfig[0]
    # If atom and site original states are not yet saved, do it now:
    for at in sl.atlist:
        at.storeOriState()
    for site in sl.sitelist:
        if site.oriState is None:
            tmp = copy.deepcopy(site)
            site.oriState = tmp

    sl.getCartesianCoordinates()
    uci = np.linalg.inv(sl.ucell)
    for at in sl.atlist:
        # make list of searchpars addressing this atom:
        sps = [sp for sp in rp.searchpars if sp.atom == at and sp.el != "vac"]
                                                # ignore vacancies
        if len(sps) > 0:
            # first deal with occupation:
            par = [sp for sp in sps if sp.mode == "occ"][0] # exactly one
            i = rp.searchpars.index(par)
            newocc = {}
            for el in at.disp_occ:
                # store as offset for now, recalculate for sites later
                at.offset_occ[el] = (at.disp_occ[el][parinds[i]-1] 
                                     - at.site.oriState.occ[el])
                newocc[el] = at.disp_occ[el][parinds[i]-1]
                                # will be used for re-calculating position
            # now vibrations:
            vibpars = [sp for sp in sps if sp.mode == "vib"]
                                        # can be empty, or one per element
            for par in vibpars:
                i = rp.searchpars.index(par)
                el = par.el
                if el not in at.disp_vib:
                    dict_el = "all"
                else: 
                    dict_el = el
                # store as offset for now, recalculate for sites later
                at.offset_vib[el] = at.disp_vib[dict_el][parinds[i]-1]
            # now position - recalculate right away:
            geopars = [sp for sp in sps if sp.mode == "geo"]
                                        # can be empty, or one per element
            if not geopars:
                continue
            totalocc = 0
            newpos = {}     # new position per element
            totalpos = np.array([0.,0.,0.])
            for par in geopars:
                i = rp.searchpars.index(par)
                el = par.el                
                if el not in at.disp_geo:
                    dict_el = "all"
                else:
                    dict_el = el
                disp = np.copy(at.disp_geo[dict_el][parinds[i]-1])
                disp[2] *= -1
                rdisp = np.dot(uci, disp)
                newpos[el] = at.oriState.pos + rdisp
                totalpos = totalpos + (newpos[el] * newocc[el])
                totalocc += newocc[el]
            if totalocc > 0:
                totalpos = totalpos / totalocc
            at.pos = totalpos
            for el in newpos:
                rel_off = newpos[el] - at.pos
                off = np.dot(sl.ucell, rel_off)
                off[2] *= -1
                at.offset_geo[el] = off
    sl.collapseFractionalCoordinates()
    sl.getCartesianCoordinates()
    sl.updateLayerCoordinates()
    # now update site occupations and vibrations:
    for site in sl.sitelist:
        siteats = [at for at in sl.atlist if at.site == site 
                   and not at.layer.isBulk]
        for el in site.occ:
            total_occ = 0
            # n_occ = 0
            total_vib = 0
            # n_vib = 0
            for at in siteats:
                if el in at.offset_occ:
                    total_occ += at.offset_occ[el]
                # n_occ += 1
                if el in at.offset_vib:
                    total_vib += at.offset_vib[el]
                # n_vib += 1
            # if n_occ > 0:
            offset_occ = total_occ/len(siteats)
            site.occ[el] = site.oriState.occ[el] + offset_occ
            # else:
            #     offset_occ = 0.0
            # if n_vib > 0:
            offset_vib = total_vib/len(siteats)
            site.vibamp[el] = site.oriState.vibamp[el] + offset_vib
            # else:
            #     offset_vib = 0.0
            for at in siteats:
                if el in at.offset_occ:
                    at.offset_occ[el] -= offset_occ
                if el in at.offset_vib:
                    at.offset_vib[el] -= offset_vib
    # write POSCAR_OUT
    # fn = "POSCAR_OUT_N={}_R={:.4f}".format(popcount[0], pops[0][0])
    fn = "POSCAR_OUT_"+rp.timestamp
    tmpslab = copy.deepcopy(sl)
    tmpslab.sortOriginal()
    try:
        writeCONTCAR(tmpslab, filename=fn, comments="all", silent=silent)
    except:
        logger.error("Exception occured while writing POSCAR_OUT: ", 
                      exc_info=True)
        rp.setHaltingLevel(2)
    # write VIRBOCC_OUT
    # fn = "VIBROCC_OUT_N={}_R={:.4f}".format(popcount[0], pops[0][0])
    fn = "VIBROCC_OUT_"+rp.timestamp
    try:
        writeVIBROCC_OUT(sl, rp, filename=fn, silent=silent)
    except:
        logger.error("Exception occured while writing VIBROCC_OUT: ", 
                      exc_info=True)
        rp.setHaltingLevel(2)
    return 0

def writeRfactorPdf(beams, colsDir='', outName='Rfactor_plots.pdf', 
                    plotcolors = None, analysisFile='', v0i=0.):
    '''
    Creates a single PDF file containing the plots of R-factors
    
    Parameters
    ----------
    beams: list of (name, R) tuples
           name: str
                 formatted fractional index of the beam
           R: float
              R-factor value
    colsDir: kwarg, str
             path to folder containing the files theo.column and exp.column
             generated by the R-factor routine of TensErLEED.
             default: current path
    outName: kwarg, str
             name of the file (with or without extension) to which the plots
             will be saved.
             default: 'Rfactor_plots.pdf'
    analysisFile: kwarg, string
             if not empty, a more extensive R-factor analysis pdf with 
             calculated Y-functions and absolute errors will be written to the 
             given file name.
    v0i: kwarg, float
             imaginary part of the inner potential for calculating Y-functions.
             Should always be passed if analysisFile is passed.
    
    Returns
    -------
    None if error, 0 if successful
    
    '''
    global plotting
    if not plotting:
        logger.debug("Necessary modules for plotting not found. Skipping "
                      "R-factor plotting.")
        return 0
    
    fnames = ['theo.column', 'exp.column']
    
    if not hasattr(beams, '__len__'):
        logger.error("writeRfactorPdf: First argument should be list, not "
                      +str(type(beams)))
        return None
    
    xxyy=[]
    for fname in fnames:
        try:
            f = open(os.path.join(colsDir, fname), 'r')
        except FileNotFoundError:
            logger.error("writeRfactorPdf: File {} not found. Aborting."
                          .format(fname))
            return None
        except PermissionError:
            logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                          .format(fname))
            return None
        
        cols = [[float(col) for col in line.split()] for line in f]
        
        if(np.shape(cols)[1] != 2*len(beams)):
            logger.error("writeRfactorPdf: Number of beams in file {} does "
                          "not match the input. Aborting.".format(fname))
            return None
        
        cols = np.array(cols)
        xy = np.split(cols, len(beams), axis=1)
        # xy is now a list of 2D arrays.
        # Each array has the form [[en1, intens1], ...]
        # 
        # for each beam, get rid of the points that have (en, intens) = (0, 0)
        # so that they don't screw up the plots later
        xy = [coords[~np.all(coords < 1e-3, axis=1)]  for coords in xy]
        xxyy.append(xy)
    
    xyTheo = xxyy[0]
    xyExp = xxyy[1]
    
    # find min and max values of x and y for plotting all curves
    # on the same horizontal scale and leaving a little y space for the legend
    xmin = min(min(xy[:, 0]) for xy in [*xyTheo, *xyExp])
    xmax = max(max(xy[:, 0]) for xy in [*xyTheo, *xyExp])
    
    ymin = min(min(xy[:, 1]) for xy in [*xyTheo, *xyExp])
    ymax = max(max(xy[:, 1]) for xy in [*xyTheo, *xyExp])
    dy = ymax - ymin
    
    # Set up stuff needed for the plots
    # (will pack two plots per each pdf page)
    
    try:
        pdf = PdfPages(outName)
    except PermissionError:
        logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                      .format(outName))
        return None
    
    figsize = (7, 7)
    # set ticks spacing to 50 eV and round the x limits to a multiple of it
    tick = 50
    tickloc = plticker.MultipleLocator(base=tick)
    xlims = (np.round((xmin - 10)/tick)*tick,
             np.round((xmax + 10)/tick)*tick)
    dx = xlims[1] - xlims[0]
    
    ylims = (ymin - 0.02*dy, ymax + 0.22*dy)
    
    # positions of the beam name and R factor text
    namePos = (xlims[0] + 0.45*dx, ylims[1] - 0.1*dy)
    rPos = (namePos[0], namePos[1]-0.085*dy)
    figs = []
    # the following will spam the logger with debug messages; disable.
    loglevel = logger.level
    logger.setLevel(logging.INFO)
    try:
        for ct, (name, rfact, theo, exp) in enumerate(zip(*zip(*beams),
                                                           xyTheo, xyExp)):
            if ct % 2 == 0:
                # need a new figure
                # squeeze=True returns the 2 axes as a 1D (instead of 2D) array
                fig, axs = plt.subplots(2, figsize=figsize, squeeze=True)
                fig.subplots_adjust(left=0.08, right=0.92,
                                    bottom=0.07, top=0.98,
                                    wspace=0, hspace=0.07)
                figs.append(fig)
                [ax.set_xlim(*xlims) for ax in axs]
                [ax.set_ylim(*ylims) for ax in axs]
                [ax.get_yaxis().set_visible(False) for ax in axs]
                [ax.tick_params(bottom=True,
                                top=True,
                                axis='x', direction='in') for ax in axs]
                [ax.xaxis.set_major_locator(tickloc) for ax in axs]
                axs[1].set_xlabel("Energy (eV)")
            
            idx = ct%2
            if plotcolors is not None:
                if not all([matplotlib.colors.is_color_like(s) 
                            for s in plotcolors]):
                    plotcolors = None
                    logger.warning("writeRfactorPdf: Specified colors not "
                        "recognized, reverting to default colors")
            if plotcolors is None:
                axs[idx].plot(theo[:, 0], theo[:, 1], label='Theoretical',)
                axs[idx].plot(exp[:, 0], exp[:, 1], label='Experimental')
            else:
                axs[idx].plot(theo[:, 0], theo[:, 1], label='Theoretical',
                              color=plotcolors[0])
                axs[idx].plot(exp[:, 0], exp[:, 1], label='Experimental',
                              color=plotcolors[1])
            axs[idx].annotate(name, namePos, fontsize=11)
            axs[idx].annotate("R = {:.4f}".format(rfact), rPos, fontsize=11)
            axs[idx].legend()
            
        # finally, in case the last figure is empty (i.e. the number of beams 
        # is odd) turn off the last axes (but leave the blank space).
        if len(beams) % 2 == 1:
            axs[1].axis('off')
            axs[0].set_xlabel("Energy (eV)")
        
        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)
    except:
        logger.error("writeRfactorPdf: Error while writing rfactor pdf: ",
                      exc_info = True)
    finally:
        pdf.close()
        logger.setLevel(loglevel)
    
    if not analysisFile:
        return 0
    
    # write R-factor analysis
    try:
        pdf = PdfPages(analysisFile)
    except PermissionError:
        logger.error("writeRfactorPdf: Cannot open file {}. Aborting."
                      .format(analysisFile))
        return None
    
    figsize = (5.8, 8.3)
    figs = []
    # the following will spam the logger with debug messages; disable.
    loglevel = logger.level
    logger.setLevel(logging.INFO)
    try:
        for i, (name, rfact, theo, exp) in enumerate(zip(*zip(*beams),
                                                           xyTheo, xyExp)):
            fig, axs = plt.subplots(3, figsize=figsize, 
                                                   squeeze=True)
            fig.subplots_adjust(left=0.06, right=0.94,
                                bottom=0.07, top=0.98,
                                wspace=0, hspace=0.08)
            figs.append(fig)
            [ax.set_xlim(*xlims) for ax in axs]
            axs[0].set_ylim(*ylims)
            [ax.get_yaxis().set_ticks([]) for ax in axs]
            [ax.tick_params(bottom=True,
                            top=True,
                            axis='x', direction='in') for ax in axs]
            [ax.xaxis.set_major_locator(tickloc) for ax in axs]
            axs[0].set_ylabel("Intensity (arb. units)")
            axs[1].set_ylabel("Y")
            axs[2].set_ylabel("\u2211(\u0394Y)\u00b2")    # sum delta Y ^2
            axs[2].set_xlabel("Energy (eV)")
            
            ytheo = tl.leedbase.getYfunc(theo, v0i)
            yexp = tl.leedbase.getYfunc(exp, v0i)
            dy = np.array([(ytheo[j, 0], yexp[j, 1] - ytheo[j, 1]) 
                           for j in range(0, len(ytheo))])
            dysq = np.copy(dy)
            dysq[:,1] = dysq[:,1]**2
            idysq = np.array([dysq[0]])
            for j in range(1, len(dysq)):
                idysq = np.append(idysq, [[dysq[j,0], idysq[j-1,1]+dysq[j,1]]], 
                                  axis=0)
            
            axs[1].plot(xlims, [0., 0.], color='grey', alpha=0.2)
            if plotcolors is not None:
                if not all([matplotlib.colors.is_color_like(s) 
                            for s in plotcolors]):
                    plotcolors = None
                    logger.warning("writeRfactorPdf: Specified colors not "
                        "recognized, reverting to default colors")
            if plotcolors is None:
                axs[0].plot(theo[:, 0], theo[:, 1], label='Theoretical')
                axs[0].plot(exp[:, 0], exp[:, 1], label='Experimental')
                axs[1].plot(ytheo[:, 0], ytheo[:, 1], label='Theoretical')
                axs[1].plot(yexp[:, 0], yexp[:, 1], label="Experimental")
            else:
                axs[0].plot(theo[:, 0], theo[:, 1], label='Theoretical',
                              color=plotcolors[0])
                axs[0].plot(exp[:, 0], exp[:, 1], label='Experimental',
                              color=plotcolors[1])
                axs[1].plot(ytheo[:, 0], ytheo[:, 1], label='Theoretical',
                            color=plotcolors[0], linewidth=0.75)
                axs[1].plot(yexp[:, 0], yexp[:, 1], label="Experimental",
                            color=plotcolors[0], linewidth=0.75)
            axs[1].plot(dy[:, 0], dy[:, 1], label="\u0394Y", color="black",
                        linewidth=0.5)
            axs[2].plot(idysq[:, 0], idysq[:, 1], color="black", 
                        drawstyle="steps-mid")
            
            axs[0].annotate(name, namePos, fontsize=10)
            axs[0].annotate("R = {:.4f}".format(rfact), rPos, fontsize=10)
            axs[0].legend()
            axs[1].legend()
            
        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)
    except:
        logger.error("writeRfactorPdf: Error while writing analysis pdf: ",
                      exc_info = True)
    finally:
        pdf.close()
        logger.setLevel(loglevel)
    return 0

def writeSearchProgressPdf(rp, gens, rfacs, lastconfig, 
                           outname = "Search-progress.pdf",
                           csvname = "Search-progress.csv",
                           markers = []):
    global plotting
    if not plotting:
        return 0
    
    figsPerPage = 5
    parsPerFig = 8
    
    searchname = rp.disp_blocks[rp.search_index][1]
    
    rfacsMin = np.array([min(rfa) for rfa in rfacs])
    rfacsMax = np.array([max(rfa) for rfa in rfacs])
    rfacsMean = np.array([np.mean(rfa) for rfa in rfacs])
    rp.searchplots[-1] = (searchname, gens, rfacsMin, rfacsMax, rfacsMean)
    # rfacsStd = np.array([np.std(rfa, ddof = 1) for rfa in rfacs])
    rlastunique = [rfacs[-1][0]]
    lastpops = [1]
    for r in rfacs[-1][1:]:
        if r == rlastunique[-1]:
            lastpops[-1] += 1
        else:
            rlastunique.append(r)
            lastpops.append(1)
    allcolors = []
    if rlastunique[-1] == rlastunique[0]:
        colors = [(0.,0.,0.,1.)] * len(rlastunique)  #black
        allcolors = [(0.,0.,0.,1.)] * len(rfacs[-1])
    else:
        colors = []
        for (i, r) in enumerate(rlastunique):
            w = (r - rlastunique[0]) / (rlastunique[-1] - rlastunique[0])
            w = np.sqrt(w)   # stronger scaling towards red
            colors.append((w, 0., 0.))
            for j in range(0,lastpops[i]):
                allcolors.append((w, 0., 0., 1.))
    if (not rp.rfacscatter_all) or (rp.rfacscatter_all[-1][0] != gens[-1]):
        rp.storeRfacScatter([gens[-1]]*len(rlastunique), rlastunique, 
                            lastpops, colors)
    deltagens = [gens[n] - gens[n-1] for n in range(1, len(gens))]
    maxdgx, maxdgy = [gens[1]], [deltagens[0]]
    for i in range(0, len(deltagens)):
        if deltagens[i] >= maxdgy[-1]:
            maxdgx.append(gens[i+1])
            maxdgy.append(deltagens[i])

    sps = [sp for sp in rp.searchpars if sp.el != "vac" and sp.steps > 1]

    figsize = (5.8, 8.3)
    figs = []
    
    # CSV output
    sep = "; "
    width = 12
    titles = ["Generation", "Gen_Delta", "R_min", "R_max", "R_mean"]
    output = ""
    for t in titles:
        output += t.rjust(width)+sep
    output = output[:-len(sep)] + "\n"
    mc = markers[:]
    for i in range(0, len(gens)):
        if mc:
            if gens[i] > mc[0][0]:
                output += mc[0][1] + "\n"   # comment line for marker
                mc.pop(0)
        for l in [gens, [0]+deltagens]:
            output += str(l[i]).rjust(width)+sep
        for l in [rfacsMin, rfacsMax, rfacsMean]:
            output += "{:.4f}".format(l[i]).rjust(width)+sep
        output = output[:-len(sep)] + "\n"
    try:
        with open(csvname, "w") as wf:
            wf.write(output)
    except:
        logger.warning("Failed to write "+csvname)
    
    # the following will spam the logger with debug messages; disable.
    loglevel = logger.level
    logger.setLevel(logging.INFO)
    try:
        # R-FACTOR AND GENERATION DELTA
        # create figure
        fig, (rfp, dgp) = plt.subplots(2, 1, sharex=True, figsize = figsize)
        dgp.set_xlabel('Generations')
        rfp.set_ylabel('R-Factor')
        dgp.set_ylabel('Generation delta')
        #plot markers
        part = 0
        for (i, g) in enumerate(gens):
            if g > gens[-1] * 0.2:
                part = len(gens) - i
                break
        part = max(100, part)
        rfmin, rfmax = min(rfacsMin[-part:]), max(rfacsMax[-part:])
        rYrange = [rfmin-(rfmax-rfmin)*0.1, rfmax+(rfmax-rfmin)*0.1]
        
        labely = rYrange[0] + (rYrange[1]-rYrange[0])*0.99
        xoff = gens[-1]*0.005
        for (xpos, label) in markers:
            rfp.axvline(x = xpos, lw = 0.5, c = "black")
            dgp.axvline(x = xpos, lw = 0.5, c = "black")
            rfp.text(xpos + xoff, labely, label, rotation=-90, 
                     verticalalignment="top", size=4)
        # plot data
        rfp.plot(gens, rfacsMin, '-', color='black', label = "Best")
        rfp.fill_between(gens, rfacsMin, rfacsMax, facecolor='grey', 
                         alpha=0.2, label="Range")
        rfp.plot(gens, rfacsMean, label = "Mean")
        # rfp.fill_between(gens, rfacsMean - rfacsStd/2, rfacsMean + rfacsStd/2, 
        #                  alpha=0.2, label="Sigma")
        (x,y,s,c) = list(zip(*rp.rfacscatter))
        rfp.scatter(x, y, s=s, c=c)
        # rfp.scatter([gens[-1]]*len(rlastunique), rlastunique, 
        #             s = lastpops, c = colors)
        dgp.scatter(gens[1:], deltagens, s = 4, c='black', label="Every")
        dgp.plot(maxdgx, maxdgy, '-', color='royalblue', label="Max")
        # layout
        rfp.set_ylim(rYrange)
        dgp.set_ylim([0, max(deltagens)*1.1])
        dgp.ticklabel_format(axis="x", style="sci", scilimits=(0,4))
        fig.tight_layout()
        rfp.legend(loc="lower left")
        dgp.legend(loc="upper left")
        figs.append(fig)
        
        # SEARCH PARAMETER SCATTER
        fig, axs = plt.subplots(figsPerPage, figsize=figsize, squeeze=True)
        figcount = 0
        labels = {"geo": "GEOMETRY", "vib": "VIBRATION", 
                  "occ": "OCCUPATION"}
        offsets = []
        for mode in ["geo", "vib", "occ"]:
            spm = [sp for sp in sps if sp.mode == mode]
            while len(spm) > 0:
                if figcount >= figsPerPage:
                    fig.tight_layout()
                    figs.append(fig)
                    fig, axs = plt.subplots(figsPerPage, figsize=figsize, 
                                            squeeze=True)
                    figcount = 0
                plotpars = spm[:parsPerFig]
                spm = spm[parsPerFig:]
                title = labels[mode]
                if len(rp.disp_blocks) > 1:
                    title += " (search {})".format(searchname[:20])
                axs[figcount].set_title(title)
                pltpoints = [] # x, y, color, size, alpha
                bestpoints = []
                xlabels = []
                for (i, par) in enumerate(plotpars):
                    ind = rp.searchpars.index(par)
                    vals = []
                    for (j, conf) in enumerate(lastconfig):
                        val = (conf[ind]-1) / (par.steps-1)
                        r, g, b, alpha = allcolors[j]
                        if par.linkedTo is None and par.restrictTo is None:
                            alpha = 1.0
                            vals.append(val)
                        else:
                            alpha = 0.5
                        color = (r, g, b, alpha)
                        pltpoints.append((i+1, val, color, 1))
                        if j == 0:
                            bestpoints.append((i+1, val, alpha))
                    if mode != "occ":
                        el = par.el
                    else:
                        el = par.atom.el
                    xlabels.append("#"+str(par.atom.oriN)+"\n"+el)
                    if vals:
                        offsets.append(np.std(vals))
                        # mean = np.mean(vals)
                        # offsets.append(np.mean([abs(v - mean) for v in vals])
                        #                * 2)
                # combine duplicates:
                i = 0
                while i < len(pltpoints):
                    x, y, c, s = pltpoints[i]
                    j = i+1
                    while j < len(pltpoints):
                        if pltpoints[j][0] == x and pltpoints[j][1] == y:
                            s += pltpoints[j][3]
                            if pltpoints[j][2][0] < c[0]:
                                # use the color closer to black
                                c = (pltpoints[j][2][0], c[1], c[2], c[3])
                            pltpoints.pop(j)
                        else:
                            j += 1
                    pltpoints[i] = (x,y,c,s)
                    i += 1
                x, y = [p[0] for p in pltpoints], [p[1] for p in pltpoints]
                c, s = [p[2] for p in pltpoints], [p[3] for p in pltpoints]
                axs[figcount].plot([0, parsPerFig+2], [0.5, 0.5], color='grey', 
                                   alpha=0.2)
                axs[figcount].scatter(x, y, s=s, c=c)
                for (px, py, alpha) in bestpoints:
                    axs[figcount].annotate("", (px+0.05, py), (px+0.25, py), 
                                           arrowprops=dict(arrowstyle="wedge",
                                                           facecolor="black",
                                                           alpha=alpha))
                    axs[figcount].annotate("", (px-0.05, py), (px-0.25, py), 
                                           arrowprops=dict(arrowstyle="wedge",
                                                           facecolor="black",
                                                           alpha=alpha))
                axs[figcount].set_xlim([0, parsPerFig+1])
                axs[figcount].set_ylim([0, 1])
                axs[figcount].set_xticks(list(range(1,parsPerFig+1)))
                axs[figcount].set_xticklabels(xlabels)
                axs[figcount].tick_params(axis='y', which='both', left=False, 
                                right=False, labelleft=False)
                figcount += 1
        if offsets:
            rp.parScatter[-1].append((gens[-1],np.mean(offsets),max(offsets)))
        for i in range(figcount, figsPerPage):
            axs[i].axis('off')
        if not fig in figs:
            fig.tight_layout()
            figs.append(fig)
            
        # save
        if searchname in rp.lastParScatterFigs:
            for f in rp.lastParScatterFigs[searchname]:
                try:
                    plt.close(f)
                except:
                    pass
        rp.lastParScatterFigs[searchname] = figs[1:]
        try:
            pdf = PdfPages(outname)
            for fig in figs:
                pdf.savefig(fig)
        except PermissionError:
            logger.warning("Failed to write to " + outname
                            + ": Permission denied.")
            return 1
        except:
            logger.warning("Failed to write to "+outname)
            return 1
        finally:
            try:
                pdf.close()
            except:
                pass
    finally:
        logger.setLevel(loglevel)
        for fig in [f for f in figs if 
                    not searchname in rp.lastParScatterFigs or
                    not f in rp.lastParScatterFigs[searchname]]:
            try:
                plt.close(fig)
            except:
                pass
    return 0

def writeSearchReportPdf(rp, outname = "Search-report.pdf"):
    global plotting
    if not plotting:
        return 0
    
    allmin = []
    allmax = []
    allmean = []
    allgens = []
    markers = []
    parScatterLines = []  # list of lists [gens, mean, max] per search
    gencount = 0
    for i in range(0, len(rp.searchplots)):
        (name, gens, rmin, rmax, rmean) = rp.searchplots[i]
        markers.append((gencount, "Search "+name))
        allgens.extend([v + gencount for v in gens])
        allmin.extend(rmin)
        allmax.extend(rmax)
        allmean.extend(rmean)
        if rp.parScatter[i]:
            parScatterLines.append(list(zip(*rp.parScatter[i])))
            parScatterLines[-1][0] = [v + gencount for v in 
                                      parScatterLines[-1][0]]
        gencount = allgens[-1]
        
    figsize = (5.8, 8.3)
    figs = []
    # the following will spam the logger with debug messages; disable.
    loglevel = logger.level
    logger.setLevel(logging.INFO)
    try:
        # R-FACTORS AND MEAN SCATTER
        # create figure
        fig, (rfp, msp) = plt.subplots(2, 1, sharex=True, figsize = figsize)
        msp.set_xlabel('Generations')
        rfp.set_ylabel('R-Factor')
        msp.set_ylabel('Parameter scatter')
        #plot markers
        part = 0
        for (i, g) in enumerate(allgens):
            if g > allgens[-1] * 0.2:
                part = len(allgens) - i
                break
        part = max(50, part)
        rfmin, rfmax = min(allmin[-part:]), max(allmean[-part:])
        if rfmax == rfmin:
            rfmin *= 0.95
            rfmax *= 1.05
        rYrange = [rfmin-(rfmax-rfmin)*0.1, rfmax+(rfmax-rfmin)*0.1]
        rYrange[1] = min(rYrange[1], 2*max(allmin) - rYrange[0])
        
        labely = rYrange[0] + (rYrange[1]-rYrange[0])*0.99
        xoff = allgens[-1]*0.005
        for (xpos, label) in markers:
            rfp.axvline(x = xpos, lw = 0.5, c = "black")
            msp.axvline(x = xpos, lw = 0.5, c = "black")
            rfp.text(xpos + xoff, labely, label, rotation=-90, 
                     verticalalignment="top", size=4)
        # plot data
        
        rfp.fill_between(allgens, allmin, allmax, facecolor='grey', 
                         alpha=0.2, label="Range")
        rfp.plot(allgens, allmean, label = "Mean")
        rfp.plot(allgens, allmin, '-', color='black', label = "Best")

        labelled = False
        scattermax = 0
        for (allgens, psmean, psmax) in parScatterLines:
            meanline, = msp.plot(allgens, psmean, '-', color="tab:blue")
            maxline, = msp.plot(allgens, psmax, '-', color="black")
            scattermax = max(scattermax, max(psmax))
            if not labelled:
                meanline.set_label('Mean')
                maxline.set_label('Max')
                labelled = True

        # layout
        rfp.set_ylim(rYrange)
        msp.set_ylim([0, scattermax*1.1])
        msp.ticklabel_format(axis="x", style="sci", scilimits=(0,4))
        fig.tight_layout()
        rfp.legend(loc="lower left")
        msp.legend(loc="lower left")
        
        # ADD SCATTERS
        for k in rp.lastParScatterFigs.keys():
            figs.append(rp.lastParScatterFigs[k])
        # save
        try:
            pdf = PdfPages(outname)
            pdf.savefig(fig)
            for k in rp.lastParScatterFigs.keys():
                for f in rp.lastParScatterFigs[k]:
                    pdf.savefig(f)
        except PermissionError:
            logger.warning("Failed to write to " + outname
                            + ": Permission denied.")
            return 1
        except:
            logger.warning("Failed to write to "+outname)
            return 1
        finally:
            try:
                pdf.close()
            except:
                pass
    finally:
        logger.setLevel(loglevel)
        try:
            plt.close(fig)
        except:
            pass
    return 0

# def closePdfReportFigs(rp):           # !!! CLEANUP
#     global plotting
#     if not plotting:
#         return 0
    
#     for searchname in rp.lastParScatterFigs:
#         for f in rp.lastParScatterFigs[searchname]:
#             try:
#                 plt.close(f)
#             except:
#                 pass