# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 13:35:50 2020

@author: Florian Kraushofer

Functions for reading, processing and writing files relevant to the search
"""

import logging
import numpy as np
import fortranformat as ff
import copy
import os
import random
import shutil

from tleedmlib.files.beams import writeAUXEXPBEAMS
from tleedmlib.files.poscar import writeCONTCAR
from tleedmlib.files.vibrocc import writeVIBROCC
from tleedmlib.base import BackwardsReader, readIntLine
from tleedmlib.leedbase import getBeamCorrespondence

logger = logging.getLogger("tleedm.files.iosearch")


def readSDTL_next(filename="SD.TL", offset=0):
    """
    Reads SDTL from offset to end, returns new offset and the content in
    between as string.

    Parameters
    ----------
    filename : str, optional
        The file to read from. The default is "SD.TL".
    offset : int, optional
        Where in the file to start. The default is 0.

    Returns
    -------
    offset : int
        The position to which SD.TL was read; this should be passed as offset
        when this function is called next.
    content : str
        The content of the file found between the old offset and the new.

    """
    try:
        with open(filename, "r") as rf:
            rf.seek(offset)
            content = rf.read()
            newoffset = rf.tell()
        return(newoffset, content)
    except Exception:
        logger.error("Error reading SD.TL file")
        return (offset, "")     # return old offset, no content


def readSDTL_blocks(content, whichR=0, logInfo=False):
    """
    Attempts to interpret a given string as one or more blocks of an SD.TL
    file.

    Parameters
    ----------
    content : str
        A block of data as read from SD.TL.
    whichR : int, optional
        Which r-factor values to use (average / integer / fractional)
    logInfo : bool, optional
        Whether some basic information should be printed to logger.info

    Returns
    -------
    returnList : list
        List consists of one entry per data block. Each entry of the list is
        a tuple (gen, rfacs, configs), where gen is the generation index,
        rfacs an ordered list of R-factors, and configs the corresponding
        parameter configurations. onfig is a tuple listing (percent, dc)
        for each domain, with dc the parameter values of that domain.

    """
    returnList = []
    blocklist = content.split("CCCCCCCCCCC    GENERATION")[1:]
    for block in blocklist:
        gen = 0
        try:
            gen = int(block[:13])
            if logInfo:
                logger.info("Reading search results from SD.TL, "
                            "generation {}".format(gen))
        except ValueError:
            logger.warning("While reading SD.TL, could not interpret "
                           "generation number "+block[:13])
        lines = block.split("\n")
        rfacs = []
        configs = []
        dpars = []  # list of (list of (percent, parameters)) by domain
        for line in lines:
            if "|" in line and line[:3] != "IND":
                # this is a line with data in it
                if line.split("|")[0].strip():
                    # first line - contains R-factor
                    if dpars:
                        configs.append(tuple(dpars))
                        dpars = []
                    try:
                        rav = float(line.split("|")[2 + whichR])  # average R
                        rfacs.append(rav)
                    except ValueError:
                        logger.error("Could not read R-factor in SD.TL line:\n"
                                     + line)
                try:
                    percent = int(line.split("|")[-2].strip()[:-1])
                    valstring = line.split("|")[-1].rstrip()
                    pars = readIntLine(valstring, width=4)
                    dpars.append((percent, pars))
                except ValueError:
                    logger.error("Could not read values in SD.TL line:\n"
                                 + line)
        if dpars:
            configs.append(tuple(dpars))
        if not all([len(dp) == len(configs[0]) for dp in configs]):
            logger.warning("A line in SD.TL contains fewer values than "
                           "the others. Skipping SD.TL block.")
            continue
        if gen != 0 and len(rfacs) > 0 and len(configs) > 0:
            returnList.append((gen, rfacs, tuple(configs)))
        else:
            logger.warning("A block in SD.TL was read but not understood.")
    return returnList


def readSDTL_end(filename="SD.TL"):
    """
    Reads the last generation block from the SD.TL file, starting from the
    last line containing a GENERATION label.

    Parameters
    ----------
    filename : str, optional
        Which file to read

    Returns
    -------
    lines : str
        The contents of the file from the last line containing GENERATION to
        the end.

    """
    # get the last block from SD.TL:
    bwr = BackwardsReader(filename)
    lines = [""]
    while "CCCCCCCCCCC    GENERATION" not in lines[-1] and len(bwr.data) > 0:
        lines.append(bwr.readline().rstrip())
    lines.reverse()
    bwr.close()
    return lines


def readDataChem(rp, source):
    """
    Reads the data from a list of data.chem files, or a single file.
    Returns

    Parameters
    ----------
    rp : Rparams
        Run parameters.
    source : either a filename, or an iterable containing multiple file names.
        The files from which to read data; usually of format 'data*.chem'.

    Returns
    -------
    list
        Tuples (rfac, config), where config is a tuple listing (percent, dc)
        for each domain, with dc the parameter values of that domain.

    """
    def getPercent(domain_step, dp_vals):
        """Transforms from integer domain area values to percent"""
        dsteps = int(100 / domain_step) + 1
        shifted = [v-1 for v in dp_vals]
        norm = 1
        if sum(shifted) > 0:
            norm = (dsteps-1)/sum(shifted)
        # note that TensErLEED uses implicit floor int conversion here up to
        #  version 1.6, so results may be inconsistent with SD.TL.
        pl = [int(round(v*norm))*domain_step for v in shifted[:-1]]
        pl.append(100 - sum(pl))
        return pl

    if type(source) == str:
        source = [source]
    lines = set()
    for f in source:
        with open(f, "r") as rf:
            lines.update(rf.readlines()[1:])
    returnList = []
    parslen = 0
    for line in lines:
        try:
            rfac = float(line.split("|")[0].strip())
            valstring = line.split("|")[1].rstrip()
            pars = readIntLine(valstring, width=4)
        except ValueError:
            logger.warning("Could not read values in data.chem line:\n"+line)
            continue
        if len(pars) != parslen:
            if parslen == 0:
                parslen = len(pars)
            else:
                continue     # incomplete line, read while file was written
        if not rp.domainParams:
            returnList.append((rfac, ((100, pars[:-1]),)))
        else:
            dp_vals = pars[-len(rp.domainParams):]
            percent = getPercent(rp.DOMAIN_STEP, dp_vals)
            dpars = []
            for v in dp_vals:
                dpars.append(tuple(pars[:v]))
                pars = pars[v:]
            returnList.append((rfac, tuple(zip(percent, dpars))))
    return returnList


def writeRfInfo(sl, rp, filename="rf.info"):
    """
    Generates r-factor parameters for the search, combines them with the
    experimental beams in AUXEXPBEAMS format to make the entire input for the
    search, returns that as a string.

    Parameters
    ----------
    sl : Slab
        The Slab object used in the search. Only used for determining
        equivalent beams.
    rp : Rparams
        Run parameters.
    filename : str, optional
        Name of the output file. The default is "rf.info".

    Returns
    -------
    output : str
        Content of the output file.

    """
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
    step = min(expEnergies[1]-expEnergies[0], rp.THEO_ENERGIES[2])
    if rp.IV_SHIFT_RANGE[2] > 0:
        vincr = rp.IV_SHIFT_RANGE[2]
        step = min(step, vincr)
    else:
        vincr = step
    # find correspondence experimental to theoretical beams:
    beamcorr = getBeamCorrespondence(sl, rp)
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
    output += i3x25.write([n + 1 for n in range(0, len(rp.expbeams))]) + "\n"
    output += f31x25.write([1]*len(rp.expbeams))+"\n"
    auxexpbeams = writeAUXEXPBEAMS(rp.expbeams, header=rp.systemName,
                                   write=True, numbers=True)
    output += auxexpbeams

    try:
        with open(filename, 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write "+filename)
        raise
    logger.debug("Wrote to "+filename+" successfully")
    return output


def generateSearchInput(sl, rp, steuOnly=False, cull=False, info=True):
    """
    Generates a PARAM and a search.steu file for the search.

    Parameters
    ----------
    sl : Slab
        The Slab object used in the search. Only used if there are no domains.
    rp : Rparams
        Run parameters.
    steuOnly : bool, optional
        If True, PARAM and restrict.f will not be written. The default is
        False.
    cull : bool, optional
        If True, part of the population (given by rp.SEARCH_CULL) will be
        removed and replaced according to rp.SEARCH_CULL. The default is False.
    info : bool, optional
        If False, no info messages will be printed.

    Returns
    -------
    None

    """
    # first generate list of SearchPar objects and figure out search parameters
    rp.generateSearchPars(sl)

    # if search population is undefined, calculate a default:
    if rp.SEARCH_POPULATION == 0:
        spop = min(48, 15 + rp.indyPars)
        rp.SEARCH_POPULATION = int(np.ceil(spop / rp.N_CORES)) * rp.N_CORES

    # calculate some more things for later
    expEnergies = []
    totalrange = 0
    for b in rp.expbeams:
        expEnergies.extend([k for k in b.intens if k not in expEnergies])
        totalrange += (min(max(b.intens), rp.THEO_ENERGIES[1])
                       - max(min(b.intens), rp.THEO_ENERGIES[0]))
    if info:
        logger.info("Total energy range from experimental beams is "
                    "{:.2g} eV ({} independent fit parameters, {:.2g} eV per "
                    "parameter)"
                    .format(totalrange, rp.indyPars, totalrange / rp.indyPars))
    theoEnergies = int((rp.THEO_ENERGIES[1]-rp.THEO_ENERGIES[0])
                       / rp.THEO_ENERGIES[2]) + 1
    if theoEnergies >= len(expEnergies):
        logger.warning("Theoretical beams have more data points than "
                       "experimental beams")
        rp.setHaltingLevel(1)
        # TODO: will this crash? If so, raise here

    # merge offsets with displacement lists
    if rp.domainParams:
        attodo = [at for dp in rp.domainParams for at in dp.rp.search_atlist]
        ndom = len(rp.domainParams)
        astep = rp.DOMAIN_STEP
        nplaces = max([len(dp.rp.search_atlist) for dp in rp.domainParams])
    else:
        attodo = rp.search_atlist
        ndom = 1
        astep = 100
        nplaces = len(rp.search_atlist)
    for at in attodo:
        if at.oriState is None:
            for el in at.offset_occ:
                if el not in at.disp_occ:
                    logger.error(
                        "{} has occupation offset defined for element {}, "
                        "which does not appear in atom displacement list. "
                        "Value will be skipped, this may cause errors."
                        .format(at, el))
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
            at.disp_geo_offset = {"all": [np.array([0., 0., 0.])]}

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
      parameter (MNDOM = {})""".format(ndom)
    output += """
C MNPLACES IS NUMBER OF DIFFERENT ATOMIC SITES IN CURRENT VARIATION
      PARAMETER(MNPLACES = {})""".format(nplaces)
    output += """
C MNFILES IS MAXIMUM NUMBER OF TLEED-FILES PER SITE
      PARAMETER(MNFILES = {})""".format(rp.search_maxfiles)
    output += """
C MNCONCS IS MAX NUMBER OF CONCENTRATION STEPS PER SITE
      PARAMETER(MNCONCS = {})\n""".format(rp.search_maxconc)
    output += ("C MNPRMK IS NUMBER OF PARAMETERS - INCL ONE CONC PARAMETER "
               "PER SITE - incl. 1 per domain!\n")
    output += "      PARAMETER(MNPRMK = {})".format(len(rp.searchpars))
    output += """
C MNCSTEP IS MAX NUMBER OF VARIATIONS (geo. times therm.) IN 1 FILE
C MNATOMS IS RELICT FROM OLDER VERSIONS
      PARAMETER (MNCSTEP = {}, MNATOMS = 1)""".format(rp.mncstep)
    # write PARAM
    if not steuOnly:
        try:
            with open("PARAM", "w") as wf:
                wf.write(output)
        except Exception:
            logger.error("Failed at writing PARAM file for search.")
            raise

    # now search.steu
    output = ""
    # while producing search.steu, collect information on
    #   parameters for searchpars.info
    info = ("#".rjust(4) + "label".rjust(7) + "atom".rjust(7) + "elem".rjust(7)
            + "mode".rjust(7) + "steps".rjust(7) + "constr".rjust(7) + "\n")

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
    output += i3.write([astep]).ljust(16) + "area fraction step width (%)\n"
    output += ("SD.TL           name of search document file "
               "(max. 10 characters)\n")
    output += i3.write([ndom]).ljust(16) + ("Number of domains under "
                                            "consideration\n")
    uniquenames = []
    for k in range(0, ndom):
        if ndom == 1:
            crp = rp
            csl = sl
            frompath = ""
        else:
            crp = rp.domainParams[k].rp
            csl = rp.domainParams[k].sl
            frompath = rp.domainParams[k].workdir
            info += "---- DOMAIN {} ----\n".format(k+1)
        prev_parcount = parcount
        output += ("======= Information about Domain {}: ====================="
                   "=======================\n").format(k+1)
        output += (
            i3.write([len(crp.search_atlist)]).ljust(16) + "Number of atomic "
            "sites in variation: Domain {}\n".format(k+1))
        displistcount = len(csl.displists)+1
        surfats = csl.getSurfaceAtoms(crp)
        for (i, at) in enumerate(crp.search_atlist):
            printvac = False
            output += (
                "------- Information about site {}: -----------------------"
                "-----------------------\n".format(i+1))
            surf = 1 if at in surfats else 0
            output += (i3.write([surf]).ljust(16) + "Surface (0/1)\n")
            if at.displist in csl.displists:
                dlind = csl.displists.index(at.displist) + 1
            else:
                dlind = displistcount
                displistcount += 1
            output += (i3.write([dlind]).ljust(16) + "Atom number\n")
            output += (
                i3.write([len(at.deltasGenerated)]).ljust(16) + "No. of "
                "different files for Atom no. {}\n".format(i+1))
            for (j, deltafile) in enumerate(at.deltasGenerated):
                name = deltafile
                if frompath:  # need to get the file
                    name = "D{}_".format(k+1) + deltafile
                    if len(name) > 15:
                        un = 1
                        while "D{}_DEL_{}".format(k+1, un) in uniquenames:
                            un += 1
                        name = "D{}_DEL_{}".format(k+1, un)
                        uniquenames.append(name)
                    try:
                        shutil.copy2(os.path.join(frompath, deltafile), name)
                    except Exception:
                        logger.error("Error getting Delta file {} for search"
                                     .format(os.path.relpath(
                                         os.path.join(frompath, deltafile))))
                        raise
                output += "****Information about file {}:\n".format(j+1)
                output += (name.ljust(16) + "Name of file {} (max. 15 "
                           "characters)\n".format(j+1))
                output += "  1             Formatted(0/1)\n"
                el = deltafile.split("_")[-2]
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
                output += (i3.write([types]).ljust(16) + "Types of parameters "
                           "in file {}\n".format(j+1))
                label = i1.write([i+1])+i1.write([j+1])
                constr = {"vib": "-", "geo": "-"}
                for mode in ["vib", "geo"]:
                    spl = [s for s in crp.searchpars if s.el == el
                           and s.atom == at and s.mode == mode]
                    if spl and spl[0].restrictTo is not None:
                        sp = spl[0]
                        if type(sp.restrictTo) == int:
                            constr[mode] = str(sp.restrictTo)
                        else:
                            constr[mode] = "#"+str(crp.searchpars.index(
                                                                sp.restrictTo)
                                                   + prev_parcount + 1)
                    elif spl and spl[0].linkedTo is not None:
                        constr[mode] = "#"+str(
                            crp.searchpars.index(spl[0].linkedTo)
                            + prev_parcount + 1)
                if vib > 0:
                    output += (i3.write([vib]).ljust(16) +
                               "vibrational steps\n")
                    parcount += 1
                    info += (str(parcount).rjust(4) + ("P"+label).rjust(7)
                             + str(at.oriN).rjust(7) + el.rjust(7)
                             + "vib".rjust(7) + str(vib).rjust(7)
                             + constr["vib"].rjust(7) + "\n")
                    nsteps.append(vib)
                if geo > 0:
                    output += (i3.write([geo]).ljust(16) +
                               "geometrical steps\n")
                    parcount += 1
                    info += (str(parcount).rjust(4) + ("P"+label).rjust(7)
                             + str(at.oriN).rjust(7) + el.rjust(7)
                             + "geo".rjust(7) + str(geo).rjust(7)
                             + constr["geo"].rjust(7) + "\n")
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
            spl = [s for s in crp.searchpars if s.atom == at and
                   s.mode == "occ"]
            if spl and spl[0].restrictTo is not None:
                sp = spl[0]
                if type(sp.restrictTo) == int:
                    constr = str(sp.restrictTo)
                else:
                    constr = "#"+str(crp.searchpars.index(sp.restrictTo)+1)
            elif spl and spl[0].linkedTo is not None:
                constr = "#"+str(crp.searchpars.index(spl[0].linkedTo)+1)
            info += (str(parcount).rjust(4) + ("C"+label).rjust(7)
                     + str(at.oriN).rjust(7) + "-".rjust(7) + "occ".rjust(7)
                     + str(occsteps).rjust(7) + constr.rjust(7) + "\n")
            nsteps.append(occsteps)
    # add info for domain step parameters
    for k in range(0, ndom):
        parcount += 1
        if ndom == 1:
            astepnum = 1
        else:
            astepnum = int(100 / astep) + 1
        info += (str(parcount).rjust(4) + ("-".rjust(7)*3) + "dom".rjust(7)
                 + str(astepnum).rjust(7) + "-".rjust(7).rjust(7) + "\n")
        nsteps.append(astepnum)
    # now starting configuration
    output += ("Information about start configuration: (Parameters for each "
               "place, conc for each place)\n")
    if rp.SEARCH_START == "control":
        controlpath = ""
        if os.path.isfile("control.chem"):
            controlpath = "control.chem"
        elif os.path.isfile(os.path.join("AUX", "control.chem")):
            controlpath = os.path.join("AUX", "control.chem")
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
        except Exception:
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
                        logger.warning(
                            "SEARCH_CULL parameter too large: would cull "
                            "entire population. Culling will be skipped.")
                        ncull = 0
                if any([sp.parabolaFit["min"] is not None
                        for sp in rp.searchpars]):
                    # replace one by predicted best
                    getPredicted = True
                else:
                    getPredicted = False
                nsurvive = rp.SEARCH_POPULATION - ncull
                clines = controllines[2:]
                csurvive = []
                if (rp.SEARCH_CULL_TYPE == "genetic" or
                        getPredicted):  # prepare readable clines
                    try:
                        csurvive = [readIntLine(s, width=3)
                                    for s in clines[:nsurvive]]
                    except ValueError:
                        if rp.SEARCH_CULL_TYPE == "genetic":
                            logger.warning(
                                "SEARCH_CULL: Failed to read old "
                                "configuration from control.chem, cannot run "
                                "genetic algorithm. Defaulting to cloning.")
                            rp.setHaltingLevel(1)
                            csurvive = []
                for (i, line) in enumerate(clines):
                    if i < nsurvive:
                        output += line
                    elif (rp.SEARCH_CULL_TYPE == "random" or
                          (rp.SEARCH_CULL_TYPE == "genetic" and csurvive)
                          or getPredicted):
                        if getPredicted:
                            bc = None
                            if csurvive:
                                bc = csurvive[0]
                            nc = rp.getPredictConfig(
                                best_config=bc,
                                mincurv=rp.PARABOLA_FIT["mincurv"])
                            getPredicted = False
                        elif rp.SEARCH_CULL_TYPE == "random":
                            nc = rp.getRandomConfig()
                        else:  # "genetic"
                            nc = rp.getOffspringConfig(csurvive)
                        ol = ""
                        for v in nc:
                            ol += i3.write([v])
                        output += ol + "\n"
                    else:    # "cloning"
                        output += clines[random.randrange(0, nsurvive)]
            else:  # don't cull
                for (i, line) in enumerate(controllines):
                    if i > 1:  # first 2 lines are empty
                        output += line
        elif rp.SEARCH_START == "centered":
            for i in range(0, rp.SEARCH_POPULATION):
                ol = ""
                for n in nsteps:
                    ol += i3.write([(n+1)/2])
                output += ol + "\n"
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
    except Exception:
        logger.error("Failed to write search.steu file for search.")
        raise
    # write searchpars.info
    try:
        with open("searchpars.info", "w") as wf:
            wf.write(info)
    except Exception:
        logger.error(
            "Failed to write searchpars.info file. This file is "
            "only for user information and will not affect operation.")
    if steuOnly:
        logger.debug("Wrote search.steu input for search.")
        return None
    # now restrict.f
    output = ("""
C  This subroutine serves to restrict individual parameters to certain values
C  inside the search. These values are stored in the PARIND(IPARAM,IPOP) array.
C  Parameters are counted as listed in control file search.steu, including
C  a concentration parameter after each atomic site, and the domain weight
C  parameters for each domain at the end of the PARIND array.
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
            output += "      PARIND({},IPOP) = PARIND({},IPOP)\n".format(
                i + 1, rp.searchpars.index(sp.restrictTo) + 1)
    output += ("""
C  end restrictions

      RETURN

      END
""")
    # write restrict.f
    try:
        with open("restrict.f", "w") as wf:
            wf.write(output)
    except Exception:
        logger.error("Failed to write restrict.f file for search.")
        raise
    logger.debug("Wrote input files for search: PARAM, search.steu, "
                 "restrict.f")
    return None


def writeSearchOutput(sl, rp, parinds=None, silent=False, suffix=""):
    """
    Modifies data in sl and rp to reflect the search result given by
    parinds, then writes POSCAR_OUT and VIBROCC_OUT.

    Parameters
    ----------
    sl : Slab
        The Slab object to be modified.
    rp : Rparams
        The run parameters object to be modified.
    parinds : iterable, optional
        The final parameter values to be applied. If None,
        rp.searchResultConfig will be used instead. Can be float, in which
        case linear interpolation between the two adjacent displacements will
        be used.
    silent : bool, optional
        Suppresses output to log. The default is False.
    suffix : str, optional
        String to be appended to the POSCAR_OUT and VIBROCC_OUT file names.

    Returns
    -------
    None

    """
    def interpolate(var, ind):
        if type(ind) == int:
            return var[ind]
        o = ind % 1
        return (o * var[int(np.ceil(ind))]
                + (1 - o) * var[int(np.floor(ind))])

    if parinds is None:
        if rp.searchResultConfig is None:
            logger.error("Failed to write search output: No configuration "
                         "passed.")
            raise RuntimeError("Failed to write search output")
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
        if len(sps) > 0:
            # first deal with occupation:
            par = [sp for sp in sps if sp.mode == "occ"][0]  # exactly one
            i = rp.searchpars.index(par)
            newocc = {}
            for el in at.disp_occ:
                # store as offset for now, recalculate for sites later
                at.offset_occ[el] = (interpolate(at.disp_occ[el], parinds[i]-1)
                                     - at.site.oriState.occ[el])
                newocc[el] = interpolate(at.disp_occ[el], parinds[i]-1)
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
                at.offset_vib[el] = interpolate(at.disp_vib[dict_el],
                                                parinds[i]-1)
            # now position - recalculate right away:
            geopars = [sp for sp in sps if sp.mode == "geo"]
            # can be empty, or one per element
            if not geopars:
                continue
            totalocc = 0
            newpos = {}     # new position per element
            totalpos = np.array([0., 0., 0.])
            for par in geopars:
                i = rp.searchpars.index(par)
                el = par.el
                if el not in at.disp_geo:
                    dict_el = "all"
                else:
                    dict_el = el
                disp = np.copy(interpolate(at.disp_geo[dict_el], parinds[i]-1))
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
    fn = "POSCAR_OUT" + suffix + "_" + rp.timestamp
    tmpslab = copy.deepcopy(sl)
    tmpslab.sortOriginal()
    try:
        writeCONTCAR(tmpslab, filename=fn, comments="all", silent=silent)
    except Exception:
        logger.error("Exception occured while writing POSCAR_OUT" + suffix,
                     exc_info=rp.LOG_DEBUG)
        rp.setHaltingLevel(2)
    if not np.isclose(rp.SYMMETRY_CELL_TRANSFORM, np.identity(2)).all():
        tmpslab = sl.makeSymBaseSlab(rp)
        fn = "POSCAR_OUT_mincell" + suffix + "_" + rp.timestamp
        try:
            writeCONTCAR(tmpslab, filename=fn, silent=silent)
        except Exception:
            logger.warning(
                "Exception occured while writing POSCAR_OUT_mincell" + suffix,
                exc_info=rp.LOG_DEBUG)
    fn = "VIBROCC_OUT" + suffix + "_" + rp.timestamp
    try:
        writeVIBROCC(sl, rp, filename=fn, silent=silent)
    except Exception:
        logger.error("Exception occured while writing VIBROCC_OUT" + suffix,
                     exc_info=rp.LOG_DEBUG)
        rp.setHaltingLevel(2)
    return
