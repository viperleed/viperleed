"""Module iosearch of viperleed.calc.files.

Functions for reading, processing and writing files relevant
to the search.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import copy
import logging
import os
from pathlib import Path
import random
import shutil
import time

import fortranformat as ff
import numpy as np

from viperleed.calc.files import poscar
from viperleed.calc.files.beams import writeAUXEXPBEAMS
from viperleed.calc.files.iorfactor import largest_nr_grid_points
from viperleed.calc.files.iorfactor import prepare_rfactor_energy_ranges
from viperleed.calc.files.vibrocc import writeVIBROCC
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.base import BackwardsReader, readIntLine

logger = logging.getLogger(__name__)


class SearchIORaceConditionError(Exception):
    """Raised if reading of control.chem does not return the expected number
    of lines"""
    def __init__(self, message):
        super().__init__(message)


class SearchIOEmptyFileError(Exception):
    """Raised if file read for the search has no content"""
    def __init__(self, message):
        super().__init__(message)


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


def readSDTL_blocks(content, whichR=0, print_info=False, n_expect=0):
    """
    Attempts to interpret a given string as one or more blocks of an SD.TL
    file.

    Parameters
    ----------
    content : str
        A block of data as read from SD.TL.
    whichR : int, optional
        Which r-factor values to use (average / integer / fractional)
    print_info : bool, optional
        Whether some basic information should be printed to logger.info and
        logger.warning. The default is False.
    n_expect : int, optional
        Number of configurations expected per block. If a block contains
        fewer, it will be ignored.

    Returns
    -------
    returnList : list
        List consists of one entry per data block. Each entry of the list is
        a tuple (gen, rfacs, configs), where gen is the generation index,
        rfacs an ordered list of R-factors, and configs the corresponding
        parameter configurations. config is a tuple listing (percent, dc)
        for each domain, with dc the parameter values of that domain.

    """
    returnList = []
    blocklist = content.split("CCCCCCCCCCC    GENERATION")[1:]
    for block in blocklist:
        gen = 0
        try:
            gen = int(block[:13])
            if print_info:
                logger.info("Reading search results from SD.TL, "
                            "generation {}".format(gen))
        except ValueError:
            if print_info:
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
                        if print_info:
                            logger.error("Could not read R-factor in SD.TL line:\n"
                                        + line)
                try:
                    percent = int(line.split("|")[-2].strip()[:-1])
                    valstring = line.split("|")[-1].rstrip()
                    pars = readIntLine(valstring, width=4)
                    dpars.append((percent, pars))
                except ValueError:
                    if print_info:
                        logger.error("Could not read values in SD.TL line:\n"
                                    + line)
        if dpars:
            configs.append(tuple(dpars))
        if not all([len(dp) == len(configs[0]) for dp in configs]):
            if print_info:
                logger.warning("A line in SD.TL contains fewer values than "
                            "the others. Skipping SD.TL block.")
            continue
        if gen != 0 and len(rfacs) > 0 and len(configs) > 0:
            returnList.append((gen, rfacs, tuple(configs)))
        elif len(configs) < n_expect:
            if print_info:
                logger.warning("A block in SD.TL contains fewer configurations "
                            "than expected.")
        else:
            if print_info:
                logger.warning("A block in SD.TL was read but not understood.")
    return returnList


def repeat_fetch_SDTL_last_block(which_beams,
                                 expected_params,
                                 final,
                                 max_repeats=2000,
                                 wait_time=5):
    print_info = final
    content = None
    for repeat in range(max_repeats):
        content = _fetch_SDTL_last_block(which_beams,
                                         expected_params,
                                         print_info)
        try:  # try again if content is not complete
            n_search_params_found = len(content[0][2])
        except IndexError:
            continue
        if n_search_params_found == expected_params:
            logger.debug("Read complete block from SD.TL file after "
                         f"{repeat+1} read attempt(s).")
            return content
        elif n_search_params_found >= expected_params:
            raise RuntimeError(f"Expected a maximum of {expected_params} "
                               f"parameters in SD.TL, but found "
                               "{n_search_params_found}.")
        # wait and try again
        time.sleep(wait_time/1000)
        print_info = False
    if not content:
        raise SearchIOEmptyFileError("No data found in SD.TL file.")
    # if we get here, we have exceeded the maximum number of repeats
    raise SearchIORaceConditionError(f"Could not read complete block from "
                                     "SD.TL file.")

def _fetch_SDTL_last_block(which_beams, n_expect, print_info=False):
    try:
        lines = readSDTL_end(n_expect=n_expect)
    except FileNotFoundError:
        logger.error("Could not process Search results: SD.TL file not "
                     "found.")
        raise
    except Exception:
        logger.error("Failed to get last block from SD.TL file.")
        raise
    # check if the last block contains the expected number of lines
    lines_str = "\n".join(lines)
    sdtl_content = readSDTL_blocks("\n".join(lines),
                                   whichR=which_beams,
                                   print_info=print_info,
                                   n_expect=n_expect)
    return sdtl_content


def readSDTL_end(filename="SD.TL", n_expect=0):
    """
    Reads the last generation block from the SD.TL file, starting from the
    last line containing a GENERATION label.

    Parameters
    ----------
    filename : str, optional
        Which file to read
    n_expect : int, optional
        Number of configurations expected per block. If the last block contains
        fewer, it will be ignored, instead reading the second-to-last. Set to
        0 to read the last block irrespective of length.

    Returns
    -------
    lines : str
        The contents of the file from the last line containing GENERATION to
        the end.

    """
    # get the last block from SD.TL:
    bwr = BackwardsReader(filename)
    try:
        lines = [""]
        while len(bwr.data) > 0:
            while ("CCCCCCCCCCC    GENERATION" not in lines[-1]
                   and len(bwr.data) > 0):
                lines.append(bwr.readline().rstrip())
            if len([line for line in lines if "|" in line]) > n_expect:
                break
            lines = [""]
        lines.reverse()
    finally:
        bwr.close()
    return lines


def readDataChem(rp, source, cutoff=0, max_configs=0):
    """
    Reads the data from a list of data.chem files, or a single file.

    Parameters
    ----------
    rp : Rparams
        Run parameters.
    source : either a filename, or an iterable containing multiple file names.
        The files from which to read data; usually of format 'data*.chem'.
    cutoff : float
        0 to read all, else specifies the maximum R-factor to be stored.
        All structures with higher R will be discarded.
    max_configs : int
        0 to read all, else specifies the maximum number of configurations to
        read in. Will sort configurations and read in the ones with lowest
        R-factors first, discard the rest.

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
    if max_configs == 0:
        max_configs = len(lines)
    if len(lines) > max_configs:
        lines = set(sorted(list(lines))[:max_configs])
    for line in lines:
        try:
            rfac = float(line.split("|")[0].strip())
            if cutoff != 0 and rfac > cutoff:
                continue
            valstring = line.split("|")[1].rstrip()
            pars = readIntLine(valstring, width=4)
        except (ValueError, IndexError):
            logger.debug("Could not read values in data.chem line:\n"+line)
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
            for v in [len(dp.rp.searchpars) for dp in rp.domainParams]:
                dpars.append(tuple(pars[:v]))
                pars = pars[v:]
            returnList.append((rfac, tuple(zip(percent, dpars))))
    return returnList


def writeRfInfo(sl, rp, file_path="rf.info"):
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
    file_path : pathlike or str
        Pathlike to or name of the output file. If str uses
        rp.workdir /filename. The default is "rf.info".

    Returns
    -------
    output : str
        Content of the output file.
    """
    _, theo_range, _, vincr = prepare_rfactor_energy_ranges(rp)

    # find correspondence experimental to theoretical beams:
    beamcorr = leedbase.getBeamCorrespondence(sl, rp)
    # integer & fractional beams
    iorf = []
    for beam in rp.expbeams:
        if beam.hkfrac[0].denominator != 1 or beam.hkfrac[1].denominator != 1:
            iorf.append(2)
        else:
            iorf.append(1)

    if rp.TL_VERSION < 1.7:
        formatter = {'energies': ff.FortranRecordWriter('F7.2'),
                     'int': ff.FortranRecordWriter('I3'),
                     'beams': ff.FortranRecordWriter("25I3"),
                     'weights': ff.FortranRecordWriter("25F3.1")
                     }
    else:
        formatter = {'energies': ff.FortranRecordWriter('F7.2'),
                     'int': ff.FortranRecordWriter('I4'),
                     'beams': ff.FortranRecordWriter("25I4"),
                     'weights': ff.FortranRecordWriter("25F4.1")
                     }
    output = (formatter['energies'].write([theo_range.min]).ljust(16)
              + "EMIN\n")
    output += (
        formatter['energies'].write([theo_range.max + 0.1 * vincr]).ljust(16)
        + "EMAX\n"
        )
    # !!! BULLSHIT RESULTS WHEN EMAX > MAX ENERGY IN DELTA FILES                # TODO: Issue #138
    output += (formatter['energies'].write([vincr]).ljust(16) + "EINCR\n")
    output += (formatter['int'].write([0]).ljust(16)
               + "IPR - determines amount of output to stdout\n")
    output += (formatter['energies'].write([rp.V0_IMAG]).ljust(16) + "VI\n")
    output += (formatter['energies'].write([0.]).ljust(16) + "V0RR\n")
    output += (formatter['energies'].write([rp.IV_SHIFT_RANGE.start]).ljust(16)
               + "V01\n")
    output += (formatter['energies'].write([rp.IV_SHIFT_RANGE.stop]).ljust(16)
               + "V02\n")
    output += (formatter['energies'].write([vincr]).ljust(16) + "VINCR\n")
    output += (formatter['int'].write([rp.R_FACTOR_SMOOTH]).ljust(16)
               + "ISMOTH\n")
    output += (formatter['int'].write([0]).ljust(16)
               + "EOT - 0: exp format, 1: van Hove format\n")
    output += ((formatter['int'].write([len(rp.ivbeams)])
                + formatter['int'].write([len(rp.expbeams)])).ljust(16)
               + "NTH NEX\n")
    # numbers of theoretical beams, as they correspond to experimental beams
    output += (formatter['beams'].write([n+1 for n in beamcorr]) + "\n")
    output += " DATA MITTEL (integer & fractional beams) :\n"
    output += formatter['beams'].write(iorf) + "\n"
    output += " exp - th relationship IBP, beam weights WB\n"
    output += formatter['beams'].write([n + 1 for n
                                        in range(0, len(rp.expbeams))]) + "\n"
    output += formatter['weights'].write([1]*len(rp.expbeams))+"\n"
    auxexpbeams = writeAUXEXPBEAMS(rp.expbeams, header=rp.systemName,
                                   write=True, numbers=True,
                                   version=rp.TL_VERSION)
    output += auxexpbeams

    if isinstance(file_path, str):
        _file_path = rp.workdir / file_path
    else:
        _file_path = file_path
    try:
        with open(_file_path, 'w') as wf:
            wf.write(output)
    except Exception:
        logger.error(f"Failed to write {_file_path}")
        raise
    logger.debug(f"Wrote to {_file_path} successfully")
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

    if rp.indyPars == 0:
        logger.warning("The current search includes no variations.")
        rp.indyPars = 1

    # if search population is undefined, calculate a default:
    if rp.SEARCH_POPULATION == 0:
        spop = min(48, 15 + rp.indyPars)
        rp.SEARCH_POPULATION = int(np.ceil(spop / rp.N_CORES)) * rp.N_CORES

    # calculate some more things for later
    expEnergies = set([k for b in rp.expbeams for k in b.intens])
    if info:
        totalrange = rp.total_energy_range()
        logger.info("Total energy range from experimental beams is "
                    "{:.2g} eV ({} independent fit parameters, {:.2g} eV per "
                    "parameter)"
                    .format(totalrange, rp.indyPars, totalrange / rp.indyPars))
    n_theo_energies = rp.THEO_ENERGIES.n_energies
    if n_theo_energies >= len(expEnergies):
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
                    del at.disp_occ[el] # TODO: why? we write and then delete? Should be offset_occ? If yes, loop should go over .copy
            for (disp, offset) in [(at.disp_geo, at.offset_geo),
                           (at.disp_vib, at.offset_vib),
                           (at.disp_geo, at.disp_geo_offset)]:
                dl = []
                for el in offset:
                    if el not in disp:
                        disp[el] = copy.copy(disp["all"])
                    disp[el] = [v + offset[el] for v in disp[el]]
                    dl.append(el)
                for el in dl:
                    if offset != at.disp_geo_offset:
                        del offset[el]
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
    max_nr_data_points = largest_nr_grid_points(rp, rp.theobeams['refcalc'],
                                                False)
    output += """
C MNDATA IS MAX. NUMBER OF DATA POINTS IN EXPERIMENTAL BEAMS
      PARAMETER(MNDATA = {})""".format(max_nr_data_points)
    # array size parameter only - add some buffer -> *1.1
    output += """
C MNDATT IS NUMBER OF THEORETICAL DATA POINTS IN EACH BEAM
      PARAMETER(MNDATT = {})""".format(max_nr_data_points)
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
    # i3 = ff.FortranRecordWriter("I3")
    i1 = ff.FortranRecordWriter("I1")
    maxgen = rp.SEARCH_MAX_GEN
    if rp.TL_VERSION < 1.7:
        formatter = {'int': ff.FortranRecordWriter('I3'),
                     'gens': ff.FortranRecordWriter('I6'),
                     'rmut': ff.FortranRecordWriter('F7.4'),
                     'ctrl': ff.FortranRecordWriter('I3'),
                     }
        if maxgen > 999999:
            maxgen = 999999
        ctrl_width = 3
    else:
        formatter = {'int': ff.FortranRecordWriter('I5'),
                     'gens': ff.FortranRecordWriter('I7'),
                     'rmut': ff.FortranRecordWriter('F7.4'),
                     'ctrl': ff.FortranRecordWriter('I4'),
                     }
        if maxgen > 9999999:
            maxgen = 9999999
        ctrl_width = 4
    output += (formatter['int'].write([rp.indyPars]).ljust(16)
               + "number of independent parameters\n")
    f74 = ff.FortranRecordWriter('F7.4')
    output += (formatter['rmut'].write([rp.GAUSSIAN_WIDTH]).ljust(16)
               + "gaussian width control parameter RMUT\n")
    output += (formatter['int'].write([0]).ljust(16) + "initialisation for "
               "random number generator - 0: use system time, 1,...: use"
               " init\n")
    output += (formatter['int'].write([rp.R_FACTOR_TYPE]).ljust(16)
               + "1: use RPe -- 2: use R2\n")
    output += (formatter['int'].write([rp.SEARCH_BEAMS]).ljust(16)
               + "Optimization of which beam group do you want? "
               "(0=Aver,1=Int,2=Half)\n")
    if rp.TL_VERSION >= 1.71:
        outdata = 0
        if rp.PARABOLA_FIT["type"] != "none":
            outdata = 1
        output += (formatter['int'].write([outdata]).ljust(16)
                   + "Store configurations and write data.chem files for "
                   "parabola fit (0=false)\n")
    output += formatter['gens'].write([rp.output_interval]).ljust(16) + "output interval\n"
    output += (formatter['gens'].write([maxgen]).ljust(16)
               + "desired number of generations to be performed\n")
    output += (formatter['int'].write([astep]).ljust(16)
               + "area fraction step width (%)\n")
    output += ("SD.TL           name of search document file "
               "(max. 10 characters)\n")
    output += (formatter['int'].write([ndom]).ljust(16)
               + "Number of domains under consideration\n")
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
            formatter['int'].write([len(crp.search_atlist)]).ljust(16)
            + "Number of atomic sites in variation: Domain {}\n".format(k+1))
        displistcount = len(csl.displists)+1
        surfats = csl.getSurfaceAtoms(crp)
        for (i, at) in enumerate(crp.search_atlist):
            printvac = False
            output += (
                "------- Information about site {}: -----------------------"
                "-----------------------\n".format(i+1))
            surf = 1 if at in surfats else 0 # Flag that goes into variable NSURF used by search in GetInt
            output += (formatter['int'].write([surf]).ljust(16)
                       + "Surface (0/1)\n")
            if at.displist in csl.displists:
                dlind = csl.displists.index(at.displist) + 1
            else:
                dlind = displistcount
                displistcount += 1
            output += (formatter['int'].write([dlind]).ljust(16)
                       + "Atom number\n")
            output += (
                formatter['int'].write([len(at.known_deltas)]).ljust(16)
                + "No. of different files for Atom no. {}\n".format(i+1))
            for (j, deltafile) in enumerate(at.known_deltas):
                name = deltafile
                if frompath:  # need to get the file; if True frompath is Path
                    name = "D{}_".format(k+1) + deltafile
                    if len(name) > 15:
                        un = 1
                        while "D{}_DEL_{}".format(k+1, un) in uniquenames:
                            un += 1
                        name = "D{}_DEL_{}".format(k+1, un)
                        uniquenames.append(name)
                    try:
                        shutil.copy2(frompath / deltafile, name)
                    except Exception:
                        logger.error("Error getting Delta file {} for search"
                                     .format(os.path.relpath(
                                         frompath / deltafile)))
                        raise
                output += "****Information about file {}:\n".format(j+1)
                output += (name.ljust(16) + "Name of file {} (max. 15 "
                           "characters)\n".format(j+1))
                output += (formatter['int'].write([1]).ljust(16)
                           + "Formatted(0/1)\n")
                el = deltafile.split("_")[-2]
                if el.lower() == "vac":
                    geo = 1
                    vib = 0
                    types = 1
                    printvac = True
                else:
                    # geo0 = True
                    vib0 = True
                    for (mode, disp) in [(1, at.disp_geo), (2, at.disp_vib)]:
                        if el in disp:
                            dl = disp[el]
                        else:
                            dl = disp["all"]
                        if mode == 1:
                            geo = len(dl)
                            # if geo == 1 and np.linalg.norm(dl[0]) >= 1e-4:
                            #     geo0 = False
                        else:
                            vib = len(dl)
                            if vib == 1 and dl[0] != 0.:
                                vib0 = False
                    if vib == 1 and vib0:
                        vib = 0
                    # elif geo == 1 and geo0:
                    #     geo = 0
                    if geo > 0 and vib > 0:
                        types = 2
                    else:
                        types = 1
                output += (formatter['int'].write([types]).ljust(16)
                           + "Types of parameters in file {}\n".format(j+1))
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
                    output += (formatter['int'].write([vib]).ljust(16) +
                               "vibrational steps\n")
                    parcount += 1
                    info += (str(parcount).rjust(4) + ("P"+label).rjust(7)
                             + str(at.num).rjust(7) + el.rjust(7)
                             + "vib".rjust(7) + str(vib).rjust(7)
                             + constr["vib"].rjust(7) + "\n")
                    nsteps.append(vib)
                if geo > 0:
                    output += (formatter['int'].write([geo]).ljust(16) +
                               "geometrical steps\n")
                    parcount += 1
                    info += (str(parcount).rjust(4) + ("P"+label).rjust(7)
                             + str(at.num).rjust(7) + el.rjust(7)
                             + "geo".rjust(7) + str(geo).rjust(7)
                             + constr["geo"].rjust(7) + "\n")
                    nsteps.append(geo)
            output += "****concentration steps for site no. {}\n".format(i+1)
            occsteps = len(next(iter(at.disp_occ.values())))
            output += (formatter['int'].write([occsteps]).ljust(16)
                       + "no. of concentration steps - sum must equal 1 !\n")
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
                     + str(at.num).rjust(7) + "-".rjust(7) + "occ".rjust(7)
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
        elif os.path.isfile(os.path.join("SUPP", "control.chem")):
            controlpath = os.path.join("SUPP", "control.chem")
            logger.warning("No control.chem file found in working folder, "
                           "using SUPP/control.chem")
            rp.setHaltingLevel(1)
        else:
            logger.warning("No control.chem file found. Defaulting to random "
                           "starting configuration.")
            rp.SEARCH_START = "random"
            rp.setHaltingLevel(2)
    if rp.SEARCH_START == "control":
        # try reading control.chem
        # length = population +1 for header; +1 for empty line at the end
        control_chem_expected_length = rp.SEARCH_POPULATION + 2
        try:
            controllines = _read_control_chem(controlpath,
                                              control_chem_expected_length)
        except SearchIORaceConditionError:
            # Could not read last generation from control.chem
            # use backup from SD.TL instead
            if not rp.controlChemBackup:
                logger.error("Information in control.chem is incomplete "
                             "and no SD.TL data is stored. "
                             "Defaulting to random starting configuration.")
                rp.SEARCH_START = "random"
                rp.setHaltingLevel(2)
            else:
                logger.debug("Failed to read current configuration from "
                             "control.chem. "
                             "Loading last known configuraiton from stored "
                             "SD.TL data instead...")
                controllines = [s+"\n"
                                for s in rp.controlChemBackup.split("\n")]
            pass
        except OSError:
            logger.error("Error reading control.chem file. Defaulting to "
                         "random starting configuration.")
            rp.SEARCH_START = "random"
            rp.setHaltingLevel(2)

    if rp.SEARCH_START == "random":
        output += (formatter['int'].write([0]).ljust(16) +
                   "Certain start position (1) or random configuration (0)\n")
    else:
        output += (formatter['int'].write([1]).ljust(16) +
                   "Certain start position (1) or random configuration (0)\n")
        if rp.SEARCH_START == "control":
            if cull and rp.SEARCH_CULL:
                try:
                    ncull = rp.SEARCH_CULL.nr_individuals(rp.SEARCH_POPULATION)
                except ValueError:  # Too small population
                    ncull = 0
                    logger.warning(
                        "SEARCH_CULL parameter too large: would cull "
                        "entire population. Culling will be skipped."
                        )
                if any([sp.parabolaFit["min"] is not None
                        for sp in rp.searchpars]):
                    # replace one by predicted best
                    getPredicted = True
                else:
                    getPredicted = False
                nsurvive = rp.SEARCH_POPULATION - ncull
                clines = controllines[2:]
                csurvive = []
                if rp.SEARCH_CULL.type_.is_genetic or getPredicted:
                    # prepare readable clines
                    try:
                        csurvive = [readIntLine(s, width=ctrl_width)
                                    for s in clines[:nsurvive]]
                    except ValueError:
                        if rp.SEARCH_CULL.type_.is_genetic:
                            logger.warning(
                                "SEARCH_CULL: Failed to read old "
                                "configuration from control.chem, cannot run "
                                "genetic algorithm. Defaulting to cloning."
                                )
                            rp.setHaltingLevel(1)
                            csurvive = []
                for (i, line) in enumerate(clines):
                    if i < nsurvive:
                        output += line
                    elif (rp.SEARCH_CULL.type_.is_random or
                          (rp.SEARCH_CULL.type_.is_genetic and csurvive)
                          or getPredicted):
                        if getPredicted:
                            bc = None
                            if csurvive:
                                bc = csurvive[0]
                            nc = rp.getPredictConfig(
                                best_config=bc,
                                mincurv=rp.PARABOLA_FIT["mincurv"])
                            getPredicted = False
                        elif rp.SEARCH_CULL.type_.is_random:
                            nc = rp.getRandomConfig()
                        else:  # "genetic"
                            nc = rp.getOffspringConfig(csurvive)
                        ol = ""
                        for v in nc:
                            ol += formatter['ctrl'].write([v])
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
                    ol += formatter['ctrl'].write([(n+1)/2])
                output += ol + "\n"
        else:    # rp.SEARCH_START == "crandom"
            pop = [rp.getCenteredConfig()]
            for i in range(1, rp.SEARCH_POPULATION):
                pop.append(rp.getRandomConfig())
            for p in pop:
                ol = ""
                for v in p:
                    ol += formatter['ctrl'].write([v])
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

def _read_control_chem(control_chem_path,
                       expected_lines,
                       max_repeats=2000,
                       sleep_time=5):
    """
    Read the content of the control.chem file located at the specified path and
    return its lines.

    Parameters
    ----------
    control_chem_path : path-like
        The path to the control.chem file.
    expected_lines : int
        The expected number of lines in the control.chem file.
    max_repeats : int, optional
        The maximum number of repetitions to read the file before giving up.
        Defaults is 1000.
    sleep_time : int, optional
        The time in milliseconds to sleep between each attempt. Default is 1.

    Returns
    -------
    list of str
        The lines of the control.chem file.

    Raises
    ------
    ValueError
        If max_repeats is less than 1 or sleep_time is less than 0.
    RuntimeError
        If the number of lines in the control.chem file exceeds the
        expected_lines. This should not happen and indicates that expected_lines
        was set incorrectly.
    SearchIORaceConditionError
        If the complete control.chem file could not be read after the
        maximum number of repetitions.

    """
    if max_repeats < 1:
        raise ValueError("max_repeats must be >= 1")
    if sleep_time < 0:
        raise ValueError("sleep_time must be >= 0 ms")
    _control_chem_path = Path(control_chem_path)
    for repeat in range(max_repeats):
        with open(_control_chem_path, "r") as rf:
            control_lines = rf.readlines()
        n_control_lines = len(control_lines)
        if n_control_lines == expected_lines:
            logger.debug(f"Read complete control.chem file after {repeat+1} "
                         "read attempt(s).")
            return control_lines
        elif n_control_lines > expected_lines:
            raise RuntimeError(f"Expected at maximum {expected_lines} lines "
                               f"in control.chem, but found {n_control_lines}.")
        time.sleep(sleep_time/1000)  # in milliseconds
    raise SearchIORaceConditionError("Could not read complete "
                                     "control.chem file")

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
    for at in sl:
        at.storeOriState()
    for site in sl.sitelist:
        if site.oriState is None:
            tmp = copy.deepcopy(site)
            site.oriState = tmp

    sl.update_cartesian_from_fractional()
    uci = np.linalg.inv(sl.ucell)
    for at in sl:
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
    sl.collapse_fractional_coordinates()
    sl.update_cartesian_from_fractional()
    sl.update_layer_coordinates()
    # now update site occupations and vibrations:
    for site in sl.sitelist:
        siteats = [at for at in sl if at.site == site and not at.is_bulk]
        if not siteats: # site is only found in bulk
            continue
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
    tmpslab.sort_original()
    try:
        poscar.write(tmpslab, filename=fn, comments="all", silent=silent)
    except OSError:
        logger.error("Exception occurred while writing POSCAR_OUT" + suffix,
                     exc_info=rp.is_debug_mode)
        rp.setHaltingLevel(2)
    if not np.isclose(rp.SYMMETRY_CELL_TRANSFORM, np.identity(2)).all():
        tmpslab = sl.make_subcell(rp, rp.SYMMETRY_CELL_TRANSFORM)
        fn = "POSCAR_OUT_mincell" + suffix + "_" + rp.timestamp
        try:
            poscar.write(tmpslab, filename=fn, silent=silent)
        except OSError:
            logger.warning(
                "Exception occurred while writing POSCAR_OUT_mincell" + suffix,
                exc_info=rp.is_debug_mode)
    fn = "VIBROCC_OUT" + suffix + "_" + rp.timestamp
    try:
        writeVIBROCC(sl, rp, filename=fn, silent=silent)
    except Exception:
        logger.error("Exception occured while writing VIBROCC_OUT" + suffix,
                     exc_info=rp.is_debug_mode)
        rp.setHaltingLevel(2)
    return
