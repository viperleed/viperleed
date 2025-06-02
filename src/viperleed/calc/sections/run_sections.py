"""Module run_sections of viperleed.calc.sections.

Wrapper functions for running a single section as well as
a section loop.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-19'
__license__ = 'GPLv3+'

import logging
import numpy as np
import os
import shutil

from viperleed.calc.classes.state_recorder import CalcStateRecorder
from viperleed.calc.classes.rparams.defaults import NO_VALUE
from viperleed.calc.constants import SKIP_IN_DOMAIN_MAIN
from viperleed.calc.files import beams as iobeams
from viperleed.calc.files import parameters
from viperleed.calc.files.displacements import readDISPLACEMENTS
from viperleed.calc.files.phaseshifts import readPHASESHIFTS
from viperleed.calc.files.vibrocc import readVIBROCC, writeVIBROCC
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.lib.time_utils import ExecutionTimer
from viperleed.calc.sections import deltas
from viperleed.calc.sections import errorcalc
from viperleed.calc.sections import fd_optimization
from viperleed.calc.sections import initialization                              # TODO: perhaps use CalcSection
from viperleed.calc.sections import refcalc
from viperleed.calc.sections import rfactor
from viperleed.calc.sections import search
from viperleed.calc.sections import superpos
from viperleed.calc.sections.cleanup import cleanup, move_oldruns


logger = logging.getLogger(__name__)
R_MAX = np.inf


def run_section(index, sl, rp):
    """Run a specific viperleed.calc section.

    Parameters
    ----------
    index : int
        Index of the section, see sectionNames variable below.
    sl : Slab
        Slab object containing atom information.
    rp : Rparams
        The run parameters.

    Returns
    -------
    None.
    """
    sectionNames = {0: "INITIALIZATION",
                    1: "REFERENCE CALCULATION",
                    2: "DELTA-AMPLITUDES",
                    3: "SEARCH",
                    11: "R-FACTOR CALCULATION",
                    12: "R-FACTOR CALCULATION",
                    31: "SUPERPOS",
                    5:  "ERROR CALCULATION",
                    6:  "FULL-DYNAMIC OPTIMIZATION"}
    # files that need to be there for the different parts to run
    requiredFiles = {0: ["POSCAR", "PARAMETERS", "VIBROCC", "IVBEAMS"],
                     1: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC"],
                     2: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS"],
                     3: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS", "EXPBEAMS"],
                     11: ["BEAMLIST", "PARAMETERS", "IVBEAMS", "EXPBEAMS"],
                     12: ["BEAMLIST", "PARAMETERS", "IVBEAMS", "EXPBEAMS"],
                     31: ["BEAMLIST", "POSCAR", "PARAMETERS", "IVBEAMS",
                          "VIBROCC", "DISPLACEMENTS"],
                     5: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "DISPLACEMENTS", "EXPBEAMS"],
                     6: ["BEAMLIST", "PHASESHIFTS", "POSCAR", "PARAMETERS",
                         "IVBEAMS", "VIBROCC", "EXPBEAMS"]}

    checkfiles = requiredFiles[index][:]
    o = "\nSTARTING SECTION: "+sectionNames[index]
    if index == 3 and rp.disp_blocks and rp.disp_blocks[rp.search_index][1]:
        o += " "+rp.disp_blocks[rp.search_index][1]  # displacement block name
    if rp.domainParams or rp.DOMAINS:
        o += " (DOMAINS)"
        for fn in SKIP_IN_DOMAIN_MAIN:
            try:
                checkfiles.remove(fn)
            except Exception:
                pass
    logger.info(o)
    since_section_started = ExecutionTimer()
    rp.runHistory.append(index)
    for dp in rp.domainParams:
        dp.rpars.runHistory = rp.runHistory
    i = 0
    while i < len(checkfiles):
        filename = checkfiles[i]
        ignoreError = False
        if rp.fileLoaded[filename]:
            i += 1
            continue
        # try loading files
        if filename == "EXPBEAMS":
            rp.try_loading_expbeams_file()
            if index != 0:
                iobeams.checkEXPBEAMS(sl, rp)
        elif filename == "IVBEAMS":
            try:
                rp.ivbeams = iobeams.readIVBEAMS()
                rp.ivbeams_sorted = False
                rp.fileLoaded["IVBEAMS"] = True
            except FileNotFoundError:
                if (os.path.isfile("EXPBEAMS")
                        or os.path.isfile("EXPBEAMS.csv")):
                    checkfiles.insert(i+1, "EXPBEAMS")
                    logger.warning("IVBEAMS file not found. Will attempt "
                                   "generating IVBEAMS from EXBEAMS.")
                    ignoreError = True
                else:
                    logger.error("Neither IVBEAMS nor EXPBEAMS file "
                                 "found.")
            except Exception as e:
                logger.error("Error while reading required file IVBEAMS",
                             exc_info=(type(e) != FileNotFoundError))
        elif filename == "BEAMLIST":
            try:
                rp.beamlist = iobeams.readBEAMLIST()
                rp.fileLoaded["BEAMLIST"] = True
            except Exception as e:
                logger.error("Error while reading required file "
                             "BEAMLIST", exc_info=(type(e) !=
                                                   FileNotFoundError))
        elif filename == "VIBROCC":
            changeVIBROCC = False
            try:
                changeVIBROCC = readVIBROCC(rp, sl)
                rp.fileLoaded["VIBROCC"] = True
            except Exception as exc:
                _not_found = not isinstance(exc, FileNotFoundError)
                logger.error("Error while reading required file VIBROCC",
                             exc_info=_not_found)
            sl.full_update(rp)
            if changeVIBROCC:
                if os.path.isfile("VIBROCC"):
                    os.rename("VIBROCC", "VIBROCC_user")
                    rp.manifest.add("VIBROCC_user")
                    logger.info(
                        "VIBROCC file was modified with automatically "
                        "generated vibration amplitudes."
                        )
                writeVIBROCC(sl, 'VIBROCC')
                rp.files_to_out.add('VIBROCC')
                shutil.copy2('VIBROCC', 'VIBROCC_generated')
            if rp.T_EXPERIMENT is not None:
                parameters.comment_out(rp, "T_EXPERIMENT")
            if rp.T_DEBYE is not None:
                parameters.comment_out(rp, "T_DEBYE")
            if len(rp.VIBR_AMP_SCALE) > 0:
                parameters.comment_out(rp, "VIBR_AMP_SCALE")
        elif filename == "PHASESHIFTS":
            try:
                (rp.phaseshifts_firstline, rp.phaseshifts,
                 newpsGen, newpsWrite) = readPHASESHIFTS(sl, rp)
                if newpsGen:
                    logger.error(
                        "PHASESHIFTS file generation is only supported "
                        "during initialization. Stopping execution...")
                    raise RuntimeError("Inconsistent PHASESHIFT file")
                elif newpsWrite:
                    logger.warning(
                        "Writing a new PHASESHIFTS file is "
                        "only supported during initialization. The "
                        "data in the provided file will be used, but "
                        "running the initialization is recommended.")
                    rp.fileLoaded["PHASESHIFTS"] = True
                else:
                    rp.fileLoaded["PHASESHIFTS"] = True
            except Exception as e:
                logger.error("Error while reading required file "
                             "PHASESHIFTS", exc_info=(type(e) !=
                                                      FileNotFoundError))
        elif filename == "DISPLACEMENTS":
            try:
                readDISPLACEMENTS(rp)
                rp.fileLoaded["DISPLACEMENTS"] = True
            except Exception as e:
                logger.error("Error while reading required file "
                             "DISPLACEMENTS", exc_info=(type(e) !=
                                                        FileNotFoundError))
        if not rp.fileLoaded[filename] and not ignoreError:
            # and if that didn't work, stop:
            logger.error(f"Step '{sectionNames[index]}' requires "
                         f"file {filename}. Stopping execution...")
            rp.setHaltingLevel(3)
            return
        i += 1
    try:
        if index == 0:
            initialization.initialization(sl, rp)
        elif index == 1:
            refcalc.refcalc(sl, rp)
        elif index in [11, 12]:
            rfactor.rfactor(sl, rp, index)
        elif index == 2:
            deltas.deltas(sl, rp)
        elif index == 3:
            search.search(sl, rp)
        elif index == 31:
            superpos.superpos(sl, rp)
        elif index == 5:
            errorcalc.errorcalc(sl, rp)
        elif index == 6:
            fd_optimization.fd_optimization(sl, rp)
    except Exception:
        logger.error(f"Error in section {sectionNames[index]}")
        raise
    logger.info(
        f'Finishing section at {DateTimeFormat.TIME.now()}. '
        f'Section took {since_section_started.how_long(as_string=True)}.'
        )
    if index == 1:
        rp.last_refcalc_time = since_section_started.how_long()


def _check_exceeds_tl_limit(rpars, slab):
    """Return whether any atom in `slab` is too far from its refcalc state."""
    return any(rpars.MAX_TL_DISPLACEMENT.is_too_far(atom) for atom in slab)


def section_loop(rp, sl):
    """
    Executes sections as specified in rp.RUN, may loop if required.

    Parameters
    ----------
    rp : Rparams
        The run parameters.
    sl : Slab
        Slab object containing atom information.

    Returns
    -------
    exit_code : int
        0: exit without errors.
        1: clean exit through KeyboardInterrupt
        3: exit due to Exception during main loop
    state_recorder : CalcStateRecorder
        A record of the states of sl and rp at the end of each
        executed section.
    """
    # Create a record for states at end of each section
    state_recorder = CalcStateRecorder()

    sectionorder = [0, 1, 6, 11, 2, 3, 31, 12, 4, 5]
    search_loop_R = {}   # dict of R for each loop
    initHalt = False
    while rp.RUN:
        try:
            sec = rp.RUN.pop(0)
            if rp.runHistory and (sectionorder.index(sec)
                                  < sectionorder.index(rp.runHistory[-1])):
                logger.info("\nExecution repeats. Moving old output to "
                            "workhistory folder.")
                try:
                    move_oldruns(rp)
                except Exception:
                    logger.warning(
                        "Exception while trying to clean up earlier segments. "
                        "Program will proceed, but old files may be lost.",
                        exc_info=True
                        )
            run_section(sec, sl, rp)
            if rp.domainParams and sl is None:
                sl = rp.pseudoSlab
            if rp.domainParams:
                rp.setHaltingLevel(max(d.rpars.halt for d in rp.domainParams))

            # record state to the state recorder
            state_recorder.record(sl, rp, sec)

            # Decide how to proceed
            next_section = next(iter(rp.RUN), None)
            if (sec == 0 and not rp.domainParams and not sl.preprocessed
                    and rp.HALTING <= 2 and len(rp.RUN) > 0):
                logger.info("Initialization finished. Execution will stop. "
                            "Please check whether comments in POSCAR are "
                            "correct, then restart.")
                rp.checklist.append("Check whether comments in "
                                    "POSCAR are correct")
                rp.setHaltingLevel(2)
                initHalt = True
            elif (sec == 1 and rp.fileLoaded["EXPBEAMS"]):
                if (next_section != 11 and          # r-factor after refcalc
                        (not rp.domainParams or 3 in rp.runHistory)):
                    rp.RUN.insert(0, 11)
            elif (sec == 3 and rp.fileLoaded["EXPBEAMS"]):
                if next_section != 31:  # superpos after search
                    rp.RUN.insert(0, 31)
            elif sec == 31 and rp.fileLoaded["EXPBEAMS"]:
                if next_section != 12:   # r-factor after superpos
                    rp.RUN.insert(0, 12)
            elif sec == 12 and not rp.STOP:
                # check for max. displacement condition:
                exceeds_tl_limit = _check_exceeds_tl_limit(rp, sl)
                # check for loops:
                loops = [t for t in rp.disp_loops if t[1] == rp.search_index]
                for loop in sorted(loops, reverse=True):
                    # At least one loop ends at the index we're at now.
                    # Starting at the deepest loop, see if there was any
                    # improvement since it has last been checked. If so,
                    # repeat it; if not, go one loop out and keep checking.
                    # Automatically repeat deepest loop if exceeds_tl_limit.
                    improved = rp.last_R < search_loop_R.get(loop, R_MAX)
                    if improved or exceeds_tl_limit:
                        # Loop back
                        search_loop_R[loop] = rp.last_R
                        rp.search_index = loop[0]
                        logger.info("Search loop: repeating at block "
                                    + rp.disp_blocks[rp.search_index][1])
                        break
                else:
                    # We arrive here if there were no loops, or no loops have
                    # improved. Go to next block.
                    rp.search_index += 1
                    if loops:
                        o = "Search loop ends."
                        if len(rp.disp_blocks) > rp.search_index:
                            o += (" Continuing at block "
                                  + rp.disp_blocks[rp.search_index][1])
                        logger.info(o)
                for dp in rp.domainParams:
                    dp.rpars.search_index = rp.search_index
                stop_search = False
                if len(rp.disp_blocks) > rp.search_index and exceeds_tl_limit:
                    if rp.MAX_TL_DISPLACEMENT.action == 'ignore':
                        logger.warning(
                            'Displacements exceed MAX_TL_DISPLACEMENT, but '
                            'actions are disabled. Calculation will '
                            'proceed, please check results carefully.'
                            )
                    elif rp.MAX_TL_DISPLACEMENT.action == 'stop':
                        logger.info(
                            'Displacements exceed MAX_TL_DISPLACEMENT. '
                            'Search will stop.'
                            )
                        stop_search = True
                    elif rp.MAX_TL_DISPLACEMENT.action == 'refcalc':
                        maxtime = rp.MAX_TL_DISPLACEMENT.max_duration
                        if (maxtime != NO_VALUE
                                and rp.last_refcalc_time is not None):
                            if rp.last_refcalc_time > maxtime:
                                logger.info(
                                    'Displacements exceed '
                                    'MAX_TL_DISPLACEMENT, and reference '
                                    'calculations take longer than the '
                                    'specified limit time. Search will stop.'
                                    )
                                stop_search = True
                if len(rp.disp_blocks) > rp.search_index and not stop_search:
                    # there are more disp_blocks to do; append another search
                    rp.resetSearchConv()
                    for dp in rp.domainParams:
                        dp.rpars.resetSearchConv()
                    if (exceeds_tl_limit
                            and rp.MAX_TL_DISPLACEMENT.action == 'refcalc'):
                        logger.info(
                            'Displacements exceed MAX_TL_DISPLACEMENT. '
                            'A new reference calculation will be performed '
                            'before the next search block.'
                            )
                        search_loop_R = {}  # ignore all previously known R
                        rp.RUN = [1, 2, 3] + rp.RUN
                    else:
                        if rp.RUN[:2] != [2, 3]:
                            rp.RUN = [2, 3] + rp.RUN
                        # call restoreOriState *only* if we don't insert a
                        #  refcalc; otherwise the refcalc would run at the
                        #  original position!
                        if not rp.domainParams:
                            sl.restoreOriState()
                        for domain in rp.domainParams:
                            domain.slab.restoreOriState()
                else:
                    # The current instance of 'RUN = 2, 3' is done, all blocks
                    # in DISPLACEMENTS have been handled. However, we may want
                    # to do another search later (e.g. 'RUN = 1-3 1-3 1').
                    rp.search_index = 0
                    for domain in rp.domainParams:
                        domain.rpars.search_index = 0
        except KeyboardInterrupt:
            logger.warning("Stopped by keyboard interrupt, attempting "
                           "clean exit...")
            cleanup(rp)
            return 1, state_recorder
        except Exception:
            logger.error("Exception during viperleed.calc execution: ",
                         exc_info=True)
            cleanup(rp)
            return 3, state_recorder
        if rp.halt >= rp.HALTING:
            if not initHalt:
                logger.info(
                    "# An exception occurred that meets the halting "
                    "criteria defined by the HALTING parameter. Execution "
                    "will stop, check log for warnings and errors."
                    )
            break
        parameters.update(rp)  # Look for a user STOP
        if rp.RUN and rp.STOP and not rp.RUN[0] in [11, 12, 31]:
            logger.info("# Stopped by user STOP command.")
            break

    logger.debug("End of section loop.")
    disp_ranges_str = '\n\t'.join(str(at.disp_ranges) for at in sl)
    logger.log(1, f'Total ranges of all displacements:\n{disp_ranges_str}')     # TODO: Consider deleting, it's not even really true - outside of the loop, this is only the *final* displacement ranges.
    cleanup(rp)
    return 0, state_recorder
