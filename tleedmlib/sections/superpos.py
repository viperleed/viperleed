# -*- coding: utf-8 -*-

"""
Created on Aug 11 2020

@author: Florian Kraushofer

Tensor LEED Manager section Superpos
"""

import copy
import logging
import os
from pathlib import Path
import shutil
import subprocess

from viperleed.tleedmlib import leedbase
from viperleed.tleedmlib.checksums import validate_multiple_files
from viperleed.tleedmlib.files.beams import (writeOUTBEAMS, averageBeams,
                                             writeFdOut)
from viperleed.tleedmlib.files.displacements import readDISPLACEMENTS_block
from viperleed.tleedmlib.files.iorefcalc import readFdOut
from viperleed.tleedmlib.files.iosearch import readSDTL_end, readSDTL_blocks
import viperleed.tleedmlib.files.iosuperpos as tl_io

logger = logging.getLogger("tleedm.superpos")


def superpos(sl, rp, subdomain=False, for_error=False, only_vary=None):
    """Runs the superpos calculation."""
    # check whether there is anything to evaluate
    if rp.searchResultConfig is not None and not for_error:
        config = rp.searchResultConfig[0]
    elif not for_error:
        config = None
        # check for an SD.TL file
        sdtl = None
        if os.path.isfile("SD.TL"):
            try:
                sdtl = readSDTL_end(filename="SD.TL",
                                    n_expect=rp.SEARCH_POPULATION)
            except Exception:
                logger.error("Superpos: Error reading SD.TL")
                rp.setHaltingLevel(2)
                return
        elif os.path.isfile(os.path.join("OUT", "SD.TL")):
            try:
                sdtl = readSDTL_end(filename=os.path.join("OUT", "SD.TL"),
                                    n_expect=rp.SEARCH_POPULATION)
            except Exception:
                logger.error("Superpos: Error reading SD.TL")
                rp.setHaltingLevel(2)
                return
        else:
            logger.error("Superpos: Found no stored results from recent "
                         "search and no SD.TL file. Cancelling...")
            rp.setHaltingLevel(2)
            return
        if sdtl is None:
            logger.error("Superpos: No data found in SD.TL")
            rp.setHaltingLevel(2)
            return
        sdtlContent = readSDTL_blocks("\n".join(sdtl),
                                      whichR=rp.SEARCH_BEAMS,
                                      n_expect=rp.SEARCH_POPULATION)
        if not sdtlContent:
            logger.error("Superpos: No data found in SD.TL")
            rp.setHaltingLevel(2)
            return
        try:
            config = sdtlContent[0][2][0]  # first block, config, best only
        except IndexError:
            logger.error("Superpos: Failed to read best configuration from "
                         "SD.TL")
            rp.setHaltingLevel(2)
            return

    if rp.domainParams:
        superpos_domains(rp, config)
        return

    # read DISPLACEMENTS block and fetch deltas
    if not rp.disp_block_read:
        readDISPLACEMENTS_block(rp, sl, rp.disp_blocks[rp.search_index])
        rp.disp_block_read = True
    if not any(ind in rp.runHistory for ind in [2, 3]) and not for_error:
        leedbase.getDeltas(rp.TENSOR_INDEX, required=True)
    # make sure search parameters are initialized
    if 3 not in rp.runHistory and not subdomain and not for_error:
        logger.debug("Superpos calculation executed without search. "
                     "Search parameters will be inferred from input files.")
        try:
            rp.generateSearchPars(sl)
        except Exception:
            logger.error("Error getting search parameters. Superpos will "
                         "stop.")
            return
    # now we have configuration and parameters, create input:
    contrin = ""
    try:
        if not for_error:
            contrin = tl_io.writeSuperposInput(sl, rp, config[0][1])
        else:
            contrin = tl_io.writeSuperposInput(sl, rp, None, for_error=True,
                                               only_vary=only_vary)
        logger.debug("Wrote Superpos input successfully")
    except Exception:
        logger.error("Error getting input data for Superpos: ", exc_info=True)
        rp.setHaltingLevel(2)
        return
    if contrin == "":
        logger.error("Error getting input data for Superpos: "
                     "writeSuperposInput returned empty. Cancelling "
                     "Superpos...")
        rp.setHaltingLevel(2)
        return
    # if execution is suppressed, stop here:
    if rp.SUPPRESS_EXECUTION:
        logger.warning("SUPPRESS_EXECUTION parameter is on. Superpos "
                       "calculation will not proceed. Stopping...")
        rp.setHaltingLevel(3)
        return

    if not rp.FORTRAN_COMP[0]:
        rp.getFortranComp()
    # Get FORTRAN files                                                         # TODO: use CompileTask subclass (Issue #43)
    try:
        tldir = leedbase.getTLEEDdir(home=rp.sourcedir, version=rp.TL_VERSION)
        if not tldir:
            raise RuntimeError("TensErLEED code not found.")
        src_dir = tldir / 'src'
        src_path = next(f for f in src_dir.glob('superpos*'))                   # TODO: StopIteration; Probably also want *.f* in the glob?
        shutil.copy2(src_path, src_path.name)
        lib_dir = tldir / 'lib'
        lib_path = next(f for f in lib_dir.glob('lib.superpos*'))               # TODO: StopIteration; Probably also want *.f* in the glob?
        shutil.copy2(lib_path, lib_path.name)
        globalname = "GLOBAL"
        shutil.copy2(src_dir / globalname, globalname)
    except Exception:
        logger.error("Error getting TensErLEED files for superpos: ")
        raise

    # Validate checksums
    if not rp.TL_IGNORE_CHECKSUM:
        files_to_check = (lib_path,
                          src_path,
                          src_dir / globalname)
        validate_multiple_files(files_to_check, logger,
                                "superpos", rp.TL_VERSION_STR)

    # Compile FORTRAN files
    sposname = Path(f"superpos-{rp.timestamp}")
    logger.info("Compiling fortran input files...")

    compile_log = "compile-superpos.log"
    try:
        leedbase.fortran_compile(f"{rp.FORTRAN_COMP[0]} -o",
                                 f"{sposname} {src_path.name} {lib_path.name}",
                                 rp.FORTRAN_COMP[1],
                                 logname=compile_log)
    except Exception:
        leedbase.copy_compile_log(rp, Path(compile_log), "superpos-compile")
        logger.error("Error compiling fortran files: ", exc_info=True)
        raise

    logger.info("Starting Superpos calculation...")
    outname = "superpos-spec.out"
    err_log = ""
    try:
        with open(outname, "w") as out:
            complete = subprocess.run(sposname.resolve(),                       # TODO: Perhaps nicer to have log be a TextIO stream?
                                      input=contrin, encoding="ascii",
                                      capture_output=True)
            out.write(complete.stdout)
            err_log = complete.stderr
    except Exception:
        logger.error("Error during Superpos calculation.")
        raise
    err_log = "\n".join([line for line in err_log.split("\n")
                         if ".  CORRECT TERMINATION" not in line])
    if err_log:
        logger.warning("Superpos output contained the following warnings/"
                       "error messages:\n"+err_log)
    logger.info("Finished Superpos calculation. Processing files...")
    try:
        rp.theobeams["superpos"], rp.superpos_specout = readFdOut(
            outname, for_error=for_error)
        # if for_error, theobeams will contain only the first variation
    except FileNotFoundError:
        logger.error(outname + " not found after superpos calculation.")
        raise
    except Exception:
        logger.error("Error reading " + outname + " after superpos "
                     " calculation.")
        raise
    if not for_error:
        try:
            writeOUTBEAMS(rp.theobeams["superpos"], filename="FITBEAMS.csv")
            theobeams_norm = copy.deepcopy(rp.theobeams["superpos"])
            for b in theobeams_norm:
                b.normMax()
            writeOUTBEAMS(theobeams_norm, filename="FITBEAMS_norm.csv")
        except Exception:
            logger.error("Error writing FITBEAMS after superpos calculation.")
    # rename and move files
    try:
        os.rename('PARAM', 'superpos-PARAM')
    except Exception:
        logger.warning("Failed to rename superpos input file PARAM to "
                       "superpos-PARAM")
    try:
        os.rename('DOC', 'superpos-DOC')
    except Exception:
        pass
    return


def superpos_domains(rp, config):
    """Runs the normal superpos function for each subdomain, collects the
    results and averages over the beams to get an overall result."""
    # make sure search parameters are initialized
    if 3 not in rp.runHistory:
        logger.debug("Superpos calculation executed without search. "
                     "Search parameters will be inferred from input files.")
        try:
            rp.generateSearchPars(None)
        except Exception:
            logger.error("Error getting search parameters. Superpos will "
                         "stop.")
            return
    home = os.getcwd()
    percentages = []
    for (i, dp) in enumerate(rp.domainParams):
        (percent, params) = config[i]
        percentages.append(percent)
        dp.rp.searchResultConfig = [[(100, params)]]
        logger.info("Running superpos calculation for domain {}"
                    .format(dp.name))
        try:
            os.chdir(dp.workdir)
            superpos(dp.sl, dp.rp, subdomain=True)
        except Exception:
            logger.error("Error while running superpos calculation for domain "
                         "{}".format(dp.name))
            raise
        finally:
            os.chdir(home)
    logger.info("Getting weighted average over domain beams...")
    rp.theobeams["superpos"] = averageBeams(
        [dp.rp.theobeams["superpos"] for dp in rp.domainParams],
        weights=percentages)
    try:
        writeOUTBEAMS(rp.theobeams["superpos"], filename="FITBEAMS.csv")
        theobeams_norm = copy.deepcopy(rp.theobeams["superpos"])
        for b in theobeams_norm:
            b.normMax()
        writeOUTBEAMS(theobeams_norm, filename="FITBEAMS_norm.csv")
    except Exception:
        logger.error("Error writing FITBEAMS after superpos calculation.",
                     exc_info=rp.LOG_DEBUG)
    try:
        rp.superpos_specout = writeFdOut(rp.theobeams["superpos"], rp.beamlist,
                                         filename="superpos-spec.out",
                                         header=rp.systemName)
    except Exception:
        logger.error("Error writing averaged superpos-spec.out for R-factor "
                     "calculation.", exc_info=rp.LOG_DEBUG)
    return
