"""Section Superpos."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michael Riva (@michele-riva)'
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-11'
__license__ = 'GPLv3+'

import copy
import logging
import os
from pathlib import Path
import shutil
import subprocess

from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.files import iosuperpos
from viperleed.calc.files.beams import averageBeams
from viperleed.calc.files.beams import writeFdOut
from viperleed.calc.files.beams import writeOUTBEAMS
from viperleed.calc.files.displacements import readDISPLACEMENTS_block
from viperleed.calc.files.iorefcalc import readFdOut
from viperleed.calc.files.iosearch import readSDTL_blocks
from viperleed.calc.files.iosearch import readSDTL_end
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.checksums import validate_multiple_files

logger = logging.getLogger(__name__)


def superpos(sl, rp, subdomain=False, for_error=False, only_vary=None):
    """Runs the superpos calculation."""
    # check whether there is anything to evaluate
    if for_error:
        config = ((None, None),)
    elif not for_error and rp.searchResultConfig is not None:
        config = rp.searchResultConfig[0]
        only_vary = None  # Meant to use only for_error
    else:  # not for_error, but we don't have the search results handy
        config = _best_config_from_sdtl(rp)
        only_vary = None  # Meant to use only for_error

    if not config:
        # Errors in _best_config_from_sdtl, reported already
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
        contrin = iosuperpos.writeSuperposInput(sl, rp, config[0][1],
                                                for_error=for_error,
                                                only_vary=only_vary)
    except Exception:
        logger.error("Error getting input data for Superpos: ", exc_info=True)
        rp.setHaltingLevel(2)
        return

    if not contrin:
        logger.error("Error getting input data for Superpos: "
                     "writeSuperposInput returned empty. Cancelling "
                     "Superpos...")
        rp.setHaltingLevel(2)
        return

    logger.debug("Wrote Superpos input successfully")
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
        tl_source = rp.get_tenserleed_directory()
        tl_path = tl_source.path
        src_dir = tl_path / 'src'
        src_path = next(f for f in src_dir.glob('superpos*'))                   # TODO: StopIteration; Probably also want *.f* in the glob?
        shutil.copy2(src_path, src_path.name)
        lib_dir = tl_path / 'lib'
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
                                "superpos", rp.TL_VERSION)

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
            complete = subprocess.run(str(sposname.resolve()),                  # TODO: Perhaps nicer to have log be a TextIO stream?
                                      input=contrin, encoding="ascii",
                                      capture_output=True)
            out.write(complete.stdout)
            err_log = complete.stderr
    except Exception:
        logger.error("Error during Superpos calculation.")
        raise
    err_log = "\n".join(line for line in err_log.splitlines()
                        if ".  CORRECT TERMINATION" not in line)
    if err_log:
        logger.warning("Superpos output contained the following "
                       f"warnings/error messages:\n{err_log}")
    logger.info("Finished Superpos calculation. Processing files...")
    try:
        beams_and_spec = readFdOut(outname, for_error=for_error)
        # if for_error, theobeams will contain only the first variation
    except FileNotFoundError:
        logger.error(f"{outname} not found after superpos calculation.")
        raise
    except Exception:
        logger.error(f"Error reading {outname} after superpos calculation.")
        raise
    rp.theobeams["superpos"], rp.superpos_specout = beams_and_spec

    if not for_error:
        _write_fitbeams(rp)

    # Rename and move files
    try:
        os.rename('PARAM', 'superpos-PARAM')
    except Exception:
        logger.warning("Failed to rename superpos input file PARAM to "
                       "superpos-PARAM")
    try:
        os.rename('DOC', 'superpos-DOC')
    except Exception:
        pass


def superpos_domains(rp, configs):
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
    for (percent, params), dp in zip(configs, rp.domainParams):
        percentages.append(percent)
        dp.rp.searchResultConfig = [[(100, params)]]
        logger.info(f"Running superpos calculation for domain {dp.name}")
        try:
            os.chdir(dp.workdir)
            superpos(dp.sl, dp.rp, subdomain=True)
        except Exception:
            logger.error("Error while running superpos "
                         f"calculation for domain {dp.name}")
            raise
        finally:
            os.chdir(home)

    logger.info("Getting weighted average over domain beams...")
    rp.theobeams["superpos"] = averageBeams(
        [dp.rp.theobeams["superpos"] for dp in rp.domainParams],
        weights=percentages
        )
    _write_fitbeams(rp)

    try:
        rp.superpos_specout = writeFdOut(rp.theobeams["superpos"], rp.beamlist,
                                         filename="superpos-spec.out",
                                         header=rp.systemName)
    except Exception:
        logger.error("Error writing averaged superpos-spec.out for R-factor "
                     "calculation.", exc_info=rp.is_debug_mode)


def _best_config_from_sdtl(rp):                                                 # TODO: more doc
    """Return the best configuration from the most recent search."""
    # check for an SD.TL file
    cwd = Path()
    sdtl_files = (
        f for f in (cwd/"SD.TL", cwd/{DEFAULT_OUT}/"SD.TL") if f.is_file()
        )
    sdtl_file = next(sdtl_files, None)
    if not sdtl_file:
        logger.error("Superpos: Found no stored results from recent "
                     "search and no SD.TL file. Cancelling...")
        rp.setHaltingLevel(2)
        return None

    try:
        sdtl = readSDTL_end(filename=str(sdtl_file),
                            n_expect=rp.SEARCH_POPULATION)
    except Exception:                                                           # TODO: better exception
        logger.error("Superpos: Error reading SD.TL")
        rp.setHaltingLevel(2)
        return None

    if sdtl is not None:
        sdtl = readSDTL_blocks("\n".join(sdtl),
                               whichR=rp.SEARCH_BEAMS,
                               n_expect=rp.SEARCH_POPULATION)
    if not sdtl:
        logger.error("Superpos: No data found in SD.TL")
        rp.setHaltingLevel(2)
        return None

    try:
        return sdtl[0][2][0]  # Most recent block, config, best only
    except IndexError:
        logger.error("Superpos: Failed to read best configuration from "
                     "SD.TL")
        rp.setHaltingLevel(2)
        return None


def _write_fitbeams(rp):                                                        # TODO: Could probably be made more general by passing in which rp beams, a base filename, and a Section (or even rp and a Section)
    """Write FITBEAMS and its _norm version."""
    try:
        writeOUTBEAMS(rp.theobeams["superpos"], filename="FITBEAMS.csv")
        theobeams_norm = copy.deepcopy(rp.theobeams["superpos"])
        for beams in theobeams_norm:
            beams.normMax()
        writeOUTBEAMS(theobeams_norm, filename="FITBEAMS_norm.csv")
    except Exception:
        logger.error("Error writing FITBEAMS after superpos calculation.",
                     exc_info=rp.is_debug_mode)
