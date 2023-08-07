"""Package calc of viperleed.

This package contains all the main functionality of the TensErLEED
manager part of ViPErLEED (tleedm), i.e., the part of ViPErLEED
dedicated to the calculation of theoretical I(V) curves as well
as the optimization of structural models such that their calculated
I(V) curves best fit experimental ones.

Modules
-------
base:
    Non-LEED-related classes and functions used throughout tleedmlib.
beamgen:
    Handle input/output to the FORTRAN program that generates the list
    of beams used internally by TensErLEED
checksums:
    Functions used to prevent user tinkering with FORTRAN source files.
    The actual checksums are contained in _checksums.dat and should
    only be edited by expert contributors.
leedbase:
    LEED- and TensErLEED-related base functions used in tleedmlib
periodic_table:
    Collection of basi atom data
psgen:
    Handle input/output to the FORTRAN program(s) that generate atomic
    phase-shifts and information concerning the energy dependence of
    the inner potential of the crystal
symmetry:
    Function for detecting symmetry of 2D- and 3D-periodic slabs

Packages
--------
classes:
    Definition of LEED-related classes used throughout tleedmlib
files:
    Input/output handling for TensErLEED sections
sections:
    Functionality to run the various, logically different parts of
    TensErLEED
wrapped:
    Python extensions written in FORTRAN
"""
from pathlib import Path
import logging
import os
import shutil
import time

import viperleed
from viperleed import GLOBALS
from viperleed.calc.base import CustomLogFormatter
from viperleed.calc.classes import rparams
from viperleed.calc.files.parameter_errors import ParameterError
from viperleed.calc.files.parameters import (readPARAMETERS,
                                                  interpretPARAMETERS)
from viperleed.calc.files.poscar import readPOSCAR
from viperleed.calc.leedbase import getMaxTensorIndex
from viperleed.calc.sections._sections import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import prerun_clean, cleanup
from viperleed.calc.sections.run_sections import section_loop

__authors__ = ["Florian Kraushofer (@fkraushofer)",
               "Alexander M. Imre (@amimre)",
               "Michele Riva (@michele-riva)"]

logger = logging.getLogger("tleedm")

def run_tleedm(system_name=None,
               console_output=True,
               slab=None,
               preset_params={},
               source=Path(),
               override_log_level=None):
    """
    Runs the TensErLEED Manager. By default, a PARAMETERS and a POSCAR file
    are expected, but can be replaced by passing the 'slab' and/or
    'present_params' kwargs.
    Parameters
    ----------
    system_name : str, optional
        Used as a comment in some output file headers
    console_output : bool, optional
        If False, will not add a logging.StreamHandler. Output will only be
        printed to the log file.
    slab : Slab, optional
        Start from a pre-existing slab, instead of reading from POSCAR.
    preset_params : dict, optional
        Parameters to add to the Rparams object after PARAMETERS has been read.
        Keys should be attributes of Rparam. Values in preset_params will
        overwrite values read from the PARAMETERS file, if present in both. If
        no PARAMETERS file is read, parameters will be read exclusively from
        present_params.
    source : str, optional
        Path where the 'tensorleed' directory can be found, which contains all
        the TensErLEED source code.
    Returns
    -------
    int
        0: exit without errors.
        1: clean exit through KeyboardInterrupt
        2: exit due to Exception before entering main loop
        3: exit due to Exception during main loop
    """
    os.umask(0)
    # start logger, write to file:
    timestamp = time.strftime("%y%m%d-%H%M%S", time.localtime())
    log_name = 'tleedm-'+timestamp+'.log'
    logger = logging.getLogger("tleedm")
    logger.setLevel(logging.INFO)
    logFormatter = CustomLogFormatter()
    fileHandler = logging.FileHandler(log_name, mode="w")
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    if console_output:
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        logger.addHandler(consoleHandler)
    logger.info("Starting new log: " + log_name + "\nTime of execution (UTC): "
                + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    logger.info("This is ViPErLEED version " + GLOBALS["version"] + "\n")
    logger.info("! THIS VERSION IS A PRE-RELEASE NOT MEANT FOR PUBLIC "
                "DISTRIBUTION !\n")
    tmp_manifest = ["SUPP", "OUT", log_name]
    try:
        rp = readPARAMETERS()
    except FileNotFoundError:
        if not preset_params:
            logger.error("No PARAMETERS file found, and no preset parameters "
                         "passed. Execution will stop.")
            cleanup(tmp_manifest)
            return 2
        rp = rparams.Rparams()
    except Exception:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmp_manifest)
        return 2
    # check if this is going to be a domain search
    domains = False
    if "DOMAIN" in rp.readParams:
        domains = True
    if domains:  # no POSCAR in main folder for domain searches
        slab = None
    elif slab is None:
        poscar_file = Path("POSCAR")
        if poscar_file.is_file():
            logger.info("Reading structure from file POSCAR")
            try:
                slab = readPOSCAR(filename=str(poscar_file.resolve()))
            except Exception:
                logger.error("Exception while reading POSCAR", exc_info=True)
                cleanup(tmp_manifest)
                return 2
        else:
            logger.error("POSCAR not found. Stopping execution...")
            cleanup(tmp_manifest)
            return 2
        if not slab.preprocessed:
            logger.info("The POSCAR file will be processed and overwritten. "
                        "Copying the original POSCAR to POSCAR_user...")
            try:
                shutil.copy2(poscar_file, "POSCAR_user")
                tmp_manifest.append("POSCAR_user")
            except Exception:
                logger.error("Failed to copy POSCAR to POSCAR_user. Stopping "
                             "execution...")
                cleanup(tmp_manifest)
                return 2
    try:
        # interpret the PARAMETERS file
        interpretPARAMETERS(rp, slab=slab, silent=False)
    except ParameterError:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmp_manifest)
        return 2
    # set logging level
    if override_log_level is not None:
        rp.LOG_LEVEL = override_log_level
        logger.info("Overriding log level to {str(override_log_level)}.")
    logger.setLevel(rp.LOG_LEVEL)
    logger.debug("PARAMETERS file was read successfully")
    rp.timestamp = timestamp
    rp.manifest = tmp_manifest
    for p in preset_params:
        try:
            setattr(rp, p, preset_params[p])
        except Exception:
            logger.warning(f"Error applying preset parameter {p}: ",
                           exc_info=True)
    if not domains:
        slab.fullUpdate(rp)   # gets PARAMETERS data into slab
        rp.fileLoaded["POSCAR"] = True
    # set source directory
    _source = Path(source).resolve()
    if not _source.is_dir():
        logger.warning(f"tensorleed directory {source} not found.")
    if _source.name == "tensorleed":
        rp.source_dir = _source
    elif _source.parent.name == "tensorleed":
        logger.warning(f"tensorleed directory found in {_source.parent}, "
                       "using that instead of {_source}.")
        rp.source_dir = _source.parent
    elif (_source / "tensorleed").is_dir():
        logger.warning(f"tensorleed directory found in {_source}, using that "
                       "instead of {_source}.")
        rp.source_dir = _source / "tensorleed"
    else:
        logger.warning(f"Could not find a tensorleed directory at {_source}. "
                       "This may cause errors.")
        rp.source_dir = _source
    if system_name is not None:
        rp.systemName = system_name
    else:
        logger.info('No system name specified. Using name "unknown".')
        rp.systemName = "unknown"
    # check if halting condition is already in effect:
    if rp.halt >= rp.HALTING:
        logger.info("Halting execution...")
        cleanup(rp.manifest, rp)
        return 0
    rp.updateDerivedParams()
    logger.info(f"ViPErLEED is using TensErLEED version {rp.TL_VERSION_STR}.")
    prerun_clean(rp, log_name)
    exit_code = section_loop(rp, slab)
    # Finalize logging - if not done, will break unit testing
    logger.handlers.clear()
    logging.shutdown()
    return exit_code
