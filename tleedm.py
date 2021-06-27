# -*- coding: utf-8 -*-
"""
Created on Nov 12 2019

@author: Florian Kraushofer

Master script running the TensErLEED Manager.
"""

import sys
import time
import logging
import os
import shutil
import multiprocessing

cd = os.path.realpath(os.path.dirname(__file__))
# NB: it's necessary to add vpr_path to sys.path so that viperleed
#     can be loaded correctly at the top-level package
vpr_path = os.path.realpath(os.path.join(cd, '..'))
for import_path in (cd, vpr_path):
    if import_path not in sys.path:
        sys.path.append(import_path)

from viperleed.tleedmlib import Rparams
from viperleed.tleedmlib.base import CustomLogFormatter
from viperleed.tleedmlib.sections.run_sections import section_loop
from viperleed.tleedmlib.sections.cleanup import prerun_clean, cleanup
from viperleed import GLOBALS
from viperleed.tleedmlib.files.parameters import (readPARAMETERS,
                                                  interpretPARAMETERS)
from viperleed.tleedmlib.files.poscar import readPOSCAR

logger = logging.getLogger("tleedm")


def run_tleedm(system_name="", console_output=True, slab=None,
               preset_params={}, source="."):
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
    logname = 'tleedm-'+timestamp+'.log'
    logger = logging.getLogger("tleedm")
    logger.setLevel(logging.INFO)
    logFormatter = CustomLogFormatter()
    fileHandler = logging.FileHandler(logname, mode="w")
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    if console_output:
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        logger.addHandler(consoleHandler)

    logger.info("Starting new log: " + logname + "\nTime of execution (UTC): "
                + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    logger.info("This is ViPErLEED version " + GLOBALS["version"] + "\n")
    logger.info("! THIS VERSION IS A PRE-RELEASE NOT MEANT FOR PUBLIC "
                "DISTRIBUTION !\n")

    tmpmanifest = ["SUPP", "OUT", logname]
    try:
        rp = readPARAMETERS()
    except FileNotFoundError:
        if not preset_params:
            logger.error("No PARAMETERS file found, and no preset parameters "
                         "passed. Execution will stop.")
            cleanup(tmpmanifest)
            return 2
        rp = Rparams()
    except Exception:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmpmanifest)
        return 2

    # check if this is going to be a domain search
    domains = False
    if "DOMAIN" in rp.readParams:
        domains = True

    if domains:  # no POSCAR in main folder for domain searches
        slab = None
    elif slab is None:
        if os.path.isfile(os.path.join(".", "POSCAR")):
            poscarfile = os.path.join(".", "POSCAR")
            logger.info("Reading structure from file POSCAR")
            try:
                slab = readPOSCAR(filename=poscarfile)
            except Exception:
                logger.error("Exception while reading POSCAR", exc_info=True)
                cleanup(tmpmanifest)
                return 2
        else:
            logger.error("POSCAR not found. Stopping execution...")
            cleanup(tmpmanifest)
            return 2

        if not slab.preprocessed:
            logger.info("The POSCAR file will be processed and overwritten. "
                        "Copying the original POSCAR to POSCAR_user...")
            try:
                shutil.copy2(poscarfile, "POSCAR_user")
                tmpmanifest.append("POSCAR_user")
            except Exception:
                logger.error("Failed to copy POSCAR to POSCAR_user. Stopping "
                             "execution...")
                cleanup(tmpmanifest)
                return 2
    try:
        interpretPARAMETERS(rp, slab=slab)
        if rp.LOG_DEBUG:
            logger.setLevel(logging.DEBUG)
            logger.debug("PARAMETERS file was read successfully")
    except Exception:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmpmanifest)
        return 2
    rp.timestamp = timestamp
    rp.manifest = tmpmanifest
    for p in preset_params:
        try:
            setattr(rp, p, preset_params[p])
        except Exception:
            logger.warning("Error applying preset parameter {}: ".format(p),
                           exc_info=True)
    if not domains:
        slab.fullUpdate(rp)   # gets PARAMETERS data into slab
        rp.fileLoaded["POSCAR"] = True
        if slab.preprocessed:
            rp.SYMMETRY_FIND_ORI = False

    rp.systemName = system_name
    rp.sourcedir = source
    if not rp.systemName:
        # use name of parent folder
        rp.systemName = os.path.basename(os.path.abspath(
            os.path.join(os.getcwd(), os.pardir)))
    # check if halting condition is already in effect:
    if rp.halt >= rp.HALTING:
        logger.info("Halting execution...")
        cleanup(rp.manifest, rp)
        return 0

    rp.updateDerivedParams()

    prerun_clean(rp, logname)

    return section_loop(rp, slab)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    run_tleedm()
