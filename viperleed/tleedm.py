# -*- coding: utf-8 -*-
"""
Created on Nov 12 2019

@author: Florian Kraushofer
@author: Alexander Imre

Master script running the TensErLEED Manager.
"""

import argparse
import logging
import multiprocessing
import os
from pathlib import Path
import shutil
import sys
import time

# NB: it's necessary to add vpr_path to sys.path so that viperleed
#     can be loaded correctly at the top-level package
cd = Path(__file__).resolve().parent
vpr_path = cd.parent
for import_path in (str(cd), str(vpr_path)):
    if import_path not in sys.path:
        sys.path.append(import_path)

import viperleed
from viperleed import GLOBALS
from viperleed.tleedmlib.base import CustomLogFormatter
from viperleed.tleedmlib.classes import rparams
from viperleed.tleedmlib.files.parameters import (readPARAMETERS,
                                                  interpretPARAMETERS)
from viperleed.tleedmlib.files.parameter_errors import ParameterError
from viperleed.tleedmlib.files.poscar import readPOSCAR
from viperleed.tleedmlib.leedbase import getMaxTensorIndex
from viperleed.tleedmlib.sections.run_sections import section_loop
from viperleed.tleedmlib.sections.cleanup import prerun_clean, cleanup
from viperleed.tleedmlib.sections._sections import ALL_INPUT_FILES
from viperleed.bookkeeper import bookkeeper

logger = logging.getLogger("tleedm")


def run_tleedm(system_name="", console_output=True, slab=None,
               preset_params={}, source=".",
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
        rp = rparams.Rparams()
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
        poscarfile = Path("POSCAR")
        if poscarfile.is_file():
            logger.info("Reading structure from file POSCAR")
            try:
                slab = readPOSCAR(filename=str(poscarfile.resolve()))
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
        # interpret the PARAMETERS file
        interpretPARAMETERS(rp, slab=slab, silent=False)
    except ParameterError:
        logger.error("Exception while reading PARAMETERS file", exc_info=True)
        cleanup(tmpmanifest)
        return 2

    # set logging level
    if override_log_level is not None:
        rp.LOG_LEVEL = override_log_level
        logger.info("Overriding log level to {str(override_log_level)}.")
    logger.setLevel(rp.LOG_LEVEL)
    logger.debug("PARAMETERS file was read successfully")

    rp.timestamp = timestamp
    rp.manifest = tmpmanifest
    for p in preset_params:
        try:
            setattr(rp, p, preset_params[p])
        except Exception:
            logger.warning(f"Error applying preset parameter {p}: ",
                           exc_info=True)
    if not domains:
        slab.fullUpdate(rp)   # gets PARAMETERS data into slab
        rp.fileLoaded["POSCAR"] = True

    rp.systemName = system_name
    rp.sourcedir = str(Path(source).resolve())                                  # TODO: use Path instead of str
    if not rp.systemName:
        # use name of parent folder
        rp.systemName = str(Path.cwd().parent.name)
    # check if halting condition is already in effect:
    if rp.halt >= rp.HALTING:
        logger.info("Halting execution...")
        cleanup(rp.manifest, rp)
        return 0

    rp.updateDerivedParams()
    logger.info(f"ViPErLEED is using TensErLEED version {rp.TL_VERSION_STR}.")

    prerun_clean(rp, logname)
    exit_code = section_loop(rp, slab)

    # Finalize logging - if not done, will break unit testing
    logger.handlers.clear()
    logging.shutdown()

    return exit_code


if __name__ == "__main__":
    multiprocessing.freeze_support() # needed for Windows

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-w", "--work",
        help=("specify execution work directory"),
        type=str
        )
    parser.add_argument(
        "--delete_workdir",
        help=("delete work directory after execution"),
        action='store_true'
        )
    parser.add_argument(
        "-v", "--verbose",
        help=("increase output verbosity and print debug messages"),
        action='store_true'
        )
    parser.add_argument(
        "-vv", "--very_verbose",
        help=("increase output verbosity and prints more debug messages"),
        action='store_true'
    )
    parser.add_argument(
        "-vvv", "--very_very_verbose",
        help=("maximum output verbosity and even more debug messages"),
        action='store_true'
    )
    parser.add_argument(
        "--version",
        help=("print version information and exit"),
    )
    args, bookie_args = parser.parse_known_args()
    sys.argv = sys.argv[:1] + bookie_args

    # run bookkeeper # TODO: make this optional
    print("Running bookkeeper...")
    bookkeeper()

    if args.work:
        work_path = Path(args.work)
    else:
        work_path = Path.cwd() / "work"
    delete_workdir = args.delete_workdir

    if sum([args.verbose, args.very_verbose, args.very_very_verbose]) > 1:
        # only one verbosity level can be chosen
        logger.error("Only one verbosity level can be chosen. Stopping ")
        sys.exit(2)
    elif args.very_very_verbose:
        override_log_level = 1
    elif args.very_verbose:
        override_log_level = 10
    elif args.verbose:
        override_log_level = logging.DEBUG
    else:
        override_log_level = None

    # create work directory if necessary
    os.makedirs(work_path, exist_ok=True)

    all_tensors = False                                                         # TODO: is there any need for this?
    # !!! TODO: it would be nice if all_tensors automatically checked PARAMETERS

    # copy Tensors and Deltas to work directory
    if all_tensors:
        try:
            shutil.copytree("Tensors", os.path.join(work_path, "Tensors"),
                            dirs_exist_ok=True)
        except FileNotFoundError:
            pass
        try:
            shutil.copytree("Deltas", os.path.join(work_path, "Deltas"),
                            dirs_exist_ok=True)
        except FileNotFoundError:
            pass
    else:
        tensor_num = getMaxTensorIndex(zip_only=True)
        if tensor_num > 0:
            os.makedirs(os.path.join(work_path, "Tensors"), exist_ok=True)
            tensor_file = os.path.join("Tensors", "Tensors_{:03d}.zip"
                                      .format(tensor_num))
            shutil.copy2(tensor_file, os.path.join(work_path, tensor_file))
            delta_file = os.path.join("Deltas", "Deltas_{:03d}.zip"
                                     .format(tensor_num))
            if os.path.isfile(delta_file):
                os.makedirs(os.path.join(work_path, "Deltas"), exist_ok=True)
                shutil.copy2(delta_file, os.path.join(work_path, delta_file))

    # copy input files to work directory
    for file in ALL_INPUT_FILES:
        try:
            shutil.copy2(file, work_path / file)
        except FileNotFoundError:
            pass

    # go to work directory, execute there
    cwd = Path.cwd()
    os.chdir(work_path)
    run_tleedm(source=os.path.join(vpr_path, "viperleed"),
               override_log_level=override_log_level)

    # copy back everything listed in manifest
    manifest = []
    if os.path.isfile("manifest"):
        with open("manifest", "r") as rf:
            manifest = [s.strip() for s in rf.readlines()]
    for p in manifest:
        try:
            if os.path.isfile(p):
                shutil.copy2(p, cwd / p)
            elif os.path.isdir(p):
                shutil.copytree(p, cwd / p, dirs_exist_ok=True)
        except Exception as e:
            print(f"Error copying {p} to home directory: {str(e)}")

    # go back, clean up if requested
    os.chdir(cwd)
    if delete_workdir:
        try:
            shutil.rmtree(work_path)
        except Exception as e:
            print("Error deleting work directory: " + str(e))
