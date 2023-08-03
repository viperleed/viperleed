#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Nov 12 2019

@author: Florian Kraushofer
@author: Alexander Imre

Master script running the TensErLEED Manager.
"""

from pathlib import Path
import argparse
import logging
import multiprocessing
import os
import shutil
import sys
import time

import viperleed
from viperleed import GLOBALS
from viperleed.bookkeeper import bookkeeper, BookkeeperMode
from viperleed.lib.base import CustomLogFormatter
from viperleed.lib.classes import rparams
from viperleed.lib.files.parameter_errors import ParameterError
from viperleed.lib.files.parameters import (readPARAMETERS,
                                                  interpretPARAMETERS)
from viperleed.lib.files.poscar import readPOSCAR
from viperleed.lib.leedbase import getMaxTensorIndex
from viperleed.lib.sections._sections import ALL_INPUT_FILES
from viperleed.lib.sections.cleanup import prerun_clean, cleanup
from viperleed.lib.sections.run_sections import section_loop

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


def _parse_command_line_arguments():
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-w", "--work",
        help=("specify execution work directory"),
        type=str
        )
    parser.add_argument(
        "--all-tensors",
        help=("Copy all Tensors to the work directory. Required if using "
              "the TENSORS parameter to calculate from old tensors."),
        action='store_true'
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
        "--name", "-n",
        help = "specify system name",
        type=str
    )
    parser.add_argument(
        "--no_cont",
        help="Do not overwrite POSCAR with the new structure after a search.",
        action='store_true'
    )
    parser.add_argument(
        "--version",
        help=("print version information and exit"),
    )
    parser.add_argument(
        "--tensorleed", "-t",
        help="specify path to TensErLEED source code",
        type=str
        )
    parser.add_argument(                                                        #TODO: implement (for cont at end; warn if called with --no_cont)
        "-j", "--job_name",
        help=("defines a name for the current run. Will be appended to the name "
              "of the history folder that is created, and is logged in "
              "history.info. Passed along to the bookkeeper."),
        type=str)
    parser.add_argument(                                                        # TODO: implement
        "--history_name",
        help=("defines the name of the history folder that is created/used. "
              "Passed along to the bookkeeper. Default is 'history'."),
        type=str,
        default="history")
    parser.add_argument(                                                        # TODO: implement
        "--work_history_name",
        help=("defines the name of the workhistory folder that is created/used. "
              "Passed along to the bookkeeper. Default is 'workhistory'."),
        type=str,
        default="workhistory")
    args, _ = parser.parse_known_args()
    return args


def _interpret_tensorleed_path_flag(args):
    # if tensorleed arg is given, use that
    if args.tensorleed:
        tensorleed_path = Path(args.tensorleed)
    # else check environment variable $VIPERLEED_TENSORLEED
    try:
        tensorleed_path = Path(os.environ["VIPERLEED_TENSORLEED"])
    except KeyError:
        # environment variable not set
        raise RuntimeError(
            "TensErLEED path not specified.\n"
            "Please either pass a path to the TensErLEED source code with "
            "the --tensorleed argument, or set the environment variable "
            "$VIPERLEED_TENSORLEED."
        )
    return tensorleed_path


def main():
    multiprocessing.freeze_support() # needed for Windows

    args = _parse_command_line_arguments()

    if args.version:
        print(f"ViPErLEED version {GLOBALS['version']}")
        return 0

    if args.work:
        work_path = Path(args.work)
    else:
        work_path = Path.cwd() / "work"
    work_path = work_path.resolve()

    # tensorleed source directory
    _tensorleed_path = _interpret_tensorleed_path_flag(args)

    # verbosity flags
    if sum([args.verbose, args.very_verbose]) > 1:
        # only one verbosity level can be chosen
        logger.error("Only one verbosity level can be chosen. Stopping ")
        return 2
    elif args.very_verbose:
        override_log_level = 1
    elif args.verbose:
        override_log_level = 5
    else:
        override_log_level = None

    # system name flag
    if args.name:
        _system_name = args.name
    else:
        # use name of parent directory
        _system_name = Path.cwd().resolve().parent.name
        logger.info("No system name specified. Using name of parent directory: "
                    f"{_system_name}")


    print("Running bookkeeper...")
    # NB: job_name is None, because this is cleanup for the previous run
    bookkeeper(mode=BookkeeperMode.DEFAULT,
               job_name=None,
               history_name=args.history_name,
               work_history_name=args.work_history_name)


    # create work directory if necessary
    os.makedirs(work_path, exist_ok=True)

    all_tensors = args.all_tensors                                              # TODO: is there any need for this?
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
            tensor_file = os.path.join("Tensors",
                                       f"Tensors_{tensor_num:03d}.zip")
            shutil.copy2(tensor_file, os.path.join(work_path, tensor_file))
            delta_file = os.path.join("Deltas", f"Deltas_{tensor_num:03d}.zip")
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
    cwd = Path.cwd().resolve()
    os.chdir(work_path)
    run_tleedm(
        system_name=_system_name,
        source=_tensorleed_path,
        override_log_level=override_log_level
    )

    # copy back everything listed in manifest
    manifest = []
    if os.path.isfile("manifest"):
        with open("manifest", "r") as rf:
            manifest = [s.strip() for s in rf.readlines()]
    for _path in manifest:
        try:
            if os.path.isfile(_path):
                shutil.copy2(_path, cwd / _path)
            elif os.path.isdir(_path):
                shutil.copytree(_path, cwd / _path, dirs_exist_ok=True)
        except Exception as exception:
            print(f"Error copying {_path} to home directory: {str(exception)}")

    # go back, clean up if requested
    os.chdir(cwd)

    # call bookkeeper again to clean up unless --no_cont is set
    if not args.no_cont:
        bookkeeper(
            mode=BookkeeperMode.CONT,
            job_name=args.job_name,
            history_name=args.history_name,
            work_history_name=args.work_history_name
        )

    # delete work directory if requested
    if args.delete_workdir:
        try:
            shutil.rmtree(work_path)
        except OSError as error:
            print(f"Error deleting work directory: {str(error)}")


if __name__ == "__main__":
    main()
