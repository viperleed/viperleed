
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""__main__.py of ViPErLEED (viperleed) calc

See viperleed.calc.__init__.py for more information.
"""
from pathlib import Path
import argparse
import logging
import multiprocessing
import os
import shutil

from viperleed.bookkeeper import bookkeeper, BookkeeperMode
from viperleed.calc import run_tleedm
from viperleed.calc.lib.leedbase import getMaxTensorIndex
from viperleed.calc.sections._sections import ALL_INPUT_FILES

__authors__ = ["Florian Kraushofer (@fkraushofer)",
               "Alexander M. Imre (@amimre)",
               "Michele Riva (@michele-riva)"]

logger = logging.getLogger("tleedm")


def add_calc_parser_arguments(parser):
    parser.add_argument(
        "-w", "--work",
        help=("specify execution work directory"),
        type=str
    )
    parser.add_argument(
        "--all-tensors",
        help=("Copy all Tensors to the work directory. Required if using "
            "the TENSORS parameter to calculate from old tensors."),
        action="store_true"
        )
    parser.add_argument(
        "--delete_workdir",
        help=("delete work directory after execution"),
        action="store_true"
        )
    parser.add_argument(
        "-v", "--verbose",
        help=("increase output verbosity and print debug messages"),
        action="store_true"
        )
    parser.add_argument(
        "-vv", "--very_verbose",
        help=("increase output verbosity and prints more debug messages"),
        action="store_true"
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
    parser.add_argument(
        "--history_name",
        help=("defines the name of the history folder that is created/used. "
            "Passed along to the bookkeeper. Default is 'history'."),
        type=str,
        default="history")
    parser.add_argument(
        "--work_history_name",
        help=("defines the name of the workhistory folder that is created/used. "
            "Passed along to the bookkeeper. Default is 'workhistory'."),
        type=str,
        default="workhistory")


def _interpret_tensorleed_path_flag(args):
    # if tensorleed arg is given, use that
    if args.tensorleed:
        tensorleed_path = Path(args.tensorleed)
    # else check environment variable $VIPERLEED_TENSORLEED
    try:
        tensorleed_path = Path(os.environ["VIPERLEED_TENSORLEED"])
    except KeyError as err:
        # environment variable not set
        raise RuntimeError(
            "TensErLEED path not specified.\n"
            "Please either pass a path to the TensErLEED source code with "
            "the --tensorleed argument, or set the environment variable "
            "$VIPERLEED_TENSORLEED."
        ) from err
    return tensorleed_path


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_calc_parser_arguments(parser)
        args = parser.parse_args()

    multiprocessing.freeze_support() # needed for Windows

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

    all_tensors = args.all_tensors
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


# run main() if this file is called directly
if __name__ == "__main__":
    main()
