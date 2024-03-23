"""Module cli of viperleed.calc.

Contains the definition of the main interface used for running
viperleed.calc from the command-line as in
>>> python3 -m viperleed.calc
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import argparse
import multiprocessing
import os
from pathlib import Path
import shutil

from viperleed import GLOBALS
from viperleed.calc import LOGGER as logger
from viperleed.calc.bookkeeper import bookkeeper, BookkeeperMode
from viperleed.calc.lib.base import copytree_exists_ok
from viperleed.calc.lib.leedbase import getMaxTensorIndex
from viperleed.calc.sections._sections import ALL_INPUT_FILES
from viperleed.calc.run import get_tensorleed_path                              # TODO: let run_calc do this!
from viperleed.calc.run import run_calc


DEFAULT_HISTORY = 'history'
DEFAULT_WORK = 'work'
DEFAULT_WORK_HISTORY = 'workhistory'
LOG_VERBOSE = 5
LOG_VERY_VERBOSE = 1


def add_calc_parser_arguments(parser):
    """Add CLI arguments for viiperleed.calc to parser."""
    parser.add_argument(
        '--version',
        action='version',
        version=GLOBALS['version_message'],
        )

    parser.add_argument(
        '--name', '-n',
        help=('Specify name of the system/sample for which '
              'calculations are run (e.g., Ag(100)-p1x1)'),
        type=str
        )

    # VERBOSITY
    parser.add_argument(
        '-v', '--verbose',
        help='Increase output verbosity and print debug messages',
        action='store_true'
        )
    parser.add_argument(
        '-vv', '--very-verbose',
        help='Increase output verbosity and print more debug messages',
        action='store_true'
        )

    # PATHS
    parser.add_argument(
        '-w', '--work',
        help='Specify execution work directory',
        type=str
        )
    parser.add_argument(
        '--tensorleed', '-t',
        help=('Specify the path to the folder containing '
              'the TensErLEED and EEASISSS source codes'),
        type=str
        )

    # BOOKKEPER
    parser.add_argument(
        '--no-cont',
        help='Do not overwrite POSCAR/VIBROCC with those after a search',
        action='store_true'
        )
    parser.add_argument(                                                        # TODO: implement (for cont at end; warn if called with --no_cont)
        '-j', '--job-name',
        help=('Define a name for the current run. Will be appended to the '
              'name of the history folder that is created, and is logged '
              'in history.info. Passed along to the bookkeeper'),
        type=str
        )
    parser.add_argument(
        '--history-name',
        help=('Define the name of the history folder that is '
              'created/used. Passed along to the bookkeeper. '
              f'Default is {DEFAULT_HISTORY!r}'),
        type=str,
        default=DEFAULT_HISTORY
        )
    parser.add_argument(
        '--work-history-name',
       help=('Define the name of the workhistory folder that '
             'is created/used. Passed along to the bookkeeper. '
             f'Default is {DEFAULT_WORK_HISTORY!r}'),
        type=str,
        default=DEFAULT_WORK_HISTORY
        )

    # CREATING/DELETING DIRECTORIES
    parser.add_argument(
        '--all-tensors',
        help=('Copy all Tensors to the work directory. Required if using '
              'the TENSORS parameter to calculate from old tensors'),
        action='store_true'
        )
    parser.add_argument(
        '--delete-workdir',
        help='Delete work directory after execution',
        action='store_true'
        )


def main(args=None):
    """Read command-line arguments and execute viperleed.calc."""
    # Make sure that pre-packed 'executables'
    # work with multiprocessing on Windows
    multiprocessing.freeze_support()

    args = args or _parse_cli_args()
    work_path = _make_work_directory(args)

    presets = {}  # Replace selected PARAMETERS
    try:
        _verbosity_to_log_level(args, presets)
    except ValueError as exc:
        logger.error(f'{exc} Stopping.')
        return 2

    print('Running bookkeeper...')                                              # TODO: This is lost to stdout if we don't log it
    # NB: job_name is None, because we're cleaning up the previous run
    bookkeeper(mode=BookkeeperMode.DEFAULT,
               job_name=None,
               history_name=args.history_name,
               work_history_name=args.work_history_name)

    _copy_tensors_and_deltas_to_work(work_path, args.all_tensors)               # TODO: it would be nice if all_tensors automatically checked PARAMETERS
    _copy_input_files_to_work(work_path)

    # Go to work directory, execute there
    cwd = Path.cwd().resolve()
    os.chdir(work_path)
    try:
        exit_code = run_calc(
            source=get_tensorleed_path(args.tensorleed),
            system_name=args.name,
            preset_params=presets
            )
    finally:
        # Copy back everything listed in manifest, then go back
        _copy_files_from_manifest(cwd)
        os.chdir(cwd)

    # Call bookkeeper again to clean up unless --no_cont is set
    if not args.no_cont:
        bookkeeper(mode=BookkeeperMode.CONT,
                   job_name=args.job_name,
                   history_name=args.history_name,
                   work_history_name=args.work_history_name)

    # Finally clean up work if requested
    if args.delete_workdir:
        try:
            shutil.rmtree(work_path)
        except OSError as exc:
            print(f'Error deleting work directory: {exc}')                      # TODO: This is lost to stdout if we don't log it

    return exit_code


def _copy_input_files_to_work(work_path):
    """Copy all the known input files present here into work_path."""
    for file in ALL_INPUT_FILES:
        try:
            shutil.copy2(file, work_path / file)
        except FileNotFoundError:
            pass


def _copy_files_from_manifest(to_path):
    """Copy all files listed in file 'manifest' back to_path."""
    manifest_file = Path('manifest')
    if not manifest_file.is_file():
        return

    with manifest_file.open('r', encoding='utf-8') as file:
        manifest = [line.strip() for line in file.readlines()]
        manifest = (line for line in manifest if line)
        manifest = (Path(line) for line in manifest)

    for to_be_copied in manifest:
        _copy = shutil.copy2 if to_be_copied.is_file() else copytree_exists_ok
        try:
            _copy(to_be_copied, to_path / to_be_copied)
        except OSError as exc:
            print(f'Error copying {to_be_copied} to home directory: {exc}')     # TODO: Why no logging?


def _copy_tensors_and_deltas_to_work(work_path, all_tensors):
    """Move appropriate files from 'Tensors' and 'Deltas' to work_path."""
    if all_tensors:  # Copy all of them
        for directory in ('Tensors', 'Deltas'):
            try:
                shutil.copytree(directory, work_path / directory,
                                dirs_exist_ok=True)
            except FileNotFoundError:
                pass
        return

    # Only copy the most recent Tensors/Deltas from zip files
    tensor_num = getMaxTensorIndex(zip_only=True)
    if not tensor_num:
        return

    for local_dir in ('Tensors', 'Deltas'):
        local_file = Path(f'{local_dir}/{local_dir}_{tensor_num:03d}.zip')
        if not local_file.is_file():
            continue
        work_directory = work_path / local_dir
        work_directory.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local_file, work_path / local_file)


def _make_work_directory(cli_args):
    """Return a suitable 'work' directory from cli_args."""
    work_path = Path(cli_args.work or DEFAULT_WORK).resolve()
    work_path.mkdir(parents=True, exist_ok=True)
    return work_path


def _parse_cli_args():
    """Return parsed arguments from args."""
    parser = argparse.ArgumentParser()
    add_calc_parser_arguments(parser)
    return parser.parse_args()


def _verbosity_to_log_level(cli_args, presets):
    """Add a LOG_LEVEL to presets if cli_args have verbosity specified."""
    if sum([cli_args.verbose, cli_args.very_verbose]) > 1:
        # only one verbosity level can be chosen
        raise ValueError('Only one verbosity level can be chosen.')
    if cli_args.very_verbose:
        presets['LOG_LEVEL'] = LOG_VERY_VERBOSE
    elif cli_args.verbose:
        presets['LOG_LEVEL'] = LOG_VERBOSE
