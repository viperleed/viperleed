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

import multiprocessing
import os
from pathlib import Path
import shutil

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK
from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.lib.base import copytree_exists_ok
from viperleed.calc.lib.leedbase import getMaxTensorIndex
from viperleed.calc.run import run_calc
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.cli_base import ViPErLEEDCLI


LOG_VERBOSE = 5
LOG_VERY_VERBOSE = 1


class ViPErLEEDCalcCLI(ViPErLEEDCLI, cli_name='calc'):
    """Main command-line interface for viperleed.calc."""

    def __call__(self, args=None):
        """Call viperleed.calc."""
        # Make sure that pre-packed 'executables'
        # work with multiprocessing on Windows
        multiprocessing.freeze_support()

        args = self.parse_cli_args(args)
        work_path = _make_work_directory(args)

        presets = {}  # Replace selected PARAMETERS
        _verbosity_to_log_level(args, presets)

        # NB: job_name is None, as we're cleaning up the previous run
        bookkeeper = Bookkeeper(
            job_name=args.job_name,
            history_name=args.history_name,
            work_history_name=args.work_history_name,
            )
        bookkeeper.run(mode=BookkeeperMode.CLEAR)

        _copy_tensors_and_deltas_to_work(work_path, args.all_tensors)           # TODO: it would be nice if all_tensors automatically checked PARAMETERS
        _copy_input_files_to_work(work_path)

        # Go to work directory, execute there
        cwd = Path.cwd().resolve()
        os.chdir(work_path)
        exit_code = 2
        try:
            exit_code, _ = run_calc(
                system_name=args.name,
                source=args.tensorleed,
                preset_params=presets,
                )
        finally:
            # Copy back everything listed in manifest, then go back
            _copy_files_from_manifest(cwd)
            os.chdir(cwd)

        # update bookkeeper with new run info
        bookkeeper.update_from_cwd()
        # run bookkeeper in archive mode
        bookkeeper.run(mode=BookkeeperMode.ARCHIVE)

        # Finally clean up work if requested
        if args.delete_workdir:
            try:
                shutil.rmtree(work_path)
            except OSError as exc:
                print(f'Error deleting work directory: {exc}')                  # TODO: This is lost to stdout if we don't log it
        return exit_code

    def add_parser_arguments(self, parser):
        """Add CLI arguments for viperleed.calc to parser."""
        super().add_parser_arguments(parser)
        parser.add_argument(
            '--name', '-n',
            help=('specify name of the system/sample for which '
                  'calculations are run (e.g., Ag(100)-p1x1)'),
            type=str
            )

        # VERBOSITY
        verbosity = parser.add_mutually_exclusive_group()
        verbosity.add_argument(
            '-v', '--verbose',
            help='increase output verbosity and print debug messages',
            action='store_true'
            )
        verbosity.add_argument(
            '-vv', '--very-verbose',
            help='increase output verbosity and print more debug messages',
            action='store_true'
            )

        # PATHS
        parser.add_argument(                                                    # TODO: bookkeeper always assumes DEFAULT_WORK!
            '-w', '--work',
            help='specify execution work directory',
            type=str
            )
        parser.add_argument(
            '--tensorleed', '-t',
            help=('specify the path to the folder containing '
                  'the TensErLEED and EEASISSS source codes'),
            type=str
            )

        # BOOKKEPER
        parser.add_argument(                                                    # TODO: implement (for cont at end; warn if called with --no_cont)
            '-j', '--job-name',
            help=('define a name for the current run. Will be appended to the '
                  'name of the history folder that is created, and is logged '
                  'in history.info. Passed along to the bookkeeper'),
            type=str
            )
        parser.add_argument(
            '--history-name',
            help=('define the name of the history folder that is '
                  'created/used. Passed along to the bookkeeper. '
                  f'Default is {DEFAULT_HISTORY!r}'),
            type=str,
            default=DEFAULT_HISTORY
            )
        parser.add_argument(
            '--work-history-name',
           help=('define the name of the workhistory folder that '
                 'is created/used. Passed along to the bookkeeper. '
                 f'Default is {DEFAULT_WORK_HISTORY!r}'),
            type=str,
            default=DEFAULT_WORK_HISTORY
            )

        # CREATING/DELETING DIRECTORIES
        parser.add_argument(
            '--all-tensors',
            help=('copy all Tensors to the work directory. Required if using '
                  'the TENSORS parameter to calculate from old tensors'),
            action='store_true'
            )
        parser.add_argument(
            '--delete-workdir',
            help='delete work directory after execution',
            action='store_true'
            )


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
                copytree_exists_ok(directory, work_path / directory)
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


def _verbosity_to_log_level(cli_args, presets):
    """Add a LOG_LEVEL to presets if cli_args have verbosity specified."""
    if cli_args.very_verbose:
        presets['LOG_LEVEL'] = LOG_VERY_VERBOSE
    elif cli_args.verbose:
        presets['LOG_LEVEL'] = LOG_VERBOSE
