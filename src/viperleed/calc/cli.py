"""Module cli of viperleed.calc.

Contains the definition of the main interface used for running
viperleed.calc from the command-line as in
>>> python3 -m viperleed.calc
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import multiprocessing
from pathlib import Path
import shutil

from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.log import LOGGER as bookie_log
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.constants import DEFAULT_DELTAS
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_TENSORS
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.constants import LOG_VERBOSE
from viperleed.calc.constants import LOG_VERY_VERBOSE
from viperleed.calc.files.manifest import ManifestFile
from viperleed.calc.files.manifest import ManifestFileError
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.log_utils import close_all_handlers
from viperleed.calc.lib.fs_utils import copytree_exists_ok
from viperleed.calc.lib.leedbase import getMaxTensorIndex
from viperleed.calc.run import run_calc
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.cli_base import ViPErLEEDCLI


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

        cwd_was_empty = not _has_history()

        bookkeeper = Bookkeeper()
        bookkeeper.run(mode=BookkeeperMode.CLEAR,
                       requires_user_confirmation=not args.skip_confirmation)

        _copy_tensors_and_deltas_to_work(work_path, args.all_tensors)           # TODO: it would be nice if all_tensors automatically checked PARAMETERS
        _copy_input_files_to_work(work_path)
        cwd_was_empty &= not any(work_path.iterdir())

        cwd = Path.cwd()
        exit_code = 2
        with execute_in_dir(work_path):
            try:
                exit_code, _ = run_calc(
                    system_name=args.name,
                    source=args.tensorleed,
                    preset_params=presets,
                    home=cwd,
                    )
            finally:
                domains = _copy_files_from_manifest(cwd)

        if cwd_was_empty:
            # Remove the empty history and history.info created by
            # the first bookkeeper call.
            _remove_history()
        else:
            # Run bookkeeper in archive mode,
            # propagating to domains if needed
            bookkeeper.run(
                mode=BookkeeperMode.ARCHIVE,
                requires_user_confirmation=not args.skip_confirmation,
                domains=domains,
                )

        # Finally clean up work if requested
        keep_workdir = args.keep_workdir or (exit_code and not cwd_was_empty)
        if not keep_workdir:
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
            action='store_true',
            )
        verbosity.add_argument(
            '-vv', '--very-verbose',
            help='increase output verbosity and print more debug messages',
            action='store_true',
            )

        # PATHS
        parser.add_argument(
            '-w', '--work',
            help='specify execution work directory',
            type=str,
            )
        parser.add_argument(
            '--tensorleed', '-t',
            help=('specify the path to the folder containing '
                  'the TensErLEED and EEASISSS source codes'),
            type=str,
            )

        # CREATING/DELETING DIRECTORIES
        parser.add_argument(
            '--all-tensors',
            help=(f'copy all {DEFAULT_TENSORS} to the work directory. '
                  'Required if using the TENSOR_INDEX parameter to calculate '
                  'from old tensors'),
            action='store_true',
            )
        parser.add_argument(
            '--keep-workdir', '-k',
            help=('do not delete the work directory after execution. By '
                  'default, the work directory is also not deleted in '
                  'case of errors'),
            action='store_true',
            )

        # OTHERS
        parser.add_argument(
            '-y',
            help=('automatically reply "yes" to all '
                  'requests for user confirmation'),
            action='store_true',
            dest='skip_confirmation',
            )


def _copy_input_files_to_work(work_path):
    """Copy all the known input files present here into `work_path`."""
    for file in ALL_INPUT_FILES:
        try:
            shutil.copy2(file, work_path)
        except FileNotFoundError:
            pass


def _copy_files_from_manifest(to_path):
    """Copy all files listed in file 'manifest' back `to_path`.

    Parameters
    ----------
    to_path : Path
        Where files/folders listed in the manifest file of
        the current directory should be copied to.

    Returns
    -------
    domain_paths : list
        Absolute paths to the subfolders of `to_path` that contain
        the results for subdomains in a DOMAIN calculation.
    """
    manifest = ManifestFile()
    manifest.read()
    if manifest.has_absolute_paths:
        raise ManifestFileError('Cannot copy resources from folders that are '
                                f'not contained in {Path.cwd()}. Destination '
                                'is not well defined.')

    dest_folders = []
    for folder, contents in manifest.iter_sections(relative=True):
        dest_folder = to_path/folder
        dest_folder.mkdir(exist_ok=True)
        for item in contents:
            src = folder/item
            _copy = shutil.copy2 if src.is_file() else copytree_exists_ok
            try:
                _copy(src, dest_folder / item)
            except OSError as exc:
                print(f'Error copying {src} to home directory: {exc}')          # TODO: Why no logging?
        dest_folders.append(dest_folder)
    return dest_folders[1:]  # The first one is always to_path


def _copy_tensors_and_deltas_to_work(work_path, all_tensors):
    """Move appropriate files from 'Tensors' and 'Deltas' to `work_path`."""
    if all_tensors:  # Copy all of them
        for directory in (DEFAULT_TENSORS, DEFAULT_DELTAS):
            try:
                copytree_exists_ok(directory, work_path / directory)
            except FileNotFoundError:
                pass
        return

    # Only copy the most recent Tensors/Deltas from zip files
    tensor_num = getMaxTensorIndex(zip_only=True)
    if not tensor_num:
        return

    for local_dir in (DEFAULT_TENSORS, DEFAULT_DELTAS):
        local_file = Path(f'{local_dir}/{local_dir}_{tensor_num:03d}.zip')
        if not local_file.is_file():
            continue
        work_directory = work_path / local_dir
        work_directory.mkdir(parents=True, exist_ok=True)
        shutil.copy2(local_file, work_path / local_file)


def _has_history():
    """Return whether the current directory was processed by bookkeeper."""
    cwd = Path.cwd()
    return (
        (cwd/DEFAULT_HISTORY).is_dir()
        or
        (cwd/HISTORY_INFO_NAME).is_file()
        )


def _make_work_directory(cli_args):
    """Return a suitable 'work' directory from `cli_args`."""
    work_path = Path(cli_args.work or DEFAULT_WORK).resolve()
    work_path.mkdir(parents=True, exist_ok=True)
    # Resolve again, in case it did not exist yet
    return work_path.resolve()


def _remove_history():
    """Delete history[.info] from the current directory."""
    close_all_handlers(bookie_log)
    cwd = Path.cwd()
    shutil.rmtree(cwd/DEFAULT_HISTORY)
    (cwd/HISTORY_INFO_NAME).unlink()


def _verbosity_to_log_level(cli_args, presets):
    """Add a LOG_LEVEL to `presets` if `cli_args` have verbosity specified."""
    if cli_args.very_verbose:
        presets['LOG_LEVEL'] = LOG_VERY_VERBOSE
    elif cli_args.verbose:
        presets['LOG_LEVEL'] = LOG_VERBOSE
