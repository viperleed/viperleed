"""ViPErLEED bookkeeper module of package calc."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from collections import defaultdict
from enum import Enum
from operator import attrgetter
import os
from pathlib import Path
import re
import shutil
import time

from viperleed.calc import DEFAULT_HISTORY
from viperleed.calc import DEFAULT_WORK
from viperleed.calc import DEFAULT_WORK_HISTORY
from viperleed.calc import LOG_PREFIX
from viperleed.calc import ORIGINAL_INPUTS_DIR_NAME
from viperleed.calc.lib import leedbase
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.sections.cleanup import PREVIOUS_LABEL
from viperleed.cli_base import ViPErLEEDCLI

_CALC_LOG_PREFIXES = (
    LOG_PREFIX,
    'tleedm',   # For backwards compatibility
    )
HIST_FOLDER_RE = re.compile(
    r't(?P<tensor_num>[0-9]{3}).r(?P<job_num>[0-9]{3})_'
    )
_HISTORY_INFO_SPACING = 12  # For the leftmost field in history.info


class BookkeeperMode(Enum):
    """Enumeration of bookkeeper modes.

    Attributes
    ----------
    CONT
        Store last run in history. Overwrite POSCAR & VIBROCC from OUT.
    DEFAULT
        Store last run in history. Do not overwrite POSCAR & VIBROCC.
    DISCARD
        Discard previous run as if it never happened.
    """
    CONT = 'cont'
    DEFAULT = 'default'
    DISCARD = 'discard'

    @property
    def discard(self):
        """Return whether this is mode DISCARD."""
        return self is BookkeeperMode.DISCARD


def store_input_files_to_history(root_path, history_path):
    """Find input files in `root_path` and copy them to `history_path`.

    Parameters
    ----------
    root_path : pathlike
        Root directory from which to take files. Should be cwd,
        not ./work.
    history_path : pathlike
        Path to the history directory in which the files should
        be stored.
    """
    root_path, history_path = Path(root_path), Path(history_path)
    original_inputs_path = root_path / DEFAULT_WORK / ORIGINAL_INPUTS_DIR_NAME  # TODO: we should use root_path / 'SUPP' / ORIGINAL_INPUTS_DIR_NAME instead! bookkeeper should not access the work directory (which may have a different name, or not be in root at all!)
    if original_inputs_path.is_dir():
        input_origin_path = original_inputs_path
    else:
        input_origin_path = root_path
        print(f'Could not find directory {ORIGINAL_INPUTS_DIR_NAME!r} with '
              'unaltered input files. Files will instead be copied from the '
              'root directory.')

    # Only files, no directories
    files_to_copy = (file for file in input_origin_path.iterdir()
                     if file.is_file() and file.name in ALL_INPUT_FILES)
    for file in files_to_copy:
        try:
            shutil.copy2(file, history_path/file.name)
        except OSError as exc:
            print(f'Failed to copy file {file} to history: {exc}')


# While there are indeed many arguments to this function, packing them
# up into a dedicated dataclass seems a bit overkill. Appropriately
# documenting them should be enough.
# pylint: disable-next=too-many-arguments
def _add_entry_to_history_info_file(cwd, tensor_nums, job_nums, job_name,
                                    last_log_lines, timestamp, dirname):
    """Add information about the current run to history.info.

    Parameters
    ----------
    cwd : Path
        The current directory, containing history.info
        and the optional notes files.
    tensor_nums : Sequence of int
        The progressive numbers of the tensors that were used
        or created during this run (and were moved to history).
    job_nums : Sequence of int
        The new job numbers for the current run (one per tensor).
    job_name : str
        Alternative name, typically given by the user as
        a command-line argument.
    last_log_lines : Sequence of str
        The lines read from the most recent viperleed.calc
        log file.
    timestamp : str
        The time-stamp for this run. Typically derived
        from the name of the viperleed.calc log file.
    dirname : str
        The name of the new history directory that was
        created during this bookkeeper run.

    Returns
    -------
    None.
    """
    history_info = cwd / 'history.info'
    contents = [] if not history_info.is_file() else ['', '']
    contents.append(
        '# TENSORS '.ljust(_HISTORY_INFO_SPACING)
        + ('None' if tensor_nums == {0} else str(tensor_nums)[1:-1])            # TODO: perhaps we should sort the tensor numbers?
        )
    contents.append('# JOB ID '.ljust(_HISTORY_INFO_SPACING)
                    + str(job_nums)[1:-1])                                      # TODO: perhaps we should sort the job numbers?
    if job_name is not None:
        contents.append(f'# JOB NAME {job_name}')
    contents.extend(_infer_run_info_from_log(last_log_lines))
    contents.append(
        '# TIME '.ljust(_HISTORY_INFO_SPACING)
        + f'{_translate_timestamp(timestamp)}'
        )
    contents.append('# FOLDER '.ljust(_HISTORY_INFO_SPACING) + f'{dirname}')
    contents.append(f'Notes: {_read_and_clear_notes_file(cwd)}')
    contents.append('\n###########')
    try:  # pylint: disable=too-many-try-statements
        with history_info.open('a', encoding='utf-8') as hist_info_file:
            hist_info_file.write('\n'.join(contents))
    except OSError:
        print('Error: Failed to append to history.info file.')


def _collect_log_files(cwd):
    """Return two lists of log files in `cwd`: 'calc' and others."""
    calc_logs, other_logs = [], []
    for file in cwd.glob('*.log'):
        if not file.is_file():
            continue
        container = (calc_logs if file.name.startswith(_CALC_LOG_PREFIXES)
                     else other_logs)
        container.append(file)
    return calc_logs, other_logs


def _collect_supp_and_out(cwd):
    """Return paths to 'SUPP' and 'OUT' if they are present in `cwd`."""
    return [d for d in cwd.iterdir()
            if d.is_dir() and d.name in {'OUT', 'SUPP'}]


def _create_new_history_directory(history_path, tensor_number,
                                  job_num, suffix, should_mkdir):
    """Return the path to a new subfolder of `history_path`.

    Parameters
    ----------
    history_path : Path
        Path to the history folder. The new
        subdirectory is created here.
    tensor_number : int
        Progressive index of the tensor for the new subdirectory.
    job_num : int
        Progressive number identifying the run for `tensor_number`
        for which the new directory is created.
    suffix : str
        Suffix to append to the directory name. Typically includes
        the time-stamp and, optionally, a user-defined job name.
    should_mkdir : bool
        Whether the new history folder should also be created on the
        filesystem as a new directory. If not given, only the path is
        returned, without any check for existence.

    Returns
    -------
    new_history_dir : Path
        Path to the new subdirectory of `history_path`.

    Raises
    ------
    OSError
        If `should_mkdir` and creation of `new_history_dir` fails.
    """
    dirname = f't{tensor_number:03d}.r{job_num:03d}_{suffix}'
    new_history_dir = history_path / dirname

    if not should_mkdir:
        return new_history_dir

    if new_history_dir.is_dir():
        suffix = suffix.replace('moved-', '')
        dirname_moved = f'{dirname}_moved-{suffix}'
        print(f'Error: Target directory {dirname} already '
              f'exists. Will use {dirname_moved} instead.')
        new_history_dir = history_path / dirname_moved
    try:
        new_history_dir.mkdir()
    except OSError:
        print('Error: Could not create target directory '
              f'{new_history_dir}\n Stopping...')
        raise
    return new_history_dir


def _discard_tensors_and_deltas(cwd, tensor_number):
    """Delete tensor and delta files with tensor_num in cwd."""
    tensor_file = cwd / 'Tensors' / f'Tensors_{tensor_number:03d}.zip'
    delta_file = cwd / 'Deltas' / f'Deltas_{tensor_number:03d}.zip'
    _move_or_discard_files((tensor_file, delta_file), cwd, True)


def _discard_workhistory_previous(work_history_path):
    """Remove 'previous'-labelled directories in `work_history_path`."""
    work_hist_prev = (
        d for d in work_history_path.glob(f'*{PREVIOUS_LABEL}*')
        if d.is_dir()
        and HIST_FOLDER_RE.match(d.name)
        )
    for directory in work_hist_prev:
        try:
            shutil.rmtree(directory)
        except OSError:
            print(f'Failed to delete {directory} directory '
                  f'from {work_history_path}')


def _find_max_run_per_tensor(history_path):
    """Return maximum run numbers for all directories in `history_path`.

    Parameters
    ----------
    history_path : Path
        The path to the history folder in which
        tensor-run directories are looked up.

    Returns
    -------
    max_nums : defaultdict(int)
        The maximum run number currently stored in `history_path`
        for each of the tensors in `history_path`.
    """
    max_nums = defaultdict(int)  # max. job number per tensor number
    for directory in history_path.iterdir():
        match = HIST_FOLDER_RE.match(directory.name)
        if not directory.is_dir() or not match:
            continue
        tensor_num = int(match['tensor_num'])
        job_num = int(match['job_num'])
        max_nums[tensor_num] = max(max_nums[tensor_num], job_num)
    return max_nums


def _infer_run_info_from_log(last_log_lines):
    """Yield history.info lines about executed segments and R-factors."""
    if not last_log_lines:
        return
    for i, line in enumerate(reversed(last_log_lines)):
        # try to read what was executed
        if line.startswith('Executed segments: '):
            run_info = line.split('Executed segments: ')[1].strip()
            yield '# RUN '.ljust(_HISTORY_INFO_SPACING) + run_info
            break
    else:
        return

    # Now try to read final R-factors
    _section_to_shorhand = {
        'refcalc': 'REF',
        'superpos': 'SUPER',
        }
    for j, line in enumerate(reversed(last_log_lines)):
        if not line.startswith('Final R'):
            continue
        if j >= i:  # R factors are after the segments
            return
        try:
            shorthand = next(
                short
                for section, short in _section_to_shorhand.items()
                if section in line
                )
        except StopIteration:
            continue
        yield (f'# R {shorthand} '.ljust(_HISTORY_INFO_SPACING)
               + line.split(':', maxsplit=1)[1].strip())


def _move_and_cleanup_workhistory(work_history_path,
                                  history_path,
                                  timestamp,
                                  max_job_for_tensor,
                                  discard):
    """Move files from `work_history_path` to `history_path`, then clean up.

    If `work_history_path` is empty, it is always deleted. Otherwise
    only if `discard` is True.

    Parameters
    ----------
    work_history_path : Path
        Path to the workhistory directory, moved in root by
        sections.cleanup, as it is included in manifest in
        move_oldruns.
    history_path : Path
        The target history folder where files should be copied.
    timestamp : str
        The time-stamp selecting the subfolders to be moved. Only
        those folders whose name contains timestamp are copied over.
    max_job_for_tensor : defaultdict(int)
        The maximum job number known for each tensor number. Used
        to generate new folder names (incrementing by one unit)
        in `history_path`.
    discard : bool
        Whether files should only be discarded without copying them.

    Returns
    -------
    tensor_nums : set of int
        Indices of tensors found in `work_history_path` that have
        been moved to `history_path` as new history entries.
    """
    tensor_nums = set()
    if not work_history_path.is_dir():
        return tensor_nums

    # Always remove any 'previous'-labelled folders
    _discard_workhistory_previous(work_history_path)
    if not discard:
        tensor_nums = _move_workhistory_folders(work_history_path,
                                                history_path,
                                                timestamp,
                                                max_job_for_tensor)
    is_empty = not any(work_history_path.iterdir())
    if is_empty or discard:
        try:
            shutil.rmtree(work_history_path)
        except OSError as exc:
            err_ = f'Failed to {{}} {work_history_path.name} folder: {exc}'
            print(err_.format('discard' if discard else 'delete empty'))
    return tensor_nums


def _move_or_discard_files(file_paths, target_folder, discard):
    """Move file_paths to target_folder or delete them."""
    for file in file_paths:
        _move_or_discard_one_file(file, target_folder, discard)


def _move_or_discard_one_file(file, target_folder, discard):
    """Move file to target_folder or delete it."""
    if not file.exists():
        return
    if discard and file.is_file():
        try:
            file.unlink()
        except OSError:
            print(f'Failed to discard file {file.name}.')
        return
    if discard:  # Should be a directory
        try:
            shutil.rmtree(file)
        except OSError:
            print(f'Failed to discard directory {file.name}.')
        return
    # Move it
    try:
        shutil.move(file, target_folder / file.name)
    except OSError:
        print(f'Error: Failed to move {file.name}.')


def _move_workhistory_folders(work_history_path, history_path,
                              timestamp, max_job_for_tensor):
    """Move timestamp-labelled folders in work_history_path to history_path.

    Parameters
    ----------
    work_history_path : Path
        Base path to the workhistory folder from
        which stuff should be moved.
    history_path : Path
        Destination path for the subfolders of
        `work_history_path` that will be moved.
    timestamp : str
        Selects which sub-folders of `work_history_path` should be
        moved. Only those whose name contains `timestamp` are moved.
    max_job_for_tensor : defaultdict(int)
        The maximum job number known for each tensor number. Used
        to generate new folder names (incrementing by one unit)
        in `history_path` for each of the directories moved from
        `work_history_path`.

    Returns
    -------
    tensor_nums : set
        The tensor numbers of the folders that have been moved.
    """
    tensor_nums = set()
    work_history_dirs = (
        d for d in work_history_path.glob(f'*{timestamp}*')
        if d.is_dir()
        and PREVIOUS_LABEL not in d.name
        )
    for directory in work_history_dirs:
        match = HIST_FOLDER_RE.match(directory.name)
        if not match:
            continue
        tensor_num = int(match['tensor_num'])
        try:
            search_num = int(directory.name[6:9])                               # TODO: is this right? 6:9 should be job_num
        except (ValueError, IndexError):
            continue

        tensor_nums.add(tensor_num)
        job_num = max_job_for_tensor[tensor_num] + 1
        newname = (
            f't{tensor_num:03d}.r{job_num:03d}.{search_num:03d}'
            + directory.name[9:]
            )
        try:
            shutil.move(directory, history_path / newname)
        except OSError:
            print(f'Error: Failed to move {directory}.')
    return tensor_nums


def _read_and_clear_notes_file(cwd):
    """Return notes read from file. Clear the file contents."""
    notes_path = next(cwd.glob('notes*'), None)
    if notes_path is None:
        return ''
    try:
        notes = notes_path.read_text(encoding='utf-8')
    except OSError:
        print(f'Error: Failed to read {notes_path.name} file.')
        return ''
    try:  # pylint: disable=too-many-try-statements
        with notes_path.open('w', encoding='utf-8'):
            pass
    except OSError:
        print(f'Error: Failed to clear the {notes_path.name} '
              'file after reading.')
    return notes


def _read_most_recent_log(cwd):
    """Return timestamp and lines from the most-recent log file in cwd."""
    last_log_lines = []
    try:
        most_recent_log = max(
            (log_file
             for prefix in _CALC_LOG_PREFIXES
             for log_file in cwd.glob(f'{prefix}*.log')
             if log_file.is_file()),
            key=attrgetter('name')
            )
    except ValueError:  # No log files
        timestamp = time.strftime('%y%m%d-%H%M%S', time.localtime())
        old_timestamp = f'moved-{timestamp}'
        return old_timestamp, last_log_lines

    old_timestamp = most_recent_log.name[-17:-4]
    try:  # pylint: disable=too-many-try-statements
        with most_recent_log.open('r', encoding='utf-8') as log_file:
            last_log_lines = log_file.readlines()
    except OSError:
        pass
    return old_timestamp, last_log_lines


def _replace_input_files_from_out(cwd):
    """Replace POSCAR and VIBROCC in cwd with those from OUT."""
    out_path = cwd / 'OUT'
    if not out_path.is_dir():
        return
    for file in ('POSCAR', 'VIBROCC'):
        out_files = (f for f in out_path.glob(f'{file}_OUT_*')
                     if f.is_file()
                     and 'parabola' not in f.name)
        try:
            most_recent = max(out_files, key=attrgetter('name'))                # TODO: is this right? are the OUT files sorted by timestamp or by R?
        except ValueError:  # No files
            # Do not complain if not found, since we move
            # previous runs to the history by default
            continue
        try:
            shutil.copy2(most_recent, cwd / file)
        except OSError:
            print(f'Error: failed to copy {most_recent} as new {file}.')


def _translate_timestamp(time_stamp):
    """Return a 'DD.MM.YY hh:mm:ss' timestamp from a YYMMDD-hhmmss one."""
    if len(time_stamp) != 13:
        raise ValueError('Error translating timestamp: Invalid length '
                         f'{len(time_stamp)}. Expected 13 characters.')
    year, month, day = time_stamp[:2], time_stamp[2:4], time_stamp[4:6]
    hour, minutes, secs = time_stamp[7:9], time_stamp[9:11], time_stamp[11:13]
    return f'{day}.{month}.{year} {hour}:{minutes}:{secs}'


def _workhistory_has_dirs_to_move(work_history_path):
    """Return whether work_history_path contains any directory worth moving."""
    work_history_dirs = (d for d in work_history_path.glob('r*')                # TODO: is this correct? In _move_workhistory_folders we also use HIST_FOLDER_RE, and folders contain 'r' but do not begin with it
                         if d.is_dir()
                         and PREVIOUS_LABEL not in d.name)
    return any(work_history_dirs)


def bookkeeper(mode,
               job_name=None,
               history_name=DEFAULT_HISTORY,
               work_history_name=DEFAULT_WORK_HISTORY):
    """Archive or discard the most recent run found in the current directory.

    Parameters
    ----------
    mode : str or BookkeeperMode
        Which bookkeeper mode to use. See help(BookkeeperMode).
    job_name : str or None, optional
        Custom name to append to the stored folder and to history.info.
        If not given or None, no extra name is added. Default is None.
    history_name : str, optional
        The name of the folder in the current directory where the
        most recent run should be archived. Default is 'history'.
    work_history_name : str, optional
        The name of the workhistory subfolder of ./work where
        intermediate runs may have been stored. Results are also
        copied from here to ./history_name. Default is workhistory.

    Returns
    -------
    exit_code : int
        0 -> success
        1 -> No files to be copied

    Raises
    ------
    OSError
        If creation of the history folder or any of the subfolders
        where results are to be stored fails.
    """

    # convert mode to enum if necessary
    _mode = BookkeeperMode(mode)

    # Get paths for history and workhistory
    cwd = Path.cwd()
    history_path = cwd / history_name
    work_history_path = cwd / work_history_name

    # Make list of stuff to move:
    # Log files (calc in root, others in SUPP), and SUPP/OUT folders
    files_to_move, logs_to_move = _collect_log_files(cwd)
    files_to_move.extend(_collect_supp_and_out(cwd))

    # Exit early if there's nothing to move:
    if (not files_to_move
            and not _workhistory_has_dirs_to_move(work_history_path)):
        print('Bookkeeper: Found nothing to do. Exiting...')
        return 1

    # Make history folder if not there yet
    try:
        history_path.mkdir(exist_ok=True)
    except OSError:
        print('Error creating history folder.')
        raise

    # Figure out the number of the tensor (it's the most recent one)
    # and the highest run number currently stored for each tensor in
    # history_path
    tensor_number = leedbase.getMaxTensorIndex(home=cwd, zip_only=True)
    max_nums = _find_max_run_per_tensor(history_path)

    if tensor_number not in max_nums and _mode.discard:
        # New Tensor to be discarded.
        _discard_tensors_and_deltas(cwd, tensor_number)

    # Infer timestamp from log file, if possible
    old_timestamp, last_log_lines = _read_most_recent_log(cwd)

    # Get new history subfolder tensor_dir
    tensor_dir = _create_new_history_directory(
        history_path,
        tensor_number,
        job_num=max_nums[tensor_number] + 1,
        suffix=old_timestamp + ('' if job_name is None else f'_{job_name}'),
        should_mkdir=not _mode.discard
        )
    if not mode.discard:
        store_input_files_to_history(cwd, tensor_dir)

    if _mode is BookkeeperMode.CONT:
        _replace_input_files_from_out(cwd)

    # Move (or discard) old stuff: files go to main history, logs go
    # to SUPP (except main viperleed-calc log); collect also folders
    # from workhistory
    _move_or_discard_files(files_to_move, tensor_dir, _mode.discard)
    _move_or_discard_files(logs_to_move, tensor_dir/'SUPP', _mode.discard)
    tensor_nums = _move_and_cleanup_workhistory(work_history_path,
                                                history_path,
                                                old_timestamp,
                                                max_nums,
                                                _mode.discard)
    tensor_nums.add(tensor_number)

    if not _mode.discard:  # write history.info, including notes
        _add_entry_to_history_info_file(
            cwd,
            tensor_nums,
            [max_nums[tensor] + 1 for tensor in tensor_nums],
            job_name,
            last_log_lines,
            old_timestamp,
            tensor_dir.name
            )
    return 0


class BookkeeperCLI(ViPErLEEDCLI, cli_name='bookkeeper'):
    """The main command-line interface for the bookkeeper utility."""

    def add_parser_arguments(self, parser):
        """Add bookkeeper arguments to parser."""
        super().add_parser_arguments(parser)
        what_next = parser.add_mutually_exclusive_group()
        what_next.add_argument(
            '-c', '--cont',
            help=('overwrite POSCAR and VIBROCC with POSCAR_OUT and '
                  'VIBROCC_OUT from the OUT folder, if present.'),
            action='store_true'
            )
        what_next.add_argument(
            '-d', '--discard',
            help=('discard all results from the last run, as if it had '
                  'not happened, and do not add anything to history or '
                  'history.info. Note that this will not necessarily '
                  'restore modified input files in the main folder.'),
            action='store_true',
            )
        parser.add_argument(
            '-j', '--job-name',
            help=('define a string to be appended to the name of the history '
                  'folder that is created, and is logged in history.info'),
            type=str
            )
        parser.add_argument(
            '--history-name',
            help=('define the name of the history folder that is '
                  f'created/used. Default is {DEFAULT_HISTORY!r}'),
            type=str,
            default=DEFAULT_HISTORY
            )
        parser.add_argument(
            '--work-history-name',
            help=('define the name of the workhistory folder that is '
                  f'created/used. Default is {DEFAULT_WORK_HISTORY!r}'),
            type=str,
            default=DEFAULT_WORK_HISTORY
            )

    def __call__(self, args=None):
        """Call the bookkeeper with command-line args."""
        parsed_args = self.parse_cli_args(args)

        # Select mode
        if parsed_args.cont:
            mode = BookkeeperMode.CONT
        elif parsed_args.discard:
            mode = BookkeeperMode.DISCARD
        else:   # default
            mode = BookkeeperMode.DEFAULT
        return bookkeeper(mode,
                          job_name=parsed_args.job_name,
                          history_name=parsed_args.history_name,
                          work_history_name=parsed_args.work_history_name,)


if __name__ == '__main__':
    BookkeeperCLI.run_as_script()
