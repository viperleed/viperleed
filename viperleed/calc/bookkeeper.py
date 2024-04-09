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
from viperleed.calc.sections.calc_section import ALL_INPUT_FILES
from viperleed.calc.lib import leedbase
from viperleed.cli_base import ViPErLEEDCLI

_CALC_LOG_PREFIXES = (
    LOG_PREFIX,
    'tleedm',   # For backwards compatibility
    )
HIST_FOLDER_RE = re.compile(
    r't(?P<tensor_num>[0-9]{3}).r(?P<job_num>[0-9]{3})_'
    )
_PREVIOUS_LABEL = 'previous'


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
                         and _PREVIOUS_LABEL not in d.name)
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
    # to SUPP (except main viperleed-calc log)
    _move_or_discard_files(files_to_move, tensor_dir, _mode.discard)
    _move_or_discard_files(logs_to_move, tensor_dir/'SUPP', _mode.discard)

    # if there is a workhist folder, go through it and move contents as well
    tensor_nums = {tensor_number}

    if work_history_path.is_dir() and not _mode.discard:
        work_hist_prev = [d for d in os.listdir(work_history_name) if
                        os.path.isdir(os.path.join(work_history_name, d))
                        and HIST_FOLDER_RE.match(d) and ("previous" in d)]
        for dir in work_hist_prev:
            try:
                shutil.rmtree(work_history_path / dir)
            except Exception:
                print(f"Failed to delete {dir} directory from "
                      f"{work_history_path}")
        work_history_dirs = [dir for dir in work_history_path.iterdir() if
                        (work_history_path / dir).is_dir()
                        and HIST_FOLDER_RE.match(dir.name)
                        and not ("previous" in dir.name)
                        and old_timestamp in dir.name]
        for dir in work_history_dirs:
            try:
                tensor_num_2 = int(dir.name[1:4])
                search_num = int(dir.name[6:9])
            except (ValueError, IndexError):
                pass
            else:
                if tensor_num_2 not in max_nums:
                    num = 1
                else:
                    num = max_nums[tensor_num_2] + 1
                newname = (f"t{tensor_num_2:03d}.r{num:03d}.{search_num:03d}"
                           + dir.name[9:])
                try:
                    shutil.move(os.path.join(work_history_name, dir),
                                os.path.join(history_name, newname))
                except OSError:
                    print(f"Error: Failed to move {work_history_path / dir}.")
                tensor_nums.add(tensor_num_2)
    if work_history_path.is_dir():
        if (len(list(work_history_path.iterdir())) == 0
            or _mode.discard):
            try:
                shutil.rmtree(work_history_name)
            except OSError as error:
                if _mode.discard:
                    print(f"Failed to discard workhistory folder: {error}")
                else:
                    print(f"Failed to delete empty {work_history_name} "
                          f"directory: {str(error)}")
    if _mode.discard:  # all done
        return 0
    job_nums = []
    for tensor_number in tensor_nums:
        if tensor_number not in max_nums:
            job_nums.append(1)
        else:
            job_nums.append(max_nums[tensor_number] + 1)
    # look for notes file
    notes_name = ""
    notes = ""
    for fn in ("notes", "notes.txt"):
        if os.path.isfile(fn):
            notes_name = fn
            break
    if notes_name:
        try:
            with open(notes_name, 'r') as read_file:
                notes = read_file.read()
        except Exception:
            print(f"Error: Failed to read {notes_name} file.")
    if notes:
        try:
            with open(notes_name, 'w'):
                pass
        except Exception:
            print(f"Error: Failed to clear the {notes_name} file after "
                  "reading.")
    # write history.info
    spacing = 12
    hist = ""
    if os.path.isfile("history.info"):
        hist += "\n\n"
    if tensor_nums == {0}:
        hist += "# TENSORS ".ljust(spacing) + "None\n"
    else:
        hist += "# TENSORS ".ljust(spacing) + str(tensor_nums)[1:-1] + "\n"
    hist += "# JOB ID ".ljust(spacing) + str(job_nums)[1:-1] + "\n"
    if job_name is not None:
        hist += f"# JOB NAME {job_name} \n"
    if len(last_log_lines) > 0:
        # try to read what was executed
        run_info = ""
        i = len(last_log_lines) - 1
        while i > 0:
            if last_log_lines[i].startswith("Executed segments: "):
                run_info = (last_log_lines[i].split("Executed segments: ")[1]
                           .strip())
                break
            i -= 1
        if run_info:
            hist += "# RUN ".ljust(spacing) + run_info + "\n"
        # now try to read final R-factors
        for j in range(i+1, len(last_log_lines)):
            line = last_log_lines[j]
            if line.startswith("Final R"):
                r_fac_line = ""
                if "refcalc" in line:
                    r_fac_line = "# R REF ".ljust(spacing)
                elif "superpos" in line:
                    r_fac_line = "# R SUPER ".ljust(spacing)
                if r_fac_line:
                    hist += r_fac_line + line.split(":", maxsplit=1)[1].strip() + "\n"

    hist += "# TIME ".ljust(spacing) + f"{_translate_timestamp(old_timestamp)} \n"
    hist += "# FOLDER ".ljust(spacing) + f"{tensor_dir.name} \n"
    hist += f"Notes: {notes}\n"
    hist += "\n###########\n"
    try:
        with open("history.info", "a") as hist_info_file:
            hist_info_file.write(hist)
    except Exception:
        print("Error: Failed to append to history.info file.")
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
