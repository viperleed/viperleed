"""ViPErLEED bookkeeper module of package calc."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from enum import Enum
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
from viperleed.cli_base import ViPErLEEDCLI

_CALC_LOG_PREFIXES = (
    LOG_PREFIX,
    'tleedm',   # For backwards compatibility
    )


class BookkeeperMode(Enum):
    DEFAULT = 'default'  # store last run, but do not overwrite POSCAR, VIBROCC
    CONT = 'cont'        # store last run and overwrite POSCAR, VIBROCC from OUT
    DISCARD = 'discard'  # discard previous run as if it never happened


def _translate_timestamp(s):
    """Takes a timestamp YYMMDD-hhmmss and translates it to format DD.MM.YY
    hh:mm:ss"""
    if len(s) != 13:
        print("Error translating timestamp: Invalid length")
        return s
    return "{}.{}.{} {}:{}:{}".format(s[4:6], s[2:4], s[0:2],
                                      s[7:9], s[9:11], s[11:13])


def store_input_files_to_history(root_path, history_path):
    """Finds and copies input files to history.

    Parameters
    ----------
    root_path : pathlike
        Root directory from which to take files. Should be cwd, not ./work.
    history_path : pathlike
        Path to the history directory in which the files should be stored.
    """
    _root_path, _history_path = Path(root_path), Path(history_path)
    original_inputs_path = _root_path / DEFAULT_WORK / ORIGINAL_INPUTS_DIR_NAME
    if original_inputs_path.is_dir():
        input_origin_path = original_inputs_path
    else:
        input_origin_path = _root_path
        print("Could not find directory 'original_inputs' with unaltered "
              "input files. Files will instead be copied from the root "
              "directory.")

    # only files, no dirs
    files_in_dir =  (f for f in Path(input_origin_path).iterdir()
                     if f.is_file())
    files_to_copy = (file for file in files_in_dir
                     if file.name in ALL_INPUT_FILES)
    # copy files to history
    for file in files_to_copy:
        try:
            shutil.copy2(file, _history_path/file.name)
        except OSError as error_msg:
            print(f"Failed to copy file {file} to history: {error_msg}")


def bookkeeper(mode,
               job_name=None,
               history_name=DEFAULT_HISTORY,
               work_history_name=DEFAULT_WORK_HISTORY,):

    # convert mode to enum if necessary
    _mode = BookkeeperMode(mode)

    # get paths for history and workhistory
    history_path = Path(history_name).resolve()
    work_history_path = Path(work_history_name).resolve()
    tensors_path = Path("Tensors").resolve()
    deltas_path = Path("Deltas").resolve()
    out_path = Path("OUT").resolve()

    # make list of stuff to move
    files_to_move = [d for d in os.listdir() if os.path.isdir(d)
              and (d == "OUT" or d == "SUPP")]
    # logs will be saved to history; calc in root, others in SUPP
    logs_to_move = []
    for file in os.listdir():
        if os.path.isfile(file) and file.endswith(".log"):
            if file.startswith(_CALC_LOG_PREFIXES):
                files_to_move.append(file)
            else:
                logs_to_move.append(file)
    # if there's nothing to move, return.
    if len(files_to_move) == 0:
        found = False
        # check workhist folder:
        if work_history_path.is_dir():
            work_history_dirs = [d for d in work_history_path.iterdir() if
                                (work_history_path / d).is_dir()
                                and d.name.startswith("r")
                                and not ("previous" in d.name)]
            if len(work_history_dirs) > 0:
                found = True
        if not found:
            print("Bookkeeper: Found nothing to do. Exiting...")
            return 1
    # check whether history folder is there. If not, make one
    if not history_path.is_dir():
        try:
            os.mkdir(history_name)
        except Exception:
            print("Error creating history folder.")
            raise
    # figure out the number of the tensor
    tensor_number = 0
    if tensors_path.is_dir():
        indlist = []
        rgx = re.compile(r'Tensors_[0-9]{3}\.zip')
        for file in [f for f in os.listdir("Tensors")
                  if ((tensors_path / f).is_file()
                      and rgx.match(f))]:
            match = rgx.match(file)
            if match.span()[1] == 15:  # exact match
                indlist.append(int(match.group(0)[-7:-4]))
        if indlist:
            tensor_number = max(indlist)
    # figure out the number of the run
    dir_list = [d for d in history_path.iterdir()
                if (history_path / d).is_dir()]
    max_nums = {}  # max. job number per tensor number
    rgx = re.compile(r't[0-9]{3}.r[0-9]{3}_')
    for dir in dir_list:
        match = rgx.match(dir.name)
        if match:
            try:
                r_fac_line = int(dir.name[1:4])
                i = int(dir.name[6:9])
                if r_fac_line not in max_nums:
                    max_nums[r_fac_line] = i
                else:
                    max_nums[r_fac_line] = max(max_nums[r_fac_line], i)
            except (ValueError, IndexError):
                pass
    if tensor_number not in max_nums:
        num = 1  # Tensor is new - if discard: delete
        if _mode is BookkeeperMode.DISCARD:
            tensor_file = tensors_path / "Tensors_{tensor_number:03d}.zip"
            delta_file = deltas_path / "Deltas_{tensor_number:03d}.zip"
            for file in (tensor_file, delta_file):
                if os.path.isfile(file):
                    try:
                        os.remove(file)
                    except Exception:
                        print(f"Failed to discard file {file}")
    else:
        num = max_nums[tensor_number] + 1
    # find old timestamp, if possible
    old_log_files = sorted(
        f for f in os.listdir() if os.path.isfile(f) and
        f.endswith(".log") and f.startswith(_CALC_LOG_PREFIXES)
        )
    last_log_lines = []
    if len(old_log_files) > 0:
        old_timestamp = old_log_files[-1][-17:-4]
        try:
            with open(old_log_files[-1], "r") as read_file:
                last_log_lines = read_file.readlines()
        except Exception:
            pass
    else:
        timestamp = time.strftime("%y%m%d-%H%M%S", time.localtime())
        old_timestamp = f"moved-{timestamp}"
    if _mode is not BookkeeperMode.DISCARD:
        # get dirname
        dirname = f"t{tensor_number:03d}.r{num:03d}_{old_timestamp}"
        if job_name is not None:
            dirname += "_" + job_name
        tensor_dir = Path(history_name).resolve() / dirname
        if os.path.isdir(tensor_dir):
            tensor_dir_2 = f"{tensor_dir}_moved-{timestamp}"
            print(f"Error: Target directory {tensor_dir} already exists. Will "
                  "use {tensor_dir_2} instead.")
            tensor_dir = Path(tensor_dir_2).resolve()
            dirname = os.path.basename(tensor_dir)
        try:
            os.mkdir(tensor_dir)
        except Exception:
            print(f"Error: Could not create target directory {tensor_dir}"
                  "\n Stopping...")
            raise
        # copy (or discard) files
        cwd_path = Path(".")
        store_input_files_to_history(cwd_path, tensor_dir)

    # if CONT, check for POSCAR_OUT / VIBROCC_OUT
    # do not complain if not found, since we move previous runs to the history
    # by default
    if _mode is BookkeeperMode.CONT:
        if os.path.isdir("OUT"):
            for file in ("POSCAR", "VIBROCC"):

                files_in_out = sorted(
                    [f for f in out_path.iterdir()
                    if (out_path / f).is_file()
                    and f.name.startswith(f"{file}_OUT_")
                    and "parabola" not in f.name]
                )
                if len(files_in_out) > 0:
                    path = out_path / files_in_out[-1]
                    try:
                        shutil.copy2(path, file)
                    except Exception:
                        print(f"Error: failed to copy {path} as new {file}.")

    # move old stuff
    for file in files_to_move:
        if _mode is not BookkeeperMode.DISCARD:
            try:
                shutil.move(file, tensor_dir / file)
            except Exception:
                print(f"Error: Failed to move {file}.")
        else:  # delete instead
            try:
                shutil.rmtree(file)
            except NotADirectoryError:
                try:
                    os.remove(file)
                except Exception:
                    print(f"Failed to discard file {file}.")
            except Exception:
                print(f"Failed to discard directory {file}.")
    # move log files to SUPP (except for general viperleed-calc log)
    for log_file in logs_to_move:
        if _mode is not BookkeeperMode.DISCARD:
            try:
                hist_supp_path = tensor_dir / 'SUPP'
                shutil.move(log_file, hist_supp_path / log_file)
            except Exception:
                print(f"Error: Failed to move {log_file}.")
        else:   # delete instead
            try:
                os.remove(log_file)
            except Exception:
                print("Failed to discard file {log_file}.")

    # if there is a workhist folder, go through it and move contents as well
    tensor_nums = {tensor_number}

    if work_history_path.is_dir() and _mode is not BookkeeperMode.DISCARD:
        work_hist_prev = [d for d in os.listdir(work_history_name) if
                        os.path.isdir(os.path.join(work_history_name, d))
                        and rgx.match(d) and ("previous" in d)]
        for dir in work_hist_prev:
            try:
                shutil.rmtree(work_history_path / dir)
            except Exception:
                print(f"Failed to delete {dir} directory from "
                      f"{work_history_path}")
        work_history_dirs = [dir for dir in work_history_path.iterdir() if
                        (work_history_path / dir).is_dir()
                        and rgx.match(dir.name)
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
            or _mode is BookkeeperMode.DISCARD):
            try:
                shutil.rmtree(work_history_name)
            except OSError as error:
                if _mode is BookkeeperMode.DISCARD:
                    print(f"Failed to discard workhistory folder: {error}")
                else:
                    print(f"Failed to delete empty {work_history_name} "
                          f"directory: {str(error)}")
    if _mode is BookkeeperMode.DISCARD:  # all done
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
    hist += "# FOLDER ".ljust(spacing) + f"{dirname} \n"
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

