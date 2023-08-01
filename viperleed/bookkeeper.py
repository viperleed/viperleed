# -*- coding: utf-8 -*-
"""
ViPErLEED bookkeeper module.

Created on Thu Jan 30 11:12:30 2020

@author: Florian Kraushofer
@author: Alexander M. Imre
"""

from enum import Enum, auto
from pathlib import Path
import argparse
import os
import re
import shutil
import time

from viperleed.tleedmlib.sections._sections import ALL_INPUT_FILES

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
    original_inputs_path = _root_path / "work" / "original_inputs"
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
               history_name="history",
               work_history_name="workhistory",):

    # convert mode to enum if necessary
    _mode = BookkeeperMode(mode)

    # get paths for history and workhistory
    history_path = Path(history_name).resolve()
    work_history_path = Path(work_history_name).resolve()
    tensors_path = Path("Tensors").resolve()
    deltas_path = Path("Deltas").resolve()

    # make list of stuff to move
    files_to_move = [d for d in os.listdir() if os.path.isdir(d)
              and (d == "OUT" or d == "SUPP")]
    # logs will be saved to history; tleedm in root, others in SUPP
    logs_to_move = []
    for file in os.listdir():
        if os.path.isfile(file) and file.endswith(".log"):
            if file.startswith("tleedm"):
                files_to_move.append(file)
            else:
                logs_to_move.append(file)
    # if there's nothing to move, return.
    if len(files_to_move) == 0:
        found = False
        # check workhist folder:
        if os.path.isdir(work_history_name):
            work_history_dirs = [d for d in os.listdir(work_history_name) if
                            os.path.isdir(os.path.join(work_history_name, d))
                            and d.startswith("r") and not ("previous" in d)]
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
        for f in [f for f in os.listdir("Tensors")
                  if ((tensors_path / f).is_file()
                      and rgx.match(f))]:
            m = rgx.match(f)
            if m.span()[1] == 15:  # exact match
                indlist.append(int(m.group(0)[-7:-4]))
        if indlist:
            tensor_number = max(indlist)
    # figure out the number of the run
    dl = [n for n in os.listdir(history_name)
          if os.path.isdir(os.path.join(history_name, n))]
    max_nums = {}  # max. job number per tensor number
    rgx = re.compile(r't[0-9]{3}.r[0-9]{3}_')
    for d in dl:
        m = rgx.match(d)
        if m:
            try:
                t = int(d[1:4])
                i = int(d[6:9])
                if t not in max_nums:
                    max_nums[t] = i
                else:
                    max_nums[t] = max(max_nums[t], i)
            except (ValueError, IndexError):
                pass
    if tensor_number not in max_nums:
        num = 1  # Tensor is new - if discard: delete
        if _mode is BookkeeperMode.DISCARD:
            tensor_file = tensors_path / "Tensors_{tensor_number:03d}.zip"
            delta_file = deltas_path / "Deltas_{tensor_number:03d}.zip"
            for f in (tensor_file, delta_file):
                if os.path.isfile(f):
                    try:
                        os.remove(f)
                    except Exception:
                        print(f"Failed to discard file {f}")
    else:
        num = max_nums[tensor_number] + 1
    # find old timestamp, if possible
    old_log_files = sorted([f for f in os.listdir() if os.path.isfile(f) and
                            f.endswith(".log") and f.startswith("tleedm-")])
    lastLogLines = []
    if len(old_log_files) > 0:
        old_timestamp = old_log_files[-1][7:20]
        try:
            with open(old_log_files[-1], "r") as rf:
                lastLogLines = rf.readlines()
        except Exception:
            pass
    else:
        timestamp = time.strftime("%y%m%d-%H%M%S", time.localtime())
        old_timestamp = "moved-" + timestamp
    if _mode is not BookkeeperMode.DISCARD:
        # get dirname
        dirname = f"t{tensor_number:03d}.r{num:03d}_{old_timestamp}"
        if job_name is not None:
            dirname += "_" + job_name
        tensor_dir = os.path.join(history_name, dirname)
        if os.path.isdir(tensor_dir):
            tensor_dir_2 = tensor_dir+"_moved-"+timestamp
            print(f"Error: Target directory {tensor_dir} already exists. Will "
                  "use {tensor_dir_2} instead.")
            tensor_dir = tensor_dir_2
            dirname = os.path.basename(tensor_dir)
        try:
            os.mkdir(tensor_dir)
        except Exception:
            print("Error: Could not create target directory " + tensor_dir
                  + "\n Stopping...")
            raise 1
        # copy (or discard) files
        cwd_path = Path(".")
        store_input_files_to_history(cwd_path, tensor_dir)

    # if CONT, check for POSCAR_OUT / VIBROCC_OUT
    if _mode is BookkeeperMode.CONT:
        if os.path.isdir("OUT"):
            file_out = sorted([f for f in os.listdir("OUT")
                           if os.path.isfile(os.path.join("OUT", f))
                           and f.startswith("POSCAR_OUT_")
                           and "parabola" not in f])
            if len(file_out) > 0:
                path = os.path.join("OUT", file_out[-1])
                try:
                    shutil.copy2(path, "POSCAR")
                except Exception:
                    print("Error: failed to copy "+path+" as new POSCAR.")
            else:
                print("Error: Flag --cont was set, but no POSCAR_OUT was "
                      "found in OUT.")
            file_out = sorted([f for f in os.listdir("OUT")
                           if os.path.isfile(os.path.join("OUT", f))
                           and f.startswith("VIBROCC_OUT_")
                           and "parabola" not in f])
            if len(file_out) > 0:
                path = os.path.join("OUT", file_out[-1])
                try:
                    shutil.copy2(path, "VIBROCC")
                except Exception:
                    print("Error: failed to copy "+path+" as new VIBROCC.")
            else:
                print("Error: Flag --cont was set, but no VIBROCC_OUT was "
                      "found in OUT.")
        else:
            print("Error: Flag --cont was set, but no OUT folder exists.")
    # move old stuff
    for f in files_to_move:
        if _mode is not BookkeeperMode.DISCARD:
            try:
                shutil.move(f, os.path.join(tensor_dir, f))
            except Exception:
                print("Error: Failed to move "+f)
        else:  # delete instead
            try:
                shutil.rmtree(f)
            except NotADirectoryError:
                try:
                    os.remove(f)
                except Exception:
                    print("Failed to discard file " + f)
            except Exception:
                print("Failed to discard directory " + f)
    # move log files to SUPP (except for general log tleedm....log)
    for log_file in logs_to_move:
        if _mode is not BookkeeperMode.DISCARD:
            try:
                supp_path = os.path.join(tensor_dir, 'SUPP')
                shutil.move(log_file, os.path.join(supp_path, log_file))
            except Exception:
                print("Error: Failed to move " + log_file)
        else:   # delete instead
            try:
                os.remove(log_file)
            except Exception:
                print("Failed to discard file " + log_file)

    # if there is a workhist folder, go through it and move contents as well
    tensor_nums = {tensor_number}
    if os.path.isdir(work_history_name) and _mode is not BookkeeperMode.DISCARD:
        work_hist_prev = [d for d in os.listdir(work_history_name) if
                        os.path.isdir(os.path.join(work_history_name, d))
                        and rgx.match(d) and ("previous" in d)]
        for d in work_hist_prev:
            try:
                shutil.rmtree(os.path.join(work_history_name, d))
            except Exception:
                print(f"Failed to delete {d} directory from "
                      f"{work_history_name}")
        work_history_dirs = [d for d in os.listdir(work_history_name) if
                        os.path.isdir(os.path.join(work_history_name, d))
                        and rgx.match(d) and not ("previous" in d)
                        and old_timestamp in d]
        for d in work_history_dirs:
            try:
                tensor_num_2 = int(d[1:4])
                search_num = int(d[6:9])
            except (ValueError, IndexError):
                pass
            else:
                if tensor_num_2 not in max_nums:
                    num = 1
                else:
                    num = max_nums[tensor_num_2] + 1
                newname = (f"t{tensor_num_2:03d}.r{num:03d}.{search_num:03d}"
                           + d[9:])
                try:
                    shutil.move(os.path.join(work_history_name, d),
                                os.path.join(history_name, newname))
                except Exception:
                    print("Error: Failed to move "
                          + os.path.join(work_history_name, d))
                tensor_nums.add(tensor_num_2)
    if os.path.isdir(work_history_name):
        if (len(os.listdir(work_history_name)) == 0
            or _mode is BookkeeperMode.DISCARD):
            try:
                shutil.rmtree(work_history_name)
            except Exception as e:
                if _mode is BookkeeperMode.DISCARD:
                    print(f"Failed to discard workhistory folder: {e}")
                else:
                    print(f"Failed to delete empty {work_history_name} "
                          f"directory: {str(e)}")
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
            with open(notes_name, 'r') as rf:
                notes = rf.read()
        except Exception:
            print("Error: Failed to read " + notes_name + " file")
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
    if len(lastLogLines) > 0:
        # try to read what was executed
        run_info = ""
        i = len(lastLogLines) - 1
        while i > 0:
            if lastLogLines[i].startswith("Executed segments: "):
                run_info = (lastLogLines[i].split("Executed segments: ")[1]
                           .strip())
                break
            i -= 1
        if run_info:
            hist += "# RUN ".ljust(spacing) + run_info + "\n"
        # now try to read final R-factors
        for j in range(i+1, len(lastLogLines)):
            line = lastLogLines[j]
            if line.startswith("Final R"):
                t = ""
                if "refcalc" in line:
                    t = "# R REF ".ljust(spacing)
                elif "superpos" in line:
                    t = "# R SUPER ".ljust(spacing)
                if t:
                    hist += t + line.split(":", maxsplit=1)[1].strip() + "\n"

    hist += "# TIME ".ljust(spacing) + f"{_translate_timestamp(old_timestamp)} \n"
    hist += "# FOLDER ".ljust(spacing) + f"{dirname} \n"
    hist += f"Notes: {notes}\n"
    hist += "\n###########\n"
    try:
        with open("history.info", "a") as wf:
            wf.write(hist)
    except Exception:
        print("Error: Failed to append to history.info file.")
    return 0

def _bookkeeper_cli_options():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--cont",
        help=("overwrite POSCAR and VIBROCC with POSCAR_OUT and VIBROCC_OUT "
              "from the OUT folder, if present."),
        action="store_true")
    parser.add_argument(
        "-d", "--discard",
        help=("discard all results from the last run, as if it had not "
              "happened, and do not add anything to history or history.info. "
              "Note that this will not necessarily restore modified input "
              "files in the main folder."),
        action="store_true")
    parser.add_argument(
        "-j", "--job_name",
        help=("defines a string to be appended to the name of the history "
              "folder that is created, and is logged in history.info"),
        type=str)
    parser.add_argument(
        "--history_name",
        help=("defines the name of the history folder that is created/used. "
              "Default is 'history'"),
        type=str,
        default="history")
    parser.add_argument(
        "--work_history_name",
        help=("defines the name of the workhistory folder that is created/used. "
              "Default is 'workhistory'"),
        type=str,
        default="workhistory")
    return parser.parse_args()


if __name__ == '__main__':
    # parse command line arguments

    args = _bookkeeper_cli_options()

    # select mode
    if args.cont and args.discard:
        raise RuntimeError("Flags --cont and --discard are incompatible.")
    elif args.cont:
        mode = BookkeeperMode.CONT
    elif args.discard:
        mode = BookkeeperMode.DISCARD
    else:   # default
        mode = BookkeeperMode.DEFAULT

    bookkeeper(mode,
               job_name=args.job_name,
               history_name=args.history_name,
               work_history_name=args.work_history_name,)
