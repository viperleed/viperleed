# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:12:30 2020

@author: Florian Kraushofer, Alexander Imre
"""

import os
import re
import time
import argparse
import shutil
from pathlib import Path

_INPUT_FILES = ["POSCAR", "PHASESHIFTS", "PARAMETERS", "IVBEAMS",
                 "DISPLACEMENTS", "VIBROCC", "EXPBEAMS.csv", "EXPBEAMS"]

def translate_timestamp(s):
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
                     if file.name in _INPUT_FILES)
    # copy files to history
    for file in files_to_copy:
        try:
            shutil.copy2(file, _history_path/file.name)
        except OSError as error_msg:
            print(f"Failed to copy file {file} to history: {error_msg}")

def bookkeeper():
    history_name = "history"  # name of history folder in home dir
    workhistname = "workhistory"  # name of history folder in work dir
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
        "-n", "--name",
        help=("defines a string to be appended to the name of the history "
              "folder that is created, and is logged in history.info"),
        type=str)
    args = parser.parse_args()
    if args.cont and args.discard:
        print("Bookkeeper ERROR: Flags --cont and --discard are incompatible."
              "Bookkeeper will stop.")
        return 1
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
        if os.path.isdir(workhistname):
            workhistdirs = [d for d in os.listdir(workhistname) if
                            os.path.isdir(os.path.join(workhistname, d))
                            and d.startswith("r") and not ("previous" in d)]
            if len(workhistdirs) > 0:
                found = True
        if not found:
            print("Bookkeeper: Found nothing to do. Exiting...")
            return 1
    # check whether history folder is there. If not, make one
    if not os.path.isdir(history_name):
        try:
            os.mkdir(history_name)
        except Exception:
            print("Error creating history folder.")
            raise
    # figure out the number of the tensor
    tnum = 0
    if os.path.isdir("Tensors"):
        indlist = []
        rgx = re.compile(r'Tensors_[0-9]{3}\.zip')
        for f in [f for f in os.listdir("Tensors")
                  if (os.path.isfile(os.path.join("Tensors", f))
                      and rgx.match(f))]:
            m = rgx.match(f)
            if m.span()[1] == 15:  # exact match
                indlist.append(int(m.group(0)[-7:-4]))
        if indlist:
            tnum = max(indlist)
    # figure out the number of the run
    dl = [n for n in os.listdir(history_name)
          if os.path.isdir(os.path.join(history_name, n))]
    maxnums = {}  # max. job number per tensor number
    rgx = re.compile(r't[0-9]{3}.r[0-9]{3}_')
    for d in dl:
        m = rgx.match(d)
        if m:
            try:
                t = int(d[1:4])
                i = int(d[6:9])
                if t not in maxnums:
                    maxnums[t] = i
                else:
                    maxnums[t] = max(maxnums[t], i)
            except (ValueError, IndexError):
                pass
    if tnum not in maxnums:
        num = 1  # Tensor is new - if discard: delete
        if args.discard:
            tensorfile = os.path.join("Tensors", "Tensors_{:03d}.zip"
                                      .format(tnum))
            deltafile = os.path.join("Deltas", "Deltas_{:03d}.zip"
                                     .format(tnum))
            for f in (tensorfile, deltafile):
                if os.path.isfile(f):
                    try:
                        os.remove(f)
                    except Exception:
                        print("Failed to discard file " + f)
    else:
        num = maxnums[tnum] + 1
    # find old timestamp, if possible
    oldlogfiles = sorted([f for f in os.listdir() if os.path.isfile(f) and
                          f.endswith(".log") and f.startswith("tleedm-")])
    lastLogLines = []
    if len(oldlogfiles) > 0:
        oldTimeStamp = oldlogfiles[-1][7:20]
        try:
            with open(oldlogfiles[-1], "r") as rf:
                lastLogLines = rf.readlines()
        except Exception:
            pass
    else:
        timestamp = time.strftime("%y%m%d-%H%M%S", time.localtime())
        oldTimeStamp = "moved-" + timestamp
    if not args.discard:
        # get dirname
        dirname = "t{:03d}.r{:03d}_".format(tnum, num) + oldTimeStamp
        if args.name:
            dirname += "_" + args.name
        tdir = os.path.join(history_name, dirname)
        if os.path.isdir(tdir):
            tdir2 = tdir+"_moved-"+timestamp
            print("Error: Target directory " + tdir + " already exists. Will "
                  "use " + tdir2 + " instead.")
            tdir = tdir2
            dirname = os.path.basename(tdir)
        try:
            os.mkdir(tdir)
        except Exception:
            print("Error: Could not create target directory " + tdir
                  + "\n Stopping...")
            raise 1
        # copy (or discard) files
        cwd_path = Path(".")
        store_input_files_to_history(cwd_path, tdir)

    # if CONT, check for POSCAR_OUT / VIBROCC_OUT
    if args.cont and not args.discard:
        if os.path.isdir("OUT"):
            fout = sorted([f for f in os.listdir("OUT")
                           if os.path.isfile(os.path.join("OUT", f))
                           and f.startswith("POSCAR_OUT_")
                           and "parabola" not in f])
            if len(fout) > 0:
                path = os.path.join("OUT", fout[-1])
                try:
                    shutil.copy2(path, "POSCAR")
                except Exception:
                    print("Error: failed to copy "+path+" as new POSCAR.")
            else:
                print("Error: Flag --cont was set, but no POSCAR_OUT was "
                      "found in OUT.")
            fout = sorted([f for f in os.listdir("OUT")
                           if os.path.isfile(os.path.join("OUT", f))
                           and f.startswith("VIBROCC_OUT_")
                           and "parabola" not in f])
            if len(fout) > 0:
                path = os.path.join("OUT", fout[-1])
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
        if not args.discard:
            try:
                shutil.move(f, os.path.join(tdir, f))
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
        if not args.discard:
            try:
                supp_path = os.path.join(tdir, 'SUPP')
                shutil.move(log_file, os.path.join(supp_path, log_file))
            except Exception:
                print("Error: Failed to move " + log_file)
        else:   # delete instead
            try:
                os.remove(log_file)
            except Exception:
                print("Failed to discard file " + log_file)

    # if there is a workhist folder, go through it and move contents as well
    tensornums = {tnum}
    if os.path.isdir(workhistname) and not args.discard:
        workhistprev = [d for d in os.listdir(workhistname) if
                        os.path.isdir(os.path.join(workhistname, d))
                        and rgx.match(d) and ("previous" in d)]
        for d in workhistprev:
            try:
                shutil.rmtree(os.path.join(workhistname, d))
            except Exception:
                print("Failed to delete "+d+" directory from "+workhistname)
        workhistdirs = [d for d in os.listdir(workhistname) if
                        os.path.isdir(os.path.join(workhistname, d))
                        and rgx.match(d) and not ("previous" in d)
                        and oldTimeStamp in d]
        for d in workhistdirs:
            try:
                tnum2 = int(d[1:4])
                snum = int(d[6:9])
            except (ValueError, IndexError):
                pass
            else:
                if tnum2 not in maxnums:
                    num = 1
                else:
                    num = maxnums[tnum2] + 1
                newname = ("t{:03d}.r{:03d}.{:03d}".format(tnum2, num, snum)
                           + d[9:])
                try:
                    shutil.move(os.path.join(workhistname, d),
                                os.path.join(history_name, newname))
                except Exception:
                    print("Error: Failed to move "
                          + os.path.join(workhistname, d))
                tensornums.add(tnum2)
    if os.path.isdir(workhistname):
        if len(os.listdir(workhistname)) == 0 or args.discard:
            try:
                shutil.rmtree(workhistname)
            except Exception as e:
                if args.discard:
                    print(f"Failed to discard workhistory folder: {e}")
                else:
                    print(f"Failed to delete empty {workhistname} directory: "
                          f"{str(e)}")
    if args.discard:  # all done
        return 0
    jobnums = []
    for tnum in tensornums:
        if tnum not in maxnums:
            jobnums.append(1)
        else:
            jobnums.append(maxnums[tnum] + 1)
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
            print("Error: Failed to clear the " + notes_name + " file after "
                  "reading.")
    # write history.info
    spacing = 12
    hi = ""
    if os.path.isfile("history.info"):
        hi += "\n\n"
    if tensornums == {0}:
        hi += "# TENSORS ".ljust(spacing) + "None\n"
    else:
        hi += "# TENSORS ".ljust(spacing) + str(tensornums)[1:-1] + "\n"
    hi += "# JOB ID ".ljust(spacing) + str(jobnums)[1:-1] + "\n"
    if args.name:
        hi += "# JOB NAME " + args.name + "\n"
    if len(lastLogLines) > 0:
        # try to read what was executed
        runinfo = ""
        i = len(lastLogLines) - 1
        while i > 0:
            if lastLogLines[i].startswith("Executed segments: "):
                runinfo = (lastLogLines[i].split("Executed segments: ")[1]
                           .strip())
                break
            i -= 1
        if runinfo:
            hi += "# RUN ".ljust(spacing) + runinfo + "\n"
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
                    hi += t + line.split(":", maxsplit=1)[1].strip() + "\n"

    hi += "# TIME ".ljust(spacing) + translate_timestamp(oldTimeStamp) + "\n"
    hi += "# FOLDER ".ljust(spacing) + dirname + "\n"
    hi += "Notes: " + notes + "\n"
    hi += "\n###########\n"
    try:
        with open("history.info", "a") as wf:
            wf.write(hi)
    except Exception:
        print("Error: Failed to append to history.info file.")
    return 0


if __name__ == '__main__':
    bookkeeper()
