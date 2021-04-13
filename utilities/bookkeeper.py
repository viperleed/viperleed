# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:12:30 2020

@author: Florian Kraushofer
"""

import os
import re
import time
import argparse
import shutil


def translateTimestamp(s):
    """Takes a timestamp YYMMDD-hhmmss and translates it to format DD.MM.YY
    hh:mm:ss"""
    if len(s) != 13:
        print("Error translating timestamp: Invalid length")
        return s
    return "{}.{}.{} {}:{}:{}".format(s[4:6], s[2:4], s[0:2],
                                      s[7:9], s[9:11], s[11:13])


def main():
    histname = "history"  # name of history folder in home dir
    workhistname = "workhistory"  # name of history folder in work dir
    copyfiles = ["POSCAR", "PHASESHIFTS", "PARAMETERS", "IVBEAMS",
                 "DISPLACEMENTS", "VIBROCC", "EXPBEAMS.csv"]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--cont",
        help=("overwrite POSCAR and VIBROCC with POSCAR_OUT and VIBROCC_OUT "
              "from the OUT folder, if present"),
        action="store_true")
    parser.add_argument(
        "-n", "--name",
        help=("defines a string to be appended to the name of the history "
              "folder that is created, and is logged in history.info"),
        type=str)
    args = parser.parse_args()
    # make list of stuff to move
    tomove = [d for d in os.listdir() if os.path.isdir(d)
              and (d == "OUT" or d == "AUX")]
    tomove.extend([f for f in os.listdir() if os.path.isfile(f)
                   and f.endswith(".log")])
    # if there's nothing to move, return.
    if len(tomove) == 0:
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
    if not os.path.isdir(histname):
        try:
            os.mkdir(histname)
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
    dl = [n for n in os.listdir(histname)
          if os.path.isdir(os.path.join(histname, n))]
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
        num = 1
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
    # get dirname
    dirname = "t{:03d}.r{:03d}_".format(tnum, num) + oldTimeStamp
    if args.name:
        dirname += "_" + args.name
    tdir = os.path.join(histname, dirname)
    if os.path.isdir(tdir):
        tdir2 = tdir+"_moved-"+timestamp
        print("Error: Target directory " + tdir + " already exists. Will use "
              + tdir2 + " instead.")
        tdir = tdir2
        dirname = os.path.basename(tdir)
    try:
        os.mkdir(tdir)
    except Exception:
        print("Error: Could not create target directory " + tdir
              + "\n Stopping...")
        raise 1
    # copy files
    for f in [f for f in os.listdir() if f in copyfiles]:
        try:
            shutil.copy2(f, os.path.join(tdir, f))
        except Exception:
            print("Error: Failed to copy "+f)
    # if CONT, check for POSCAR_OUT / VIBROCC_OUT
    if args.cont:
        if os.path.isdir("OUT"):
            fout = sorted([f for f in os.listdir("OUT")
                           if os.path.isfile(os.path.join("OUT", f))
                           and f.startswith("POSCAR_OUT_")])
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
                           and f.startswith("VIBROCC_OUT_")])
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
    for f in tomove:
        try:
            shutil.move(f, os.path.join(tdir, f))
        except Exception:
            print("Error: Failed to move "+f)
    # if there is a workhist folder, go through it and move contents as well
    tensornums = {tnum}
    if os.path.isdir(workhistname):
        workhistprev = [d for d in os.listdir(workhistname) if
                        os.path.isdir(os.path.join(workhistname, d))
                        and rgx.match(d) and ("previous" in d)]
        for d in workhistprev:
            try:
                shutil.rmtree(d)
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
                                os.path.join(histname, newname))
                except Exception:
                    print("Error: Failed to move "
                          + os.path.join(workhistname, d))
                tensornums.add(tnum2)
    if os.path.isdir(workhistname):
        if len(os.listdir(workhistname)) == 0:
            try:
                shutil.rmtree(workhistname)
            except Exception:
                print("Failed to delete empty "+workhistname+" directory.")
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

    hi += "# TIME ".ljust(spacing) + translateTimestamp(oldTimeStamp) + "\n"
    hi += "# FOLDER ".ljust(spacing) + dirname + "\n"
    hi += "Notes: " + notes + "\n"
    hi += "\n###########\n"
    try:
        with open("history.info", "a") as wf:
            wf.write(hi)
    except Exception:
        print("Error: Failed to append to history.info file.")


if __name__ == '__main__':
    main()
