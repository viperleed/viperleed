# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 17:29:24 2021

@author: Florian Kraushofer

Quick-and-dirty job script to run viperleed in a given work directory.

Usage:
  - place this script in a new folder with input files
  - define the globals vpr_path and work_path manually in this script
  - place bookkeeper.py in the same folder if you want it to run automatically
  - run (precise behaviour controlled by PARAMETERS)
vpr_path is where to find the viperleed source code
work_path is where to run tleedm (note: creates lots of files!)

"""

import os
import sys
import shutil

try:
    from bookkeeper import bookkeeper
    bookie_exists = True
except ModuleNotFoundError:
    bookie_exists = False

# path to directory containing viperleed source - define explicitly
vpr_path = "/home/path/to/source/"       # without final /viperleed
work_path = os.path.join(".", "work")    # where to run

delete_workdir = False   # delete the work_path after copying back?
all_tensors = False      # copy all tensor files or just highest number?
# !!! TODO: it would be nice if all_tensors automatically checked PARAMETERS

if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

import viperleed
import viperleed.tleedmlib
from viperleed.tleedm import run_tleedm


def main():
    global vpr_path
    global work_path
    global all_tensors
    global delete_workdir

    # first, run bookkeeper (if found)
    if bookie_exists:
        print("Running bookkeeper...")
        bookkeeper()

    # create work directory if necessary
    os.makedirs(work_path, exist_ok=True)

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
        tensor_num = viperleed.tleedmlib.leedbase.getMaxTensorIndex(
            zip_only=True)
        if tensor_num > 0:
            os.makedirs(os.path.join(work_path, "Tensors"), exist_ok=True)
            tensorfile = os.path.join("Tensors", "Tensors_{:03d}.zip"
                                      .format(tensor_num))
            shutil.copy2(tensorfile, os.path.join(work_path, tensorfile))
            deltafile = os.path.join("Deltas", "Deltas_{:03d}.zip"
                                     .format(tensor_num))
            if os.path.isfile(deltafile):
                os.makedirs(os.path.join(work_path, "Deltas"), exist_ok=True)
                shutil.copy2(deltafile, os.path.join(work_path, deltafile))

    # copy input files to work directory
    for file in ["PARAMETERS", "VIBROCC", "IVBEAMS", "DISPLACEMENTS", "POSCAR",
                 "PHASESHIFTS", "EXPBEAMS.csv"]:
        try:
            shutil.copy2(file, os.path.join(work_path, file))
        except FileNotFoundError:
            pass

    # go to work directory, execute there
    home = os.path.abspath(".")
    os.chdir(work_path)
    run_tleedm(source=os.path.join(vpr_path, "viperleed"))

    # copy back everything listed in manifest
    manifest = []
    if os.path.isfile("manifest"):
        with open("manifest", "r") as rf:
            manifest = [s.strip() for s in rf.readlines()]
    for p in manifest:
        try:
            if os.path.isfile(p):
                shutil.copy2(p, os.path.join(home, p))
            elif os.path.isdir(p):
                shutil.copytree(p, os.path.join(home, p), dirs_exist_ok=True)
        except Exception as e:
            print("Error copying " + p + " to home directory: " + str(e))

    # go back, clean up if requested
    os.chdir(home)
    if delete_workdir:
        try:
            shutil.rmtree(work_path)
        except Exception as e:
            print("Error deleting work directory: " + str(e))
    return


if __name__ == "__main__":
    main()
