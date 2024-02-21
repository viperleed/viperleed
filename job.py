# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 17:29:24 2021

@author: Florian Kraushofer, Alexander Imre

Quick-and-dirty job script to run viperleed in a given work directory.

Usage:
  - place this script in a new folder with input files
  - define the globals vpr_path and work_path manually in this script
  - place bookkeeper.py in the same folder if you want it to run automatically
  - run (precise behaviour controlled by PARAMETERS)
vpr_path is where to find the viperleed source code - can be provided as command line argument
work_path is where to run tleedm (note: creates lots of files!) - can be provided as command line argument

"""

import os
import sys
import shutil
import argparse

# path to directory containing viperleed source - define explicitly here or pass as command line argument
vpr_path = None   # without final /viperleed - i.e. "/home/path/to/source/"
work_path = None    # where to run, without final /work - i.e. "."

parser = argparse.ArgumentParser()
parser.add_argument(
    "-s", "--source",
    help=("specify ViPErLEED source directory (without final /viperleed)"),
    type=str)
parser.add_argument(
    "-w", "--work",
    help=("specify work directory containing input files"),
    type=str)
args, bookie_args = parser.parse_known_args()
sys.argv = sys.argv[:1] + bookie_args

if args.source:
    vpr_path = args.source
if args.work:
    work_path = args.work

# If paths not supplied as command line argument - use explicit form if given above, otherwise raise Error
if not (vpr_path and work_path):
    raise ValueError("ViPErLEED source and/or work directory not defined!")

if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))

try:
    from viperleed.utilities.bookkeeper import bookkeeper
    bookie_exists = True
except ModuleNotFoundError:
    bookie_exists = False

delete_workdir = False   # delete the work_path after copying back?
all_tensors = False      # copy all tensor files or just highest number?
# !!! TODO: it would be nice if all_tensors automatically checked PARAMETERS


import viperleed
import viperleed.tleedmlib
from viperleed.tleedm import run_tleedm
from viperleed.tleedmlib.base import copytree_exists_ok


def main():
    global vpr_path
    global work_path
    global all_tensors
    global delete_workdir

    # first, run bookkeeper (if found)
    if bookie_exists:
        print("Running bookkeeper...")
        bookkeeper()

    work_path = os.path.join(work_path, "work")  # make /work subdirectory

    # create work directory if necessary
    os.makedirs(work_path, exist_ok=True)

    # copy Tensors and Deltas to work directory
    if all_tensors:
        try:
            copytree_exists_ok("Tensors", os.path.join(work_path, "Tensors"))
        except FileNotFoundError:
            pass
        try:
            copytree_exists_ok("Deltas", os.path.join(work_path, "Deltas"))
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
                 "PHASESHIFTS", "EXPBEAMS.csv", "EXPBEAMS"]:
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
                copytree_exists_ok(p, os.path.join(home, p))
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
