# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 17:45 2019

@author: Florian Kraushofer

Reads a PHASESHIFTS file and allows the user to copy and rearrange the blocks.
"""

import time
import os
import sys

cd = os.path.realpath(os.path.dirname(__file__))
# NB: it's necessary to add vpr_path to sys.path so that viperleed
#     can be loaded correctly at the top-level package
vpr_path = os.path.realpath(os.path.join(cd, '..', '..'))
for import_path in (cd, vpr_path):
    if import_path not in sys.path:
        sys.path.append(import_path)

from viperleed.tleedmlib.files.phaseshifts import (readPHASESHIFTS,
                                                   writePHASESHIFTS)


###############################################
#                  MAIN                       #
###############################################

def main():
    # print some info
    print("This utility reads a phaseshifts file (eg PHASESHIFTS) and allows "
          "the user to copy and rearrange the blocks, such that they fit the "
          "order of elements in POSCAR and the order of sites as defined by "
          "ELSPLIT and SITEDEF in the PARAMETERS file. "
          "Check the documentation for further information on how to arrange "
          "blocks in the PHASESHIFTS file.\n")

    # read the phaseshifts file
    try:
        (firstline, ps, _, _) = readPHASESHIFTS(None, None, check=False)
    except FileNotFoundError:
        print("PHASESHIFTS file not found.")
        filename = ""
        while filename == "":
            filename = input("Enter phaseshift file name: ")
            if not filename:
                print("Input failed. Please try again.")
            else:
                try:
                    (firstline, ps, _, _) = readPHASESHIFTS(
                        None, None, readfile=filename, check=False)
                    print("Phaseshifts file read successfully.")
                except FileNotFoundError:
                    print(filename+" not found.")
                    filename = ""
                except Exception:
                    print("Exception while reading phaseshifts file: ",
                          exc_info=True)
                    return 1
    except Exception:
        print("Exception while reading phaseshifts file: ", exc_info=True)
        return 1

    # print some info
    print("Found "+str(len(ps[0][1]))+" blocks. Enter a space-separated list "
          "of how the blocks should be arranged in the new file.\n\n"
          "Examples:\n"
          "'1 1 2 2': Print the first block twice, then the second block "
          "twice.\n"
          "'2 1': Swap the first an the second block\n"
          "'2 3': Print the second, then the third block, delete the first "
          "block.\n")

    # get input for new order
    neworder = ""
    while neworder == "":
        neworder = input("Enter new order: ")
        if not neworder:
            print("Input failed. Please try again.")
        else:
            try:
                ol = [int(s) for s in neworder.split()]
            except ValueError:
                print("Could not parse input. Please try again.")
                neworder = ""
            for i in ol:
                if not 0 < i < len(ps[0][1])+1:
                    print("Input out of bounds. Please try again.")
                    neworder = ""
                    break

    # now rearrange blocks as given in ol
    oldps = ps[:]
    ps = []
    for (en, oldenps) in oldps:
        enps = []
        for i in ol:
            enps.append(oldenps[i-1])
        ps.append((en, enps))

    # in firstline, change the block-count to the new number of blocks
    firstline = str(len(ps[0][1])).rjust(3) + firstline[3:]

    # write new file
    try:
        timestamp = time.strftime("%y%m%d-%H%M%S", time.localtime())
        fn = "PHASESHIFTS_mod_"+timestamp
        writePHASESHIFTS(firstline, ps, filename=fn)
        print("Wrote new phaseshifts file as "+fn)
    except Exception:
        print("Error writing new phaseshifts file.")
        return 1
    return 0


if __name__ == "__main__":
    main()
