"""ViPErLEED utility: rearrange phaseshifts.

Reads a PHASESHIFTS file and allows the user to copy and rearrange the blocks.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-12-16'
__license__ = 'GPLv3+'

import time

from viperleed.calc.files.phaseshifts import (readPHASESHIFTS,
                                              writePHASESHIFTS)

###############################################
#                  MAIN                       #
###############################################

def main(args=None):
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
                    print(f"{filename} not found.")
                    filename = ""
                except Exception:
                    print("Exception while reading phaseshifts file: ")
                    return 1
    except Exception:
        print("Exception while reading phaseshifts file: ")
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
        fn = f"PHASESHIFTS_mod_{timestamp}"
        writePHASESHIFTS(firstline, ps, file_path=fn)
        print(f"Wrote new phaseshifts file as {fn}")
    except Exception:
        print("Error writing new phaseshifts file.")
        raise
    return 0


if __name__ == "__main__":
    main()
