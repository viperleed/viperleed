# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 17:45 2019

@author: Florian Kraushofer

Reads a POSCAR file, asks at what c value the bulk starts, then automatically
reduces the size of the POSCAR to non-redundant bulk layers only, and outputs
the appropriate N_BULK_LAYERS and BULK_REPEAT values.
"""

import copy

import numpy as np

from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.files import poscar
from viperleed.tleedmlib.files.woods_notation import writeWoodsNotation


def main():
    # print some info
    print("This utility reads a POSCAR file with arbitrary thickness and "
          "orientation of the bulk, requiring only an input where the bulk "
          "starts. The bulk repeat unit is then automatically determined. "
          "Output is a new POSCAR with a minimum amount of bulk, and values "
          "for the BULK_REPEAT and N_BULK_LAYERS parameters.\n")

    # read the POSCAR file
    filename = ""
    while not filename:
        filename = input("Enter POSCAR file name (Default: [POSCAR]): ")
        if not filename:
            filename = "POSCAR"
        try:
            sl = poscar.read(filename=filename)
        except FileNotFoundError:
            print(f"File {filename} not found.")
            filename = ""
        except Exception:
            print("Exception while reading POSCAR file")
            return 1

    # print some info
    print("Slab POSCAR was read successfully.")

    cut = -1
    while cut < 0:
        cutstr = input("Enter c value at which bulk starts: ")
        if cutstr != "":
            try:
                cut = float(cutstr)
            except ValueError:
                print("Failed to convert input to float.")
                cut = -1
            else:
                if not (0 <= cut < 1):
                    print("Error: c value has to be between 0 and 1.")
                    cut = -1

    eps = -1
    while eps < 0:
        s = input("Enter tolerances for symmetry search (Default: [0.1 A]): ")
        if s == "":
            eps = 0.1
        else:
            try:
                eps = float(s)
            except ValueError:
                print("Could not convert value to float.")
            else:
                if eps < 0:
                    print("Value has to be greater than zero.")

    rp = Rparams()
    rp.LAYER_CUTS = [cut]
    rp.N_BULK_LAYERS = 1
    rp.SYMMETRY_EPS = eps
    rp.SYMMETRY_EPS_Z = eps
    sl.fullUpdate(rp)

    sl.bulkslab = sl.makeBulkSlab(rp)
    bsl = sl.bulkslab
    bsl.createSublayers(eps)

    print("Checking bulk unit cell...")
    changecell, mincell = bsl.getMinUnitCell(rp)
    if changecell:
        sl.changeBulkCell(rp, mincell)
        bsl = sl.bulkslab
    if not rp.superlattice_defined:
        ws = writeWoodsNotation(rp.SUPERLATTICE)
        # !!! replace the writeWoodsNotation from baselib with
        #   the one from guilib
        si = rp.SUPERLATTICE.astype(int)
        if ws:
            print("Found SUPERLATTICE = "+ws)
        else:
            print("Found SUPERLATTICE M = {} {}, {} {}".format(
                si[0, 0], si[0, 1], si[1, 0], si[1, 1]))

    bsl.getCartesianCoordinates()
    newC = bsl.getMinC(rp, z_periodic=False)

    if newC is None:
        print("No repeat vector was found inside the bulk. The bulk may "
              "already be minimal.")
        return 0

    rp.BULK_REPEAT = -newC

    bsl.atlist = [at for at in bsl
                  if at.cartpos[2] > bsl.topat_ori_z - abs(newC[2])]
    bsl.layers[0].atlist = bsl.atlist

    rp.SUPERLATTICE = np.eye(2)
    newbsl = bsl.makeBulkSlab(rp)

    # find largest layer spacing in reduced POSCAR
    newbsl.createSublayers(eps)
    bulkcut = -1
    if newbsl.n_sublayers > 1:
        maxdist = abs(newbsl.sublayers[1].cartbotz
                      - newbsl.sublayers[0].cartbotz)
        cutlayer = 0
        for i in range(1, newbsl.n_sublayers - 1):
            d = abs(newbsl.sublayers[i+1].cartbotz
                    - newbsl.sublayers[i].cartbotz)
            if d > maxdist:
                maxdist = d
                cutlayer = i
        if maxdist < 1.2 and maxdist > 0.7:
            dec = None
            while dec is None:
                r = input("Cutting the bulk into two layers is possible with "
                          "a spacing of {} A. Proceed? (y/[n]): "
                          .format(round(maxdist, 2)))
                if r == "" or r.lower() in ["n", "no"]:
                    dec = False
                elif r.lower() in ["y", "yes"]:
                    dec = True
        elif maxdist > 1.2:
            dec = True
        else:
            dec = False
        if dec:
            bulkcut = newbsl.sublayers[cutlayer].cartbotz + maxdist/2

    # write POSCAR_bulk
    newbsl.sort_original()
    try:
        writePOSCAR(newbsl, filename='POSCAR_bulk', comments='none')
        print("Wrote POSCAR_bulk. Check file to see if periodicity is "
              "correct.")
    except Exception:
        print("Exception occurred while writing POSCAR_bulk")

    # create POSCAR with reduced size
    newsl = copy.deepcopy(sl)
    newsl.sort_by_z()
    topBulkAt = [at for at in newsl if at.pos[2] <= cut][-1]
    botSlabAt = [at for at in newsl if at.pos[2] > cut][0]
    fracRepeat = np.dot(np.linalg.inv(newsl.ucell), newC)
    newZero = (topBulkAt.pos[2]
               + (abs(botSlabAt.pos[2] - topBulkAt.pos[2]) / 2)
               - abs(fracRepeat[2]))
    newsl.atlist = [at for at in newsl if at.pos[2] > newZero]
    newsl.updateElementCount()   # update the number of atoms per element
    for at in newsl:
        at.pos[2] -= newZero
    newsl.getCartesianCoordinates()
    newsl.ucell[:, 2] = newsl.ucell[:, 2] * (1 - newZero)
    newsl.getFractionalCoordinates()
    newcut = round((topBulkAt.pos[2]
                   + (abs(botSlabAt.pos[2] - topBulkAt.pos[2]) / 2)), 3)
    if bulkcut > 0:
        newbulkcut = round((topBulkAt.pos[2] - (bulkcut / newsl.ucell[2, 2])),
                           3)
    # write POSCAR_min
    newsl.sort_original()
    try:
        writePOSCAR(newsl, filename='POSCAR_min', comments='none')
        print("Wrote POSCAR_min, to be used with parameters below.")
    except Exception:
        print("Exception occurred while writing POSCAR_min")

    # print info
    print("\nParameters found:")
    print("BULK_REPEAT = xyz"+str(rp.BULK_REPEAT))
    if bulkcut > 0:
        print("N_BULK_LAYERS = 2")
        print("LAYER_CUTS = {} {}   ! etc.".format(newbulkcut, newcut))
    else:
        print("N_BULK_LAYERS = 1")
        print("LAYER_CUTS = {}      ! etc.".format(newcut))

    input("Program finished, exit with return")

    return 0


if __name__ == "__main__":
    r = main()
    if r != 0:
        print("Exit with error")
