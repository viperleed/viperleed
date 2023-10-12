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

from viperleed.tleedmlib.classes.atom_containers import AtomList
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.classes.slab import NoBulkRepeatError
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
    rp.BULK_LIKE_BELOW = cut
    rp.SYMMETRY_EPS = rp.SYMMETRY_EPS_Z = eps

    print('Checking bulk unit cell...')
    try:
        new_bulk_c, bulk_cuts, bulk_interlayer = sl.detect_bulk(rp, 0.7)
    except NoBulkRepeatError:
        _no_repeat = True
    else:
        _no_repeat = False

    ws = writeWoodsNotation(rp.SUPERLATTICE)                                    # TODO: replace writeWoodsNotation with guilib functions
    if ws:
        info = f'= {ws}'
    else:
        info = 'M = {} {}, {} {}'.format(*rp.SUPERLATTICE.astype(int).ravel())
    print(f'Found SUPERLATTICE {info}')

    if _no_repeat:
        print('No repeat vector was found inside the bulk. The bulk may '
              'already be minimal.')
        return 0

    if 0.7 < bulk_interlayer < 1.2:
        dec = None
        while dec is None:
            reply = input('Cutting the bulk into two layers is possible '
                          f'with a spacing of {round(bulk_interlayer, 2)} A. '
                          'Proceed? (y/[n]): ')
            reply = reply.lower()
            if not reply or reply.startswith('n'):
                dec = False
            elif reply.startswith('y'):
                dec = True
        if not dec:
            bulk_cuts.pop()  # Keep only the top one

    # Update with newly detected info
    rp.LAYER_CUTS = sl.create_layers(rp, bulk_cuts=bulk_cuts)
    rp.N_BULK_LAYERS = len(bulk_cuts)

    # write POSCAR_bulk
    sl.make_bulk_slab(rp)
    try:
        poscar.write(sl.bulkslab, filename='POSCAR_bulk', comments='none')
    except Exception:
        print('Exception occurred while writing POSCAR_bulk')
    else:
        print('Wrote POSCAR_bulk. Check file to see if periodicity is '
              'correct.')

    # Crop off all the atoms at the bottom, leaving only one bulk cell
    newsl = copy.deepcopy(sl)

    newsl.sort_by_z()
    topBulkAt = [at for at in newsl if at.pos[2] <= cut][-1]
    botSlabAt = [at for at in newsl if at.pos[2] > cut][0]
    fracRepeat = np.dot(np.linalg.inv(newsl.ucell), -new_bulk_c)
    newZero = (topBulkAt.pos[2]
               + (abs(botSlabAt.pos[2] - topBulkAt.pos[2]) / 2)
               - abs(fracRepeat[2]))
    newsl.atlist = AtomList(at for at in newsl if at.pos[2] > newZero)
    newsl.update_element_count()   # update the number of atoms per element
    for at in newsl:
        at.pos[2] -= newZero
    newsl.update_cartesian_from_fractional()
    newsl.ucell[:, 2] = newsl.ucell[:, 2] * (1 - newZero)
    newsl.update_fractional_from_cartesian()

    # Finally, convert the cut positions
    bulk_cuts = [round(c - newZero, 3) for c in bulk_cuts]

    # write POSCAR_min
    newsl.sort_original()
    try:
        poscar.write(newsl, filename='POSCAR_min', comments='none')
    except Exception:
        print('Exception occurred while writing POSCAR_min')
    else:
        print('Wrote POSCAR_min, to be used with parameters below.')

    # print info
    print('\nParameters found:')
    print(f'BULK_REPEAT = xyz{rp.BULK_REPEAT}')
    print(f'N_BULK_LAYERS = {len(bulk_cuts)}')
    print(f'LAYER_CUTS = {" ".join(bulk_cuts)}   ! etc.')
    input('Program finished, exit with return')
    return 0


if __name__ == "__main__":
    r = main()
    if r != 0:
        print("Exit with error")
