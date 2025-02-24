"""POSCAR utility get_bulk_repeat.

Reads a POSCAR file, asks at what c value the bulk starts, then
automatically reduces the size of the POSCAR to non-redundant bulk
layers only, and outputs the appropriate N_BULK_LAYERS and BULK_REPEAT
values.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-12-16'
__license__ = 'GPLv3+'

import copy

import numpy as np

from viperleed.calc.classes.atom_containers import AtomList
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.slab import NoBulkRepeatError
from viperleed.calc.files import poscar
from viperleed.calc.lib.woods_notation import writeWoodsNotation
from viperleed.cli_base import ViPErLEEDCLI


_BULK_DIST_ENOUGH = 1.2
_BULK_DIST_SMALL = 0.7


def ask_user_bulk_cut():
    """Return a valid fractional cut position after asking the user."""
    cut = -1.0
    while cut < 0:
        cutstr = input('Enter c value at which bulk starts: ')
        if not cutstr:
            continue
        try:
            cut = float(cutstr)
        except ValueError:
            print('Failed to convert input to float.')
            continue
        if not 0 <= cut < 1:
            print('Error: c value has to be between 0 and 1.')
            cut = -1.0
    return cut



def ask_user_confirmation_to_cut(bulk_interlayer, bulk_cuts):
    """Ask confirmation if cutting at a somewhat narrow layer spacing."""
    if not _BULK_DIST_SMALL < bulk_interlayer < _BULK_DIST_ENOUGH:
        return

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


def ask_user_symmetry_eps():
    """Return tolerance for symmetry search after asking the user."""
    eps = -1.0
    while eps < 0:
        reply = input('Enter tolerances for symmetry search '
                      '(Default: [0.1 A]): ')
        if not reply:
            eps = 0.1
            break
        try:
            eps = float(reply)
        except ValueError:
            print('Could not convert value to float.')
            continue
        if eps < 0:
            print('Value must be non-negative.')
    return eps


def print_utility_description():
    """Print some introductory information about this utility."""
    print('This utility reads a POSCAR file with arbitrary thickness and '
          'orientation of the bulk, requiring only an input where the bulk '
          'starts. The bulk repeat unit is then automatically determined. '
          'Output is a new POSCAR with a minimum amount of bulk, and values '
          'for the BULK_REPEAT, N_BULK_LAYERS and LAYER_CUTS parameters.\n')


def read_user_poscar():
    """Return a Slab read from a POSCAR file after asking the user."""
    # read the POSCAR file
    filename = ''
    while not filename:
        filename = input('Enter POSCAR file name (Default: [POSCAR]): ')
        if not filename:
            filename = 'POSCAR'
        try:
            slab = poscar.read(filename=filename)
        except FileNotFoundError:
            print(f'File {filename} not found.')
            filename = ''
    print('Slab POSCAR was read successfully.')
    return slab


def write_poscar_min(slab, rpars, cut, bulk_cuts):
    """Write to file a POSCAR with one bulk cell only at the bottom."""
    # Crop off all the atoms at the bottom, leaving only one bulk cell
    slab = copy.deepcopy(slab)
    frac_atoms_z = [at.pos[2] for at in slab]
    frac_topmost_bulk = max(f for f in frac_atoms_z if f <= cut)
    frac_bottom_slab = min(f for f in frac_atoms_z if f > cut)
    dist_to_slab = abs(frac_bottom_slab - frac_topmost_bulk)
    frac_bulk_c = abs(np.linalg.inv(slab.ucell).dot(rpars.BULK_REPEAT)[2])

    new_zero = frac_topmost_bulk + dist_to_slab / 2 - frac_bulk_c
    slab.atlist = AtomList(at for at in slab if at.pos[2] > new_zero)
    slab.update_element_count()   # update the number of atoms per element

    # Shift atoms down, then scale the unit-cell height
    for atom in slab:
        atom.pos[2] -= new_zero
    slab.update_cartesian_from_fractional()
    slab.c_vector[:] *= 1 - new_zero
    slab.update_fractional_from_cartesian()

    # Finally, convert the cut positions
    bulk_cuts = [round(c - new_zero, 3) for c in bulk_cuts]

    # write POSCAR_min
    slab.sort_original()
    try:
        poscar.write(slab, filename='POSCAR_min', comments='none')
    except OSError:
        print('Exception occurred while writing POSCAR_min')
    else:
        print('Wrote POSCAR_min, to be used with parameters below.')

    # print info
    print('\nParameters found:')
    print(f'BULK_REPEAT = xyz{rpars.BULK_REPEAT}')
    print(f'N_BULK_LAYERS = {len(bulk_cuts)}')
    print(f'LAYER_CUTS = {" ".join(bulk_cuts)}   ! etc.')
    input('Program finished, exit with return')


class DetectBulkCLI(ViPErLEEDCLI,
                    cli_name='get_bulk_repeat',
                    help_='automatically detect the bulk portion of a POSCAR'):
    """Detect bulk of a user-given POSCAR from a BULK_LIKE_BELOW."""

    def __call__(self, _=None):
        """Execute the core functionality of the utility.

        Returns
        -------
        exit_code : int
            0: successful completion, 1: errors
        """
        print_utility_description()
        try:
            slab = read_user_poscar()
        except OSError:
            print('Exception while reading POSCAR file')
            return 1

        rpars = Rparams()
        rpars.BULK_LIKE_BELOW = cut = ask_user_bulk_cut()
        rpars.SYMMETRY_EPS = rpars.SYMMETRY_EPS.from_value(
            ask_user_symmetry_eps()
            )

        print('Checking bulk unit cell...')
        try:
            bulk_cuts, bulk_interlayer = slab.detect_bulk(rpars,
                                                          _BULK_DIST_SMALL)
        except NoBulkRepeatError:
            _no_repeat = True
        else:
            _no_repeat = False

        superlattice = rpars.SUPERLATTICE
        ws = writeWoodsNotation(superlattice)                                   # TODO: replace writeWoodsNotation with guilib functions
        info = (
            f'= {ws}' if ws
            else 'M = {} {}, {} {}'.format(*superlattice.astype(int).ravel())
            )
        print(f'Found SUPERLATTICE {info}')

        if _no_repeat:
            print('No repeat vector was found inside the bulk. The bulk may '
                  'already be minimal.')
            return 0

        ask_user_confirmation_to_cut(bulk_interlayer, bulk_cuts)

        # write POSCAR_bulk
        try:
            poscar.write(slab.bulkslab,
                         filename='POSCAR_bulk',
                         comments='none')
        except OSError:
            print('Exception occurred while writing POSCAR_bulk')
        else:
            print('Wrote POSCAR_bulk. Check file to see if periodicity is '
                  'correct.')

        write_poscar_min(slab, rpars, cut, bulk_cuts)
        return 0


if __name__ == '__main__':
    DetectBulkCLI.run_as_script()
