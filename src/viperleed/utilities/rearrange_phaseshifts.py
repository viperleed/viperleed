"""ViPErLEED utility: rearrange phaseshifts.

Reads a PHASESHIFTS file and allows the
user to copy and rearrange the blocks.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2019-12-16'
__license__ = 'GPLv3+'


from viperleed.calc.files.phaseshifts import (readPHASESHIFTS,
                                              writePHASESHIFTS)
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.cli_base import ViPErLEEDCLI


_HEADER = '''
This utility reads a phase-shifts file (e.g., PHASESHIFTS) and
allows the user to copy and rearrange the blocks, such that they
fit the order of elements in POSCAR and the order of sites as
defined by ELEMENT_MIX and SITE_DEF in the PARAMETERS file.
Check the documentation for further information on how to arrange
blocks in the PHASESHIFTS file.
'''
_HELP = ('copy and rearrange blocks in a PHASESHIFTS file to fit '
         'POSCAR elements, ELEMENT_MIX, and SITE_DEF')
_INSTRUCTIONS = '''
Found {n_blocks} blocks. Enter a space-separated list of how the blocks
should be arranged in the new file.

Examples:
'1 1 2 2': Print the first block twice, then the second block twice
'2 1': Swap the first an the second block
'2 3': Print the second, then the third block, delete the first block.
'''


class RearrangePhaseShiftsCLI(ViPErLEEDCLI,
                              cli_name='rearrange-phaseshifts',
                              help_=_HELP):
    """Utility to swap around and duplicate blocks in a PHASESHIFTS file."""

    def __call__(self, args=None):
        """Call the phase-shifts rearrangement utility."""
        # Print some info, then read the phase-shifts file
        print(_HEADER)
        try:
            firstline, phaseshifts, *_ = _read_phaseshifts_file()
        except (ValueError, IndexError)  as exc:
            print(f'Exception while reading phaseshifts file: {exc}')
            return 1
        print('Phase-shifts file read successfully.')

        # Print some more instructions
        n_blocks = len(phaseshifts[0][1])
        print(_INSTRUCTIONS.format(n_blocks=n_blocks))

        # Get input for new order, then rearrange blocks
        new_blocks = _ask_user_new_order(n_blocks)
        phaseshifts = _rearrange_blocks(phaseshifts, new_blocks)

        # Update the block count in the first line (e.g.,
        # user replicated some (site, element) pairs)
        n_blocks = len(phaseshifts[0][1])
        firstline = f'{n_blocks:>3d}{firstline[3:]}'

        # write new file
        filename = f'PHASESHIFTS_mod_{DateTimeFormat.FILE_SUFFIX.now()}'
        try:
            writePHASESHIFTS(firstline, phaseshifts, file_path=filename)
        except Exception:
            print('Error writing new phase-shifts file.')
            raise
        print(f'Wrote new phaseshifts file as {filename}')
        return 0


def _ask_user_new_order(n_blocks):
    """Return the new order of blocks after asking the user."""
    while True:
        neworder_str = input('Enter new order: ')
        if not neworder_str:
            print('Input failed. Please try again.')
            continue
        try:
            neworder = [int(s) for s in neworder_str.split()]
        except ValueError:
            print('Could not parse input. Please try again.')
            continue
        if all(0 < i < n_blocks+1 for i in neworder):
            return neworder
        print('Input out of bounds. Please try again.')


def _read_phaseshifts_file():
    """Read a PHASESHIFTS file specified by the user."""
    try:
        return readPHASESHIFTS(None, None, check=False)
    except FileNotFoundError:
        print('PHASESHIFTS file not found.')

    filename = ''
    while not filename:
        filename = input('Enter phase-shifts file name: ')
        if not filename:
            print('Input failed. Please try again.')
            continue
        try:
            return readPHASESHIFTS(None, None, readfile=filename, check=False)
        except FileNotFoundError:
            print(f'{filename} not found.')
            filename = ''
    return None, None, None, None


def _rearrange_blocks(phaseshifts_vs_energy, new_blocks):
    """Return phase-shifts with blocks re-ordered as per new_blocks."""
    new_phaseshifts = []
    for energy, phaseshifts_at_energy in phaseshifts_vs_energy:
        new_phaseshifts_at_energy = [phaseshifts_at_energy[i - 1]
                                     for i in new_blocks]
        new_phaseshifts.append((energy, new_phaseshifts_at_energy))
    return new_phaseshifts


if __name__ == '__main__':
    RearrangePhaseShiftsCLI.run_as_script()
