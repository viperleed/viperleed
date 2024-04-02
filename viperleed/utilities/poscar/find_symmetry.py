"""ViPErLEED utility: Find Symmetry."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import sys

from viperleed.calc import symmetry
from viperleed.utilities.poscar.base import _PoscarSymmetryCLI


class FindSymmetryCLI(_PoscarSymmetryCLI, cli_name='find_symmetry'):
    """Detect the plane group of a the structure in a POSCAR."""

    long_name = 'find symmetry'

    def process_slab(self, slab, args):
        """Find the plane group of slab."""
        rpars = self.prepare_rpars(slab, args)
        symmetry.findSymmetry(slab, rpars)

    def write_to_stdout(self, processed_slab, _):
        """Write the detected plane group to the terminal."""
        sys.stdout.write(f'{processed_slab.foundplanegroup}\n')


if __name__ == '__main__':
    FindSymmetryCLI.run_as_script()
