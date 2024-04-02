"""ViPErLEED utility: Preparing POSCAR for VASP relaxation

This utility takes a slab in POSCAR format as used by ViPErLEED and prepares it
for relaxation in VASP. This includes adding a vacuum gap on top of the slab,
writing the "Selective dynamics" flags, and giving logical flags for each atom.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import sys

from viperleed.calc.files import poscar
from viperleed.cli_base import float_in_zero_one
from viperleed.utilities.poscar.base import _PoscarStreamCLI

# TODO: add an option to add a mirror image to the slab, so that the slab is
#       symmetric with respect to the center of the slab. This could be useful
#       when dealing with a polar surface.

MIN_VACUUM_GAP = 10  # angstrom. Emit warnings if too small.


class PrepareForVASPRelaxCLI(_PoscarStreamCLI, cli_name='vasp_relax'):
    """Add VASP relaxation tags to the atoms of a POSCAR."""

    long_name = 'preparing slab for VASP relaxation'

    def add_parser_arguments(self, parser):
        """Add arguments specifying which atoms/directions to mark."""
        super().add_parser_arguments(parser)
        parser.add_argument(
            'above_c',
            help='specify above which c fraction to relax the slab.',
            type=float_in_zero_one,
            )
        parser.add_argument(
            '--all_directions',
            help='relax all directions, not just the c direction.',
            action='store_true'
            )

    def process_slab(self, slab, args):
        """Return an unchanged slab, as we only edit comments."""
        logger = self.get_logger()
        logger.debug(f'Relaxing above c fraction {args.above_c}.')
        vacuum = slab.vacuum_gap
        if vacuum < MIN_VACUUM_GAP:
            logger.warning(
                f'Vacuum gap between top and bottom of slab is {vacuum:.2f}A, '
                f'i.e., less than {MIN_VACUUM_GAP}A. This may result in '
                'interaction between the periodic replicas of slabs.'
                )
        return slab

    def write_to_stdout(self, processed_slab, args):
        """Write VASP POSCAR to the terminal."""
        write_vasp_poscar(processed_slab, args)


def write_vasp_poscar(slab, args):
    """Pipe a Slab to stdout given some command-line arguments."""
    # What follows is very similar to poscar.write. The reason not
    # to do this there is to prevent adding a dedicated argument
    # that would only be used in this specific use case. It would
    # also complicate uselessly the code: it would need to decide
    # to use a VASPPOSCARWriter rather than a POSCARFileWriter
    slab.sort_by_element()
    relax_info = {'above_c': args.above_c,
                  'c_only': not args.all_directions}
    writer = poscar.VASPPOSCARWriter(sys.stdout, relax_info=relax_info)
    writer.write(slab)


if __name__ == '__main__':
    PrepareForVASPRelaxCLI.run_as_script()
