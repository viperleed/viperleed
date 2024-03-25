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

import argparse
import logging
import sys

from viperleed.calc.files import poscar
from viperleed.utilities.poscar import add_verbose_option

# TODO: add an option to add a mirror image to the slab, so that the slab is
#       symmetric with respect to the center of the slab. This could be useful
#       when dealing with a polar surface.


logger = logging.getLogger(
    "viperleed.utilities.poscar.prepare_for_vasp_relaxation"
    )


def add_cli_parser_arguments(parser):

    parser.add_argument(
        "above_c",
        help=("Specify above which c fraction to relax the slab."),
        type=float,
    )
    parser.add_argument(
        "--all_directions",
        help=("Relax all directions, not just the c direction."),
        action="store_true"
    )


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


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)


    logger.debug("ViPErLEED utility: Preparing slab for VASP relaxation\n")

    above_c = args.above_c
    if above_c <= 0 or above_c >= 1:
        raise ValueError("c fraction has to be in range [0,1], but was "
                         f"{above_c}.")
    logger.debug(f"Relaxing above c fraction {above_c}.")

    # read the POSCAR file from stdin
    slab = poscar.read(sys.stdin)

    if slab.vacuum_gap < 10:
        logger.warning("Gap between top and bottom of slab is less than 10 Å. "
                       "This may result in interaction between the slabs.")
    write_vasp_poscar(slab, args)


if __name__ == "__main__":
    main()
