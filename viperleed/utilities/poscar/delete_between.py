"""ViPErLEED utility: Delete atoms between a certain c fraction."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import argparse
from copy import deepcopy
import logging
import sys

from viperleed.calc.files import poscar
from viperleed.utilities.poscar import add_verbose_option

logger = logging.getLogger(__name__)


def add_cli_parser_arguments(parser):

    parser.add_argument(
        "c",
        help="delete all atoms between these c fractions",
        type=float,
        nargs=2,
    )


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: delete atoms between\n")

    for c in args.c:
        if c <= 0 or c >= 1:
            raise RuntimeError("c must be in range [0, 1]")
    if args.c[0] > args.c[1]:
        raise RuntimeError("First c fraction must be smaller than second c "
                           "fraction.")

    # read the POSCAR file
    slab = poscar.read(sys.stdin)


    # process the slab
    modified_slab = deepcopy(slab)
    modified_slab.atlist = [atom for atom in slab.atlist
                            if atom.pos[2] < args.c[0]
                            or atom.pos[2] > args.c[1]]
    modified_slab.update_element_count()

    logger.debug(f"Deleted {len(slab.atlist) - len(modified_slab.atlist)} atoms"
                 f" in the range c = [{args.c[0]:5.3},{args.c[1]:5.3}].")

    # write the output file
    poscar.write(slab=modified_slab,
                 filename=sys.stdout,
                 comments='none',
                 silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
