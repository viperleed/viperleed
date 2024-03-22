"""ViPErLEED utility: Delete atoms above a certain height."""

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
        help="delete all atoms above this c fraction",
        type=float,
    )


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    logger.info("ViPErLEED utility: delete atoms above\n")

    if args.c <= 0 or args.c >= 1:
        raise RuntimeError("c must be in range [0, 1]")

    # read the POSCAR file
    slab = poscar.read(sys.stdin)


    # process the slab
    modified_slab = deepcopy(slab)
    modified_slab.atlist = [atom for atom in slab.atlist
                            if args.c > atom.pos[2]]
    modified_slab.update_element_count()

    logger.debug(f"Deleted {len(slab.atlist) - len(modified_slab.atlist)} atoms"
                 f" above c > {args.c}.")

    # write the output file
    poscar.write(slab=modified_slab,
                 filename=sys.stdout,
                 comments='none',
                 silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
