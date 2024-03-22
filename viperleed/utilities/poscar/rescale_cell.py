"""ViPErLEED utility: Rescale cell.

This utility takes a slab in POSCAR format and rescales the unit cell.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
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
        "scaling",
        help=("One or three scaling factors for the unit cell. If three values "
              "are given, the scaling factors are applied to the a, b, and c "
              "vector, respectively. If only one value is given, an isotropic "
              "scaling is applied."),
        type=float,
        nargs="+",
    )


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: rescale unit cell\n")

    # read the POSCAR file
    slab = poscar.read(sys.stdin)

    if len(args.scaling) == 1 or len(args.scaling) == 3:
        slab.apply_scaling(*args.scaling)
    else:
        raise ValueError("The number of scaling factors must be either 1 or 3.")

    # write the output file
    poscar.write(slab=slab,
                 filename=sys.stdout,
                 comments='none',
                 silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
