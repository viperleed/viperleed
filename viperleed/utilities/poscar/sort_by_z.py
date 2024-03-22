"""ViPErLEED utility: Sort slab by z."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import argparse
import logging
import sys

from viperleed.calc.files import poscar
from viperleed.utilities.poscar import add_verbose_option

logger = logging.getLogger(__name__)


def add_cli_parser_arguments(parser):
    parser.add_argument(
        "-r", "--reversed",
        help="Reverse z ordering (bottom to top)",
        action="store_true",
    )


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.debug("ViPErLEED utility: sort slab by z\n")

    # read the POSCAR file from stdin
    slab = poscar.read(sys.stdin)
    slab.sort_by_z(bottom_to_top=args.reversed)

    # write the output file
    poscar.write(slab=slab,
                 filename=sys.stdout,
                 comments='none',
                 silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
