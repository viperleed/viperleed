"""ViPErLEED utility: Reorder elements."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import argparse
import logging
import sys
import os

import numpy as np

from viperleed.calc.files import poscar
from viperleed.calc.lib import periodic_table
from viperleed.utilities.poscar import add_verbose_option

logger = logging.getLogger(__name__)


def add_cli_parser_arguments(parser):

    parser.add_argument(
        "-a", "--alphabetical",
        help="Sort elements alphabetically",
        action="store_true",
    )
    parser.add_argument(
        "-d", "--descending",
        help="Sort elements by atomic number in descending order",
        action="store_true",
    )
    parser.add_argument(
        "-c", "--custom",
        help=("Sort elements by custom order. Provide a comma-separated list of "
             "elements in the desired order."),
        type=str,
    )


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.debug("ViPErLEED utility: reorder elements\n")

    if sum([args.alphabetical, args.descending, bool(args.custom)]) > 1:
        raise ValueError("Only one of the options -a, -d, and -c can be used.")

    # read the POSCAR file from stdin
    slab = poscar.read(sys.stdin)

    # sort the elements as desired
    if args.custom:
        _custom_order = [elem.lower() for elem in args.custom.split(",")]
        for elem in slab.elements:
            if elem.lower() not in _custom_order:
                raise ValueError(f"Element {elem} not found in custom order.")
        elem_key = lambda elem: _custom_order.index(elem.lower())
    elif args.alphabetical:
        elem_key = lambda elem: elem.lower()
    elif args.descending:
        elem_key = lambda elem: periodic_table.get_atomic_number(elem)*-1
    else:
        elem_key = lambda elem: periodic_table.get_atomic_number(elem)

    key = lambda elem_item: elem_key(elem_item[0])

    # reorder the elements
    slab.n_per_elem = dict(sorted(slab.n_per_elem.items(), key=key))
    slab.updateElementCount()

    # write the output file
    slab.sort_by_z()
    poscar.write(slab=slab,
                 filename=sys.stdout,
                 comments='none',
                 silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
