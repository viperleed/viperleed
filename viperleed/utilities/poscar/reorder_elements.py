#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Reorder elements

Created on 2023-08-03

@author: Alexander M. Imre
"""

import argparse
import logging
import sys
import os

import numpy as np

import viperleed.calc.periodic_table as periodic_table
from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import default_cli_parser

logger = logging.getLogger("viperleed.utilities.poscar.reorder_elements")


def _parse_command_line_arguments():
    parser, args, unparsed_args = default_cli_parser()
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
    parser.parse_args(args=unparsed_args, namespace=args)
    return args


def main():
    args = _parse_command_line_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.debug("ViPErLEED utility: reorder elements\n")

    if sum([args.alphabetical, args.descending, bool(args.custom)]) > 1:
        raise ValueError("Only one of the options -a, -d, and -c can be used.")

    # read the POSCAR file from stdin
    slab = readPOSCAR(sys.stdin)

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
    writePOSCAR(slab=slab,
                filename=sys.stdout,
                reorder=True,
                comments='none',
                silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
