#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Sort slab by z

Created on 2023-08-03

@author: Alexander M. Imre
"""

import argparse
import logging
import sys
import os

from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import default_cli_parser

logger = logging.getLogger("viperleed.utilities.poscar.sort_by_z")


def _parse_command_line_arguments():
    parser, args, unparsed_args = default_cli_parser()
    parser.parse_args(args=unparsed_args, namespace=args)
    return args


def main():
    args = _parse_command_line_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.debug("ViPErLEED utility: sort slab by z\n")

    # read the POSCAR file from stdin
    slab = readPOSCAR(sys.stdin)

    # write the output file
    writePOSCAR(slab=slab,
                filename=sys.stdout,
                reorder=True,
                comments='none',
                silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
