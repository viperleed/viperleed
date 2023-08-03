#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Strip comments from POSCAR file

Created on 2023-08-03

@author: Alexander M. Imre
"""

from copy import deepcopy
import logging
import sys

from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import default_cli_parser


logger = logging.getLogger("viperleed.utilities.poscar.strip_comments")


def _parse_command_line_arguments():
    parser, args, unparsed_args = default_cli_parser()
    parser.parse_args(args=unparsed_args, namespace=args)
    return args


def main():
    args = _parse_command_line_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: strip comments\n")

    # read the POSCAR files
    slab = readPOSCAR(sys.stdin)

    # write the output file without comments
    writePOSCAR(slab=slab,
                filename=sys.stdout,
                comments='none',
                silent=logger.level<=logging.DEBUG)


if __name__ == "__main__":
    main()
