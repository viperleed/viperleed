#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Strip comments from POSCAR file
"""
import argparse
from copy import deepcopy
import logging
import sys

from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import add_verbose_option

__authors__ = ["Alexander M. Imre (@amimre)",]
__created__ = "2023-08-03"

logger = logging.getLogger("viperleed.utilities.poscar.strip_comments")


def add_cli_parser_arguments(parser):
    pass


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.debug("ViPErLEED utility: strip comments\n")

    # read the POSCAR files
    slab = readPOSCAR(sys.stdin)

    # write the output file without comments
    writePOSCAR(slab=slab,
                filename=sys.stdout,
                comments='none',
                silent=logger.level<=logging.DEBUG)


if __name__ == "__main__":
    main()
