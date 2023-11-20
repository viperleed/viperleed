#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Sort slab by z
"""
import argparse
import logging
import sys
import os

from viperleed.calc.files import poscar
from viperleed.utilities.poscar import add_verbose_option

__authors__ = ["Alexander M. Imre (@amimre)",]
__created__ = "2023-08-03"

logger = logging.getLogger("viperleed.utilities.poscar.sort_by_z")


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

    logger.debug("ViPErLEED utility: sort slab by z\n")

    # read the POSCAR file from stdin
    slab = poscar.read(sys.stdin)

    # write the output file
    poscar.write(slab=slab,
                 filename=sys.stdout,
                 reorder=True,
                 comments='none',
                 silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
