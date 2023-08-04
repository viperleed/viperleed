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
from viperleed.utilities.poscar import add_verbose_option

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
    slab = readPOSCAR(sys.stdin)

    # write the output file
    writePOSCAR(slab=slab,
                filename=sys.stdout,
                reorder=True,
                comments='none',
                silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
