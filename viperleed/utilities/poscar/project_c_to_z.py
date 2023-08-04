#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Project c vector to z axis

This utility takes a POSCAR slab and projects the c vector to the z axis.

Created on 2023-08-03

@author: Alexander M. Imre
based on work by Florian Kraushofer
"""
import argparse
import logging
import os
import sys

from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR

logger = logging.getLogger("viperleed.utilities.poscar.project_c_to_z")

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

    logger.debug("ViPErLEED utility: Project c to z axis\n")

    # read the POSCAR file from stdin
    slab = readPOSCAR(sys.stdin)

    # project the c vector to the z axis
    slab.projectCToZ()

    # write the output file
    writePOSCAR(slab=slab,
                filename=sys.stdout,
                comments='none',
                silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
