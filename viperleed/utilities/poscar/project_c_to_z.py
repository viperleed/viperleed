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
import sys
import os

from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import default_cli_parser


logger = logging.getLogger("viperleed.utilities.poscar.project_c_to_z")


def _parse_command_line_arguments():
    parser = default_cli_parser()
    args, _ = parser.parse_known_args()
    return args


def main():
    args = _parse_command_line_arguments()

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
    # if executed from the terminal, send all logs to stderr because stdout is
    # used for piping out the POSCAR file
    logger.addHandler(logging.StreamHandler(sys.stderr))
    main()
