#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Rescale cell

This utility takes a slab in POSCAR format and rescales the unit cell.


Created on 2023-08-03

@author: Alexander M. Imre
"""
import argparse
from copy import deepcopy
import logging
import sys

from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import add_verbose_option

logger = logging.getLogger("viperleed.utilities.poscar.prepare_for_vasp_relaxation")


def add_cli_parser_arguments(parser):

    parser.add_argument(
        "scaling",
        help=("One or three scaling factors for the unit cell. If three values "
              "are given, the scaling factors are applied to the a, b, and c "
              "vector, respectively. If only one value is given, an isotropic "
              "scaling is applied."),
        type=float,
        nargs="+",
    )


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: rescale unit cell\n")

    # read the POSCAR file
    slab = readPOSCAR(sys.stdin)

    print(args.scaling)
    if len(args.scaling) == 1 or len(args.scaling) == 3:
        slab.apply_scaling(*args.scaling)
    else:
        raise ValueError("The number of scaling factors must be either 1 or 3.")

    # write the output file
    writePOSCAR(slab=slab,
                filename=sys.stdout,
                comments='none',
                silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
