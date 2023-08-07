#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Preparing POSCAR for VASP relaxation

This utility takes a slab in POSCAR format as used by ViPErLEED and prepares it
for relaxation in VASP. This includes writing the "Selective dynamics" flags
and giving logical flags for each atom that specify along which direction to
relax.

Created on 2023-08-02

@author: Alexander M. Imre
"""
import argparse
import logging
import sys
import os

from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import add_verbose_option


# TODO: add an option to add a mirror image to the slab, so that the slab is
#       symmetric with respect to the center of the slab. This could be useful
#       when dealing with a polar surface.

logger = logging.getLogger("viperleed.utilities.poscar.prepare_for_vasp_relaxation")


def add_cli_parser_arguments(parser):

    parser.add_argument(
        "above_c",
        help=("Specify above which c fraction to relax the slab."),
        type=float,
    )
    parser.add_argument(
        "--all_directions",
        help=("Relax all directions, not just the c direction."),
        action="store_true"
    )



def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    logger.debug("ViPErLEED utility: Preparing slab for VASP relaxation\n")

    above_c = args.above_c
    if above_c <= 0 or above_c >= 1:
        raise ValueError("c fraction has to be in range [0,1], but was "
                         f"{above_c}.")
    logger.debug(f"Relaxing above c fraction {above_c}.")

    # read the POSCAR file from stdin
    slab = readPOSCAR(sys.stdin)

    if slab.vacuum_gap < 10:
        logger.warning("Gap between top and bottom of slab is less than 10 Å. "
                       "This may result in interaction between the slabs.")

    c_only = not args.all_directions

    # write the output file
    writePOSCAR(slab=slab,
                filename=sys.stdout,
                comments='relax',
                relax_info={"above_c": above_c,
                            "c_only": c_only},
                silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
