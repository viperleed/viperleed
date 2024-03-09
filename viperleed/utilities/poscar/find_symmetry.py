#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Find Symmetry
"""
import argparse
from copy import deepcopy
import logging
import sys

from viperleed.calc.lib import symmetry
from viperleed.calc.classes import rparams
from viperleed.calc.files import poscar
from viperleed.utilities.poscar import add_verbose_option

__authors__ = ["Alexander M. Imre (@amimre)",]
__created__ = "2023-08-03"

logger = logging.getLogger("viperleed.utilities.poscar.prepare_for_vasp_relaxation")


def add_cli_parser_arguments(parser):

    parser.add_argument(
        "-e", "--symmetry-eps",
        help=("Epsilon for symmetry detection in Å. Default: 0.1Å"),
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--symmetry-eps-z",
        help=("Epsilon for symmetry detection in z in Å. If not provided, "
              "the value of --symmetry-eps is used."),
        type=float,
    )


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: find symmetry\n")

    # read the POSCAR file
    slab = poscar.read(sys.stdin)

    param = rparams.Rparams()
    slab.full_update(param)
    param.SYMMETRY_EPS = args.symmetry_eps
    param.SYMMETRY_EPS_Z = (args.symmetry_eps_z
                            if args.symmetry_eps_z is not None
                            else args.symmetry_eps)
    param.SYMMETRY_FIND_ORI = True

    # find the symmetry
    found_symmetry = symmetry.findSymmetry(slab, param)

    # write the found symmetry group
    sys.stdout.write(found_symmetry + "\n")

if __name__ == "__main__":
    main()
