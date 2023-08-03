#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Find Symmetry



Created on 2023-08-03

@author: Alexander M. Imre
"""

from copy import deepcopy
import logging
import sys


from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.calc import symmetry
from viperleed.calc.classes import rparams
from viperleed.utilities.poscar import default_cli_parser


logger = logging.getLogger("viperleed.utilities.poscar.prepare_for_vasp_relaxation")


def _parse_command_line_arguments():
    parser, args, unparsed_args = default_cli_parser()
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
    parser.parse_args(args=unparsed_args, namespace=args)
    return args


def main():
    args = _parse_command_line_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: find symmetry\n")

    # read the POSCAR file
    slab = readPOSCAR(sys.stdin)

    param = rparams.Rparams()
    slab.fullUpdate(param)
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
