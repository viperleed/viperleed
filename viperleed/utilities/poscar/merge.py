#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Merge POSCAR files

Created on 2023-08-03

@author: Alexander M. Imre
based on work by Florian Kraushofer
"""

from copy import deepcopy
import logging
import sys

import numpy as np
from scipy.spatial.distance import cdist

from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import default_cli_parser


logger = logging.getLogger("viperleed.utilities.poscar.merge")


def merge_slabs(slabs, check_collisions=True, eps=1e-5):
    merged_slab = deepcopy(slabs[0])
    for slab in slabs[1:]:
        if check_collisions:
            # check for collisions
            dists = cdist(merged_slab.cartpos, slab.cartpos)
            if np.any(dists < eps):
                raise RuntimeError("Distance between atoms in different slabs "
                                   "is smaller than eps. ")
        merged_slab.atlist.extend(slab.atlist)

    merged_slab.updateAtomNumbers()
    merged_slab.updateElementCount()
    merged_slab.getCartesianCoordinates()

    return merged_slab

def _parse_command_line_arguments():
    parser, args, unparsed_args = default_cli_parser()
    parser.add_argument(
        "files",
        nargs="+",
    )
    parser.add_argument(
        "--eps",
        help=("maximum element-wise allowed difference in lattice parameters "
              "between the unit cells of the POSCAR files to be merged"),
        type=float,
        default=1e-6,
    )
    parser.add_argument(
        "--no_check_collisions",
        help="do not check for collisions between atoms in different slabs",
        action="store_true",
    )
    parser.parse_args(args=unparsed_args, namespace=args)
    return args


def main():
    args = _parse_command_line_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: merge POSCAR files\n")


    slabs = []
    for filename in args.files:
        slabs.append(readPOSCAR(filename))

    if len(slabs) < 2:
        raise RuntimeError("Need at least two POSCAR files to merge.")

    # check that all slabs have the same lattice parameters
    for slab in slabs[1:]:
        if not np.allclose(slab.ucell, slabs[0].ucell, atol=args.eps):
            raise RuntimeError("Files do not have the same unit cell. Use "
                               "--eps to increase tolerance.")

    # merge the slabs
    merged_slab = merge_slabs(slabs,
                              check_collisions=not args.no_check_collisions,
                              eps=args.eps)

    # write the output file
    writePOSCAR(slab=merged_slab,
                filename=sys.stdout,
                comments='none',
                silent=logger.level<=logging.DEBUG)


if __name__ == "__main__":
    main()
