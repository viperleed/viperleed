"""ViPErLEED utility: Merge POSCAR files."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import argparse
from copy import deepcopy
import logging
import sys

import numpy as np
from scipy.spatial.distance import cdist

from viperleed.calc.files import poscar
from viperleed.utilities.poscar import add_verbose_option

logger = logging.getLogger(__name__)


def merge_slabs(slabs, check_collisions=True, eps=1e-5):
    merged_slab = deepcopy(slabs[0])
    merged_slab.atlist.strict = False  # won't allow duplicate atom numbers otherwise
    for slab in slabs[1:]:
        if check_collisions:
            # check for collisions
            dists = cdist([at.cartpos for at in merged_slab],
                          [at.cartpos for at in slab])
            if np.any(dists < eps):
                raise RuntimeError("Distance between atoms in different slabs "
                                   "is smaller than eps. ")
        merged_slab.atlist.extend(slab.atlist)

    merged_slab.update_atom_numbers()
    merged_slab.update_element_count()
    merged_slab.update_cartesian_from_fractional()
    merged_slab.atlist.strict = True

    return merged_slab

def add_cli_parser_arguments(parser):

    parser.add_argument(
        "files",
        nargs="+",
    )
    parser.add_argument(
        "--eps-cell",
        help=("maximum element-wise allowed difference in lattice parameters "
              "between the unit cells of the POSCAR files to be merged"),
        type=float,
        default=1e-3,
    )
    parser.add_argument(
        "--eps-collision",
        help=("minimum distance between atoms in different before a collision "
              "is detected"),
        type=float,
        default=1e-1
    )


def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: merge POSCAR files\n")


    slabs = []
    for filename in args.files:
        slabs.append(poscar.read(filename))

    if len(slabs) < 2:
        raise RuntimeError("Need at least two POSCAR files to merge.")

    # check that all slabs have the same lattice parameters
    for slab in slabs[1:]:
        if not np.allclose(slab.ucell, slabs[0].ucell, atol=args.eps_cell):
            raise RuntimeError("Files do not have the same unit cell. Use "
                               "--eps to increase tolerance.")

    # merge the slabs
    merged_slab = merge_slabs(slabs,
                              check_collisions=True,
                              eps=args.eps_collision)

    # write the output file
    poscar.write(slab=merged_slab,
                 filename=sys.stdout,
                 comments='none',
                 silent=logger.level<=logging.DEBUG)


if __name__ == "__main__":
    main()
