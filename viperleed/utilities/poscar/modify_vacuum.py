"""ViPErLEED utility: Modify vacuum gap.

This utility takes a slab in POSCAR format and modifies the vacuum gap.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import argparse
from copy import deepcopy
import logging
import sys

from viperleed.calc.files import poscar
from viperleed.utilities.poscar import add_verbose_option


logger = logging.getLogger(__name__)


def modify_vacuum(slab, vacuum_gap_size, absolute=False):
    """Modify the vacuum gap size of a slab.

    Parameters
    ----------
    slab : Slab
        Slab to prepare for relaxation,
    vacuum_gap_size : float
        Size of the vacuum gap in Angstroms. If absolute is True, this
        is the absolute size of the vacuum gap, else it is the amount
        of vacuum to add or remove from the slab.
    absolute : bool, optional
        If True, vacuum_gap_size is the absolute size of the vacuum gap.

    Returns
    -------
    Slab
        The slab with the modified vacuum gap.

    Raises
    ------
    RuntimeError
        If the resulting vacuum gap size is negative.
    """
    processed_slab = deepcopy(slab)
    processed_slab.check_a_b_in_plane()
    processed_slab.update_cartesian_from_fractional()

    slab_thickness = processed_slab.thickness
    current_gap_size = processed_slab.vacuum_gap

    if absolute:
        new_vacuum_gap_size = vacuum_gap_size
    else:
        new_vacuum_gap_size = current_gap_size + vacuum_gap_size

    logger.debug(f"Current vacuum gap size:\t{current_gap_size:9.3f}")
    logger.debug(f"New vacuum gap size:\t\t{new_vacuum_gap_size:9.3f}")

    if new_vacuum_gap_size < 0:
        raise RuntimeError("The resulting vacuum gap size would be "
                            "negative.")

    new_c_z = new_vacuum_gap_size + slab_thickness

    processed_slab.ucell[:, 2] = new_c_z / processed_slab.c_vector[2]           # TODO: @amimre is this (and the next line) correct??
    processed_slab.ucell[2, 2] = new_c_z
    processed_slab.collapse_cartesian_coordinates()
    return processed_slab


def add_cli_parser_arguments(parser):

    parser.add_argument(
        "vacuum",
        help=("Add vacuum on top of the slab. Value in Å. Default is 0 Å."),
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "-a", "--absolute",
        help=("If set, the value given to vacuum is the absolute size"
              "of the vacuum gap, else it is the amount of vacuum to add or"
              "remove from the slab."),
        action = "store_true"
    )



def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser()
        add_verbose_option(parser)
        add_cli_parser_arguments(parser)
        args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: modify vacuum gap\n")

    if args.absolute and args.vacuum < 0:
        raise ValueError("Absolute vacuum gap has to be positive, but was "
                         f"{args.vacuum}.")

    if args.absolute:
        logger.debug("Using absolute vacuum gap size.")

    # read the POSCAR file
    slab = poscar.read(sys.stdin)

    # process the slab
    processed_slab = modify_vacuum(slab, args.vacuum, absolute=args.absolute)

    # write the output file
    poscar.write(slab=processed_slab,
                 filename=sys.stdout,
                 comments='none',
                 silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
