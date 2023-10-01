#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Preparing POSCAR for VASP relaxation

This utility takes a slab in POSCAR format as used by ViPErLEED and prepares it 
for relaxation in VASP. This includes adding a vacuum gap on top of the slab,
writing the "Selective dynamics" flags, and giving logical flags for each atom.

Created on 2023-08-03

@author: Alexander M. Imre
"""

from copy import deepcopy
import logging
import sys
import os

cd = os.path.realpath(os.path.dirname(__file__))
# NB: it's necessary to add vpr_path to sys.path so that viperleed
#     can be loaded correctly at the top-level package
vpr_path = os.path.realpath(os.path.join(cd, '..', '..'))
for import_path in (cd, vpr_path):
    if import_path not in sys.path:
        sys.path.append(import_path)

from viperleed.tleedmlib.files import poscar
from viperleed.utilities import default_cli_parser


logger = logging.getLogger(
    "viperleed.utilities.poscar.prepare_for_vasp_relaxation"
    )


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

    processed_slab.ucell[:, 2] = new_c_z / processed_slab.ucell[2, 2]
    processed_slab.ucell[2, 2] = new_c_z
    processed_slab.getFractionalCoordinates()
    processed_slab.collapseFractionalCoordinates()
    return processed_slab


def _parse_command_line_arguments():
    parser = default_cli_parser()
    parser.add_argument(
        "vacuum",
        help=("Add vacuum on top of the slab. Value in Å. Default is 0 Å."),
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "-a", "--absolute",
        help=("If set, the value given to --vacuum is the absolute size"
              "of the vacuum gap, else it is the amount of vacuum to add or"
              "remove from the slab."),
        action = "store_true"
    )
    args, _ = parser.parse_known_args()
    return args


def main():
    args = _parse_command_line_arguments()

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
    # if executed from the terminal, send all logs to stderr because stdout is
    # used for piping out the POSCAR file
    logger.addHandler(logging.StreamHandler(sys.stderr))
    main()
