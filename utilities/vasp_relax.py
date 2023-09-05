#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Preparing POSCAR for VASP relaxation

Created on 2023-08-02

@author: Alexander M. Imre
@author: Michele Riva

This utility takes a slab in POSCAR format as used by ViPErLEED and prepares it
for relaxation in VASP. This includes adding a vacuum gap on top of the slab,
writing the "Selective dynamics" flags, and giving logical flags for each atom.
"""

import argparse
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

# TODO: add an option to add a mirror image to the slab, so that the slab is
#       symmetric with respect to the center of the slab. This could be useful
#       when dealing with a polar surface.

logger = logging.getLogger(
    "viperleed.utilities.poscar.prepare_for_vasp_relaxation"
    )


def _parse_command_line_arguments():
    parser = default_cli_parser()
    parser.add_argument(
        "above_c",
        help=("Specify above which c fraction to relax the slab."),
        type=float,
    )
    parser.add_argument(
        "--reorder",
        help=("Reorder the atoms in the slab by z coordinate."
              "Atoms will always be sorted by original-element order."),
        action="store_true"
    )
    parser.add_argument(
        "--all_directions",
        help=("Relax all directions, not just the c direction."),
        action="store_true"
    )
    args, _ = parser.parse_known_args()
    return args


def write_vasp_poscar(slab, args):
    """Pipe a Slab to stdout given some command-line arguments."""
    # What follows is very similar to poscar.write. The reason not
    # to do this there is to prevent adding a dedicated argument
    # that would only be used in this specific use case. It would
    # also complicate uselessly the code: it would need to decide
    # to use a VASPPOSCARWriter rather than a POSCARFileWriter
    if args.reorder:
        slab.sort_by_z()
    slab.sort_by_element()
    relax_info = {'above_c': args.above_c,
                  'c_only': not args.all_directions}
    writer = poscar.VASPPOSCARWriter(sys.stdout, relax_info=relax_info)
    writer.write(slab)


def main():
    args = _parse_command_line_arguments()
    if args.verbose:
        logger.setLevel(logging.DEBUG)

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
    write_vasp_poscar(slab, args)


if __name__ == "__main__":
    # if executed from the terminal, send all logs to stderr because stdout is
    # used for piping out the POSCAR file
    logger.addHandler(logging.StreamHandler(sys.stderr))
    main()
