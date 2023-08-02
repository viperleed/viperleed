#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Preparing POSCAR for VASP relaxation

This utility takes a slab in POSCAR format as used by ViPErLEED and prepares it 
for relaxation in VASP. This includes adding a vacuum gap on top of the slab,
writing the "Selective dynamics" flags, and giving logical flags for each atom.

Created on 2023-08-02

@author: Alexander M. Imre
"""

import argparse
import logging
from copy import deepcopy
from pathlib import Path
from warnings import warn

from viperleed.lib.files.poscar import readPOSCAR, writePOSCAR

# TODO: add an option to add a mirror image to the slab, so that the slab is
#       symmetric with respect to the center of the slab. This could be useful
#       when dealing with a polar surface.

logger = logging.getLogger("viperleed.utilities.poscar.prepare_for_vasp_relaxation")


def prepare_for_vasp_relaxation(slab, vacuum_gap_size):
    """_summary_

    Parameters
    ----------
    slab : Slab
        slab to prepare for relaxation
    vacuum_gap_size : float
        size of the vacuum gap to add on top of the slab in Angstroms

    Returns
    -------
    Slab
        the slab with the vacuum gap added on top
    """
    processed_slab = deepcopy(slab)
    processed_slab.check_a_b_out_of_plane()
    processed_slab.getCartesianCoordinates()
    top_atom_z = max([at.cartpos[2] for at in processed_slab.atlist])
    old_c_z = processed_slab.ucell[2, 2]
    new_c_z = top_atom_z + vacuum_gap_size
    processed_slab.ucell[2, 2] = new_c_z
    processed_slab.getFractionalCoordinates()
    processed_slab.collapseFractionalCoordinates()
    return processed_slab


def _parse_command_line_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "above_c",
        help=("Specify above which c fraction to relax the slab."),
        type=float,
    )
    parser.add_argument(
        "-i", "--input",
        help=("specify input file"),
        type=str,
        default="POSCAR"
        )
    parser.add_argument(
        "-o", "--output",
        help=("specify output file"),
        type=str,
        default="POSCAR_for_relaxation"
        )
    parser.add_argument(
        "-v", "--verbose",
        help=("increase output verbosity"),
        action="store_true"
    )
    parser.add_argument(
        "--vacuum",
        help=("Specify how much vacuum to add on top of the slab in Å, "
              "measured from the topmost atom. Default is 20 Å."),
        type=float,
        default=20.0,
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


def main():
    args = _parse_command_line_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: Preparing slab for VASP relaxation\n")

    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()

    logger.debug(f"Input file: {input_path}")
    logger.debug(f"Output file: {output_path}")

    if not input_path.is_file():
        raise FileNotFoundError(f"Input file {input_path} not found.")
    if output_path.is_file():
        raise FileExistsError(f"Output file {output_path} already exists.")

    above_c = args.above_c
    if above_c <= 0 or above_c >= 1:
        raise ValueError("c fraction has to be in range [0,1], but was "
                         f"{above_c}.")

    if args.vacuum <= 0:
        raise ValueError("Vacuum gap has to be positive, but was "
                         f"{args.vacuum}.")
    if args.vacuum < 15:
        logger.warn("Vacuum gap is less than 15 Å. This may lead to "
                    "interactions between the slab and its periodic images. "
                    "Consider increasing the vacuum gap.")

    c_only = not args.all_directions
    if not c_only:
        logger.warn("Relaxing all directions may break the slab symmetry. "
                    "Use with caution.")

    # read the POSCAR file
    slab = readPOSCAR(input_path)

    # TODO: process the slab
    processed_slab = prepare_for_vasp_relaxation(slab, args.vacuum)


    # TODO: write the output file
    writePOSCAR(slab=processed_slab,
                filename=output_path,
                reorder=args.reorder,
                comments='relax',
                relax_info={"above_c": above_c,
                            "c_only": c_only},
                silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    logger.addHandler(logging.StreamHandler())
    main()
