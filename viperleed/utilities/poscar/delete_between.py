#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Delete atoms between a certain c fraction

Created on 2023-08-03

@author: Alexander M. Imre
based on work by Florian Kraushofer
"""

from copy import deepcopy
import logging
import sys


from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import default_cli_parser


logger = logging.getLogger("viperleed.utilities.poscar.delete_between")

def _parse_command_line_arguments():
    parser = default_cli_parser()
    parser.add_argument(
        "c",
        help="delete all atoms between these c fractions",
        type=float,
        nargs=2,
    )
    args, _ = parser.parse_known_args()
    return args


def main():
    args = _parse_command_line_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: delete atoms between\n")

    for c in args.c:
        if c <= 0 or c >= 1:
            raise RuntimeError("c must be in range [0, 1]")
    if args.c[0] > args.c[1]:
        raise RuntimeError("First c fraction must be smaller than second c "
                           "fraction.")

    # read the POSCAR file
    slab = readPOSCAR(sys.stdin)


    # process the slab
    modified_slab = deepcopy(slab)
    modified_slab.atlist = [atom for atom in slab.atlist
                            if atom.pos[2] < args.c[0]
                            or atom.pos[2] > args.c[1]]

    logger.debug(f"Deleted {len(slab.atlist) - len(modified_slab.atlist)} atoms"
                 f" in the range c = [{args.c[0]:5.3},{args.c[1]:5.3}].")

    # write the output file
    writePOSCAR(slab=modified_slab,
                filename=sys.stdout,
                comments='none',
                silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    # if executed from the terminal, send all logs to stderr because stdout is
    # used for piping out the POSCAR file
    logger.addHandler(logging.StreamHandler(sys.stderr))
    main()
