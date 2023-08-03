#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utility: Delete atoms below a certain height

Created on 2023-08-03

@author: Alexander M. Imre
based on work by Florian Kraushofer
"""

from copy import deepcopy
import logging
import sys


from viperleed.calc.files.poscar import readPOSCAR, writePOSCAR
from viperleed.utilities.poscar import default_cli_parser


logger = logging.getLogger("viperleed.utilities.poscar.prepare_for_vasp_relaxation")

def _parse_command_line_arguments():
    parser, args, unparsed_args = default_cli_parser()
    parser.add_argument(
        "c",
        help="delete all atoms below this c fraction",
        type=float,
    )
    parser.parse_args(args=unparsed_args, namespace=args)
    return args


def main():
    args = _parse_command_line_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("ViPErLEED utility: delete atoms above\n")

    if args.c <= 0 or args.c >= 1:
        raise RuntimeError("c must be in range [0, 1]")

    # read the POSCAR file
    slab = readPOSCAR(sys.stdin)


    # process the slab
    modified_slab = deepcopy(slab)
    modified_slab.atlist = [atom for atom in slab.atlist
                            if args.c < atom.pos[2]]

    logger.debug(f"Deleted {len(slab.atlist) - len(modified_slab.atlist)} atoms"
                 f" below c < {args.c}.")

    # write the output file
    writePOSCAR(slab=modified_slab,
                filename=sys.stdout,
                comments='none',
                silent=logger.level<=logging.DEBUG)

if __name__ == "__main__":
    main()
