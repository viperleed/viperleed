#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ViPErLEED utilities.
"""
from viperleed.utilities.poscar.__main__ import main as poscar_main
from viperleed.utilities.poscar.__main__ import add_poscar_parser_arguments

def add_util_parser_arguments(parser):
    subparsers = parser.add_subparsers()

    poscar_util_parser = subparsers.add_parser(                                 # TODO: add poscar utilities
        "poscar",
        help="utilities for POSCAR files"
    )
    add_poscar_parser_arguments(poscar_util_parser)
    poscar_util_parser.set_defaults(func=poscar_main)

    aux_to_exp_parser = subparsers.add_parser(
        "AUXEXPBEAMS_to_EXPBEAMS",
        help="call utility to convert AUXEXPBEAMS to EXPBEAMS"
    )
    aux_to_exp_parser.set_defaults(func=None)

    rearrange_ps_parser = subparsers.add_parser(
        "rearrange_phaseshifts",
        help="call utility to rearrange phaseshifts"
    )
    rearrange_ps_parser.set_defaults(func=None)

def main(args=None):
    print("ViPErLEED utilities.")

if __name__ == "__main__":
    main()