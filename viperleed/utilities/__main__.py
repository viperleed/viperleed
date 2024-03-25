"""ViPErLEED utilities."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import argparse

from viperleed.utilities import SATLEED_to_EXPBEAMS
from viperleed.utilities import rearrange_phaseshifts
from viperleed.utilities.poscar.__main__ import add_poscar_parser_arguments
from viperleed.utilities.poscar.__main__ import main as poscar_main


def add_util_parser_arguments(parser):
    subparsers = parser.add_subparsers()

    poscar_util_parser = subparsers.add_parser(
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

    satleed_parser = subparsers.add_parser(
        "SATLEED_to_EXPBEAMS",
        help="call utility to convert AUXEXPBEAMS to EXPBEAMS"
    )
    satleed_parser.set_defaults(func=SATLEED_to_EXPBEAMS.main)
    SATLEED_to_EXPBEAMS.add_cli_parser_arguments(satleed_parser)

    rearrange_ps_parser = subparsers.add_parser(
        "rearrange_phaseshifts",
        help="call utility to rearrange phaseshifts"
    )
    rearrange_ps_parser.set_defaults(func=rearrange_phaseshifts.main)

def main(args=None):
    if args is None:
        parser = argparse.ArgumentParser(prog='viperleed.utilities')
        add_util_parser_arguments(parser)
        args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:
        parser.print_help()

if __name__ == "__main__":
    main()