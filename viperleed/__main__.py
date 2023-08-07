#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
=================
    ViPErLEED
=================
"""
from argparse import ArgumentParser

from viperleed import GLOBALS
from viperleed.bookkeeper import bookkeeper_cli_options
from viperleed.bookkeeper import main as bookkeeper_main
from viperleed.calc.__main__ import add_calc_parser_arguments
from viperleed.calc.__main__ import main as main_calc
from viperleed.gui import main as gui_main                                      # TODO: gui arguments
from viperleed.utilities.__main__ import add_util_parser_arguments
from viperleed.utilities.__main__ import main as utilities_main
from viperleed.utilities.poscar.__main__ import add_poscar_parser_arguments
from viperleed.utilities.poscar.__main__ import main as poscar_main


def main():
    viperleed_parser = ArgumentParser(prog="viperleed")
    viperleed_parser.add_argument(
        '--version',
        help='print version number',
        action='version',
        version=GLOBALS['version_message'],
    )
    # get subparsers
    subparsers = viperleed_parser.add_subparsers()

    # viperleed bookkeeper
    parser_bookkeeper = subparsers.add_parser("bookkeeper")
    bookkeeper_cli_options(parser_bookkeeper)
    parser_bookkeeper.set_defaults(func=bookkeeper_main)

    # viperleed calc
    parser_calc = subparsers.add_parser("calc",)
    add_calc_parser_arguments(parser_calc)
    parser_calc.set_defaults(func=main_calc)

    # viperleed gui
    parser_gui = subparsers.add_parser("gui")
    parser_gui.set_defaults(func=gui_main)

    # viperleed utilities
    parser_util = subparsers.add_parser("util")
    add_util_parser_arguments(parser_util)
    parser_util.set_defaults(func=utilities_main)

    # viperleed poscar utilities
    parser_poscar = subparsers.add_parser("poscar")
    add_poscar_parser_arguments(parser_poscar)
    parser_poscar.set_defaults(func=poscar_main)

    args = viperleed_parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
