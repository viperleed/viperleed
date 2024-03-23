"""Module cli of viperleed.

Defines the main command-line interface for running viperleed
and its subpackages.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-01'
__license__ = 'GPLv3+'

from argparse import ArgumentParser

from viperleed import GLOBALS
from viperleed.calc.bookkeeper import bookkeeper_cli_options
from viperleed.calc.bookkeeper import main as bookkeeper_main
from viperleed.calc.cli import add_calc_parser_arguments
from viperleed.calc.cli import main as main_calc
from viperleed.gui import main as gui_main                                      # TODO: gui arguments
from viperleed.utilities.__main__ import add_util_parser_arguments
from viperleed.utilities.__main__ import main as utilities_main
from viperleed.utilities.poscar.__main__ import add_poscar_parser_arguments
from viperleed.utilities.poscar.__main__ import main as poscar_main


def main():
    """Execute ViPErLEED from the command line."""
    viperleed_parser = ArgumentParser(prog='viperleed')
    viperleed_parser.add_argument(
        '--version',
        help='print version number',
        action='version',
        version=GLOBALS['version_message'],
        )
    # get subparsers
    subparsers = viperleed_parser.add_subparsers()

    # viperleed bookkeeper
    parser_bookkeeper = subparsers.add_parser('bookkeeper')
    bookkeeper_cli_options(parser_bookkeeper)
    parser_bookkeeper.set_defaults(func=bookkeeper_main)

    # viperleed calc
    parser_calc = subparsers.add_parser('calc',)
    add_calc_parser_arguments(parser_calc)
    parser_calc.set_defaults(func=main_calc)

    # viperleed gui
    parser_gui = subparsers.add_parser('gui')
    parser_gui.set_defaults(func=gui_main)

    # viperleed utilities
    parser_util = subparsers.add_parser('util')
    add_util_parser_arguments(parser_util)
    parser_util.set_defaults(func=utilities_main)

    # viperleed poscar utilities
    parser_poscar = subparsers.add_parser('poscar')
    add_poscar_parser_arguments(parser_poscar)
    parser_poscar.set_defaults(func=poscar_main)

    args = viperleed_parser.parse_args()

    try:
        return args.func(args)
    except AttributeError:
        viperleed_parser.print_help()
