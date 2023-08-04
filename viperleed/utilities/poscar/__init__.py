"""ViPErLEED POSCAR utilities.

@author: Alexander M. Imre
created: 2023-08-03
"""

import argparse
import logging
import sys

poscar_utility_logger = logging.getLogger("viperleed.utilities.poscar")

def default_cli_parser():
    # if executed from the terminal, send all logs to stderr because stdout is
    # used for piping out the POSCAR file
    poscar_utility_logger.addHandler(logging.StreamHandler(sys.stderr))

    # always have a verbose option
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-v", "--verbose",
        help=("increase output verbosity"),
        action="store_true",
    )

    args, unknown_args = parser.parse_known_args()

    # re add the help option
    parser.add_argument(
        '-h', '--help',
        action='help', default=argparse.SUPPRESS,
        help='show this help message and exit'
    )

    # verbose option shows debug messages and stacktrace on error
    if args.verbose:
        poscar_utility_logger.setLevel(logging.DEBUG)
    else:
        # don't show entire backtrace on error
        sys.tracebacklimit = 0

    return parser, args, unknown_args
