"""ViPErLEED POSCAR utilities."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__created__ = '2023-08-03'

import argparse
import logging
import pkgutil
import sys

from viperleed.utilities import poscar

POSCAR_UTILITIES = [module.name for module in
                    pkgutil.iter_modules(poscar.__path__)
                    if (not module.ispkg and
                        module.name != '__main__' and
                        module.name != 'poscar')]

poscar_utility_logger = logging.getLogger("viperleed.utilities.poscar")

def add_verbose_option(parser):
        parser.add_argument(
        "-v", "--verbose",
        help=("increase output verbosity"),
        action="store_true",
    )
