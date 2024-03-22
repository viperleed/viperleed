"""ViPErLEED POSCAR utilities."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import argparse
import logging
import pkgutil
import sys

from viperleed.utilities import poscar

POSCAR_UTILITIES = [
    module.name
    for module in pkgutil.iter_modules(poscar.__path__)
    if not module.ispkg and module.name not in {'__main__', 'poscar'}
    ]

poscar_utility_logger = logging.getLogger(__name__)

def add_verbose_option(parser):
    """Add --verbose flag to `parser`."""
    parser.add_argument(
        '-v', '--verbose',
        help='increase output verbosity',
        action='store_true',
        )
