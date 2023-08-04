"""ViPErLEED POSCAR utilities.

@author: Alexander M. Imre
created: 2023-08-03
"""
import argparse
import logging
import pkgutil
import sys

from viperleed.utilities import poscar
from viperleed.utilities.poscar import *

POSCAR_UTILITIES = [module.name for module in
                    pkgutil.iter_modules(poscar.__path__)
                    if (not module.ispkg and
                        module.name != '__main__' and
                        module.name != 'poscar')]
print(POSCAR_UTILITIES)

poscar_utility_logger = logging.getLogger("viperleed.utilities.poscar")

def add_verbose_option(parser):
        parser.add_argument(
        "-v", "--verbose",
        help=("increase output verbosity"),
        action="store_true",
    )
