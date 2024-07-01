"""Module constants of viperleed.calc.bookkeeper.

Collects the definitions of a few constants used
in multiple modules of the bookkeeper package.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

import logging
import re

from viperleed.calc import LOG_PREFIX
from viperleed.calc.lib.base import parent_name


LOGGER = logging.getLogger(parent_name(__name__))

CALC_LOG_PREFIXES = (
    LOG_PREFIX,
    'tleedm',   # For backwards compatibility
    )
HIST_FOLDER_RE = re.compile(
    r't(?P<tensor_num>[0-9]{3}).r(?P<job_num>[0-9]{3})_'
    )
HISTORY_INFO_NAME = 'history.info'

STATE_FILES = ('PARAMETERS', 'POSCAR', 'VIBROCC')
# optional input files that may be generated at runtime - do not warn if missing
RUNTIME_GENERATED_INPUT_FILES = ('IVBEAMS', 'PHASESHIFTS')
