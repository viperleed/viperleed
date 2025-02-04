"""Module constants of viperleed.calc.bookkeeper.

Collects constants used in various submodules of bookkeeper.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-10'
__license__ = 'GPLv3+'

import re

from viperleed.calc.constants import LOG_PREFIX

CALC_LOG_PREFIXES = (
    LOG_PREFIX,
    # The next ones are for backwards compatibility. Must
    # be in historical order of use: Most recent first!
    'tleedm',
    )
HISTORY_FOLDER_RE = re.compile(
    r't(?P<tensor_num>[0-9]{3})\.r(?P<job_num>[0-9]{3})(?P<rest>.*)'
    )

# Input/output files that may have _ori or _OUT suffix
STATE_FILES = ('PARAMETERS', 'POSCAR', 'VIBROCC')


# SUFFIXES FOR INPUT FILES
EDITED_SUFFIX = '_edited'  # Edited after calc and before bookkeeper
ORI_SUFFIX = '_ori'        # Non-edited input file; used at ARCHIVE
