"""Module constants of viperleed.calc.bookkeeper.

Collects constants used in various submodules of bookkeeper.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
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
# A suffix for input files that were edited by the user since
# the most-recent calc execution, and before bookkeeper ran.
EDITED_SUFFIX = '_edited'
HISTORY_FOLDER_RE = re.compile(
    r't(?P<tensor_num>\d{3,})\.r(?P<job_num>\d{3,})(?P<rest>.*)'
    )

# Input/output files that may have _ori or, before #302, _OUT suffix
STATE_FILES = ('PARAMETERS', 'POSCAR', 'VIBROCC')


# SUFFIXES FOR INPUT FILES
EDITED_SUFFIX = '_edited'  # Edited after calc and before bookkeeper
ORI_SUFFIX = '_ori'        # Non-edited input file; used at ARCHIVE
