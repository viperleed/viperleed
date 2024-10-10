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

HISTORY_FOLDER_RE = re.compile(
    r't(?P<tensor_num>[0-9]{3})\.r(?P<job_num>[0-9]{3})(?P<rest>.*)'
    )
