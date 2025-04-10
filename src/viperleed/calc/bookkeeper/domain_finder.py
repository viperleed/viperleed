"""Module domain_finder of viperleed.calc.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-04-10'
__license__ = 'GPLv3+'

import logging

from viperleed.calc.bookkeeper.constants import CALC_LOG_PREFIXES
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.utils import make_property
from viperleed.calc.constants import DEFAULT_DOMAIN_FOLDER_PREFIX
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_WORK_HISTORY
from viperleed.calc.lib.string_utils import harvard_commas


LOGGER = logging.getLogger(__name__)


class DomainFinder:
    """A class for finding domain subfolders of a root directory."""

    def __init__(self, bookkeeper):
        """Initialize an instance from a `bookkeeper`."""
        self._bookkeeper = bookkeeper

    path = make_property('_bookkeeper.cwd')

    def find_potential_domains(self):
        """Return paths to subfolders of self.path that may be domains."""
        not_a_domain = {
            DEFAULT_OUT,
            DEFAULT_SUPP,
            DEFAULT_HISTORY,
            DEFAULT_WORK_HISTORY,
            }
        subfolders = (d for d in self.path.iterdir()
                      if d.is_dir()
                      and d.name not in not_a_domain)
        return tuple(d for d in subfolders
                     if self._is_potential_domain_subfolder(d))
    def _is_potential_domain_subfolder(self, path):
        """Return whether `path` points to a (likely) domain subfolder."""
        if path.name.startswith(DEFAULT_DOMAIN_FOLDER_PREFIX):
            return True
        contents = (
            # A folder is likely a domain subfolder if it contains:
            DEFAULT_OUT,            # an OUT subfolder
            DEFAULT_SUPP,           # or a SUPP subfolder
            DEFAULT_WORK_HISTORY,   # or a workhistory subfolder
            # or it has a calc log file (this means calc
            # was run explicitly in the subfolder)
            *(f'{prefix}*.log' for prefix in CALC_LOG_PREFIXES),
            # or it has been processed before by bookkeeper
            DEFAULT_HISTORY,
            HISTORY_INFO_NAME,
            )
        return any(f for pattern in contents for f in path.glob(pattern))
