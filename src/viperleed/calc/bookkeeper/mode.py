"""Module mode of viperleed.calc.bookkeeper.

Defines BookkeeperMode, an enumeration of the known bookkeeper modes.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from enum import Enum


class BookkeeperMode(Enum):
    """Enumeration of bookkeeper modes.

    Attributes
    ----------
    ARCHIVE
        Store last run in history. Overwrite PARAMETERS, POSCAR, and
        VIBROCC from OUT. Runs after run_calc by default.
    CLEAR
        Clear the input directory of the last run, archiving files
        beforehand if they haven't been archived already. During
        archiving, it takes the '_ori'-suffixed files, if present,
        to create the history folder. Runs before run_calc by default.
    DISCARD
        Re-start from the same input as the previous run. The
        discarded run is kept in history and marked as such (the
        behavior upon archiving is identical to CLEAR). Has to
        be run manually after run_calc.
    DISCARD_FULL
        Discard previous run as if it never happened and remove
        it from history. Has to be run manually after run_calc.
    """

    ARCHIVE = 'archive'
    CLEAR = 'clear'
    DISCARD = 'discard'
    DISCARD_FULL = 'discard_full'

    @property
    def uses_ori_files_as_fallback(self):
        """Return whether '*_ori' files are used as fallback for archiving."""
        return self is BookkeeperMode.CLEAR
