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
        Store last run in history. Overwrite PARAMETERS, POSCAR & VIBROCC from
        OUT. Runs after run_calc by default.
    CLEAR
        Clear the input directory of last run.
        Runs before run_calc by default.
    DISCARD
        Re-start from the same input as the previous run. The discarded run is
        kept in history. Has to be run manually after run_calc.
    DISCARD_FULL
        Discard previous run as if it never happened and removes it from
        history. Has to be run manually after run_calc.
    """
    ARCHIVE = 'archive'
    CLEAR = 'clear'
    DISCARD = 'discard'
    DISCARD_FULL = 'discard_full'

    @property
    def discard(self):
        """Return whether this is mode DISCARD."""
        return self is BookkeeperMode.DISCARD
