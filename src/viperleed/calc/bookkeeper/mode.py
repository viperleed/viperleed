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
    FIX
        Edit information in history and history.info according
        to the most recent format implemented in Bookkeeper.
    """

    ARCHIVE = 'archive'
    CLEAR = 'clear'
    DISCARD = 'discard'
    DISCARD_FULL = 'discard_full'
    FIX = 'fix'

    @property
    def flags(self):
        """Return the CLI flags to select this mode."""
        # First modes that only have a long version
        if self is BookkeeperMode.FIX:
            return (self.long_flag,)
        short_flag = self.value[0]
        if self is BookkeeperMode.DISCARD_FULL:
            short_flag = 'df'
        short_flag = f'-{short_flag}'
        return short_flag, self.long_flag

    @property
    def long_flag(self):
        """Return the long CLI flag for this mode."""
        long_flag = self.value.replace('_', '-')
        return f'--{long_flag}'

    @property
    def uses_ori_files_as_fallback(self):
        """Return whether '*_ori' files are used as fallback for archiving.

        This is the case if an input file in original_inputs is not
        found, but it is present in the current directory with an
        _ori suffix.

        Returns
        -------
        bool
        """
        return self in (BookkeeperMode.CLEAR, BookkeeperMode.DISCARD)
