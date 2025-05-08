"""Module mode of viperleed.calc.bookkeeper.

Defines BookkeeperMode, an enumeration of the known bookkeeper modes.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

from enum import Enum


class BookkeeperMode(Enum):
    """Enumeration of bookkeeper modes.

    Modes ARCHIVE, CLEAR, DISCARD, and DISCARD_FULL emit logging
    warnings if any *_edited file is found (after potentially
    archiving).

    Attributes
    ----------
    ARCHIVE
        Store last run in history. Overwrite PARAMETERS, POSCAR, and
        VIBROCC from OUT. The original inputs in root are renamed to
        *_ori, unless they were edited after calc started and before
        bookkeeper ran. For such files: the edited version is marked
        as *_edited, no *_ori is present, and the corresponding
        "unlabeled" file is pulled from SUPP/original_inputs. Runs
        after run_calc by default.
    CLEAR
        Archive files if they haven't been archived already, then
        clean up the input directory from *_ori files and all the
        already archived outputs of the previous calc execution
        (i.e., log files, OUT, and SUPP). When archiving, *_edited
        files are marked as such, but missing input files are not
        pulled from SUPP/original_inputs. Runs before run_calc by
        default.
    DISCARD
        Re-start from the same input as the previous run. The discarded
        run is kept in history and marked as such. The behavior upon
        archiving is identical to CLEAR. *_ori files are taken as new
        inputs. Has to be run manually after run_calc.
    DISCARD_FULL
        Discard previous run as if it never happened, and remove
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
