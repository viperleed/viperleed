"""Package bookkeeper of viperleed.calc.

Defines functionality for keeping track of calculation history,
as in a logbook. This package is a result of the refactor of
the original bookkeeper.py module.

Modules
-------
__main__
    Entry point for the bookkeeper when run via python -m bookkeeper.
bookkeeper
    Defines the Bookkeeper class that handles the whole process.
cli
    Defines the command-line interface for invoking the bookkeeper
    as a standalone utility.
log
    Defines the bookkeeper logger instance and its related
    functionality.
mode
    Defines the BookkeeperMode enumeration of the possible modes
    of execution of the bookkeeper. See help(BookkeeperMode) for
    information on what each mode entails.

Packages
--------
history
    Defines the functionality for handling the history.info file
    and its contents, including reading/writing and editing formats.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'
