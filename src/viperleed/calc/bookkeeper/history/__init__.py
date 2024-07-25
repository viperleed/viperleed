"""Package history of viperleed.calc.bookkeeper.

Defines the functionality for handling the history.info file
and its contents, including reading/writing and editing formats.

Modules
-------
constants
    Collects constants that are use by the submodules/packages.
    It is intended exclusively to prevent cyclic-import issues.
errors
    Defines the HistoryInfoError and its subclasses, signaling
    unexpected situations occurring when dealing with the file.
file
    Defines the HistoryInfoFile class, the handler for the
    history.info file that takes care of reading/writing
    and consistent formatting.

Packages
--------
entry
    Defines the functionality that pertains to each separate
    'block' of a history.info file.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'
