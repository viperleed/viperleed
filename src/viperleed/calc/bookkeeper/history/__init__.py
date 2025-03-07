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
explorer
    Defines functionality to collect information from the 'history'
    folder, removing subfolders from it. It also provides an interface
    for accessing the history.info file.
folder
    Defines classes used in explorer to gather information from
    subfolders of 'history'.
info
    Defines the HistoryInfoFile class, the handler for the
    history.info file that takes care of reading/writing
    and consistent formatting.
meta
    Defines functionality for creating and reading the metadata file
    stored in history subfolders by Bookkeeper. This file is used to
    keep track of general information about the viperleed.calc run
    that was archived. It provides a different mechanism from the
    history.info file to store lower-level information that should
    not be of interest to the user.
workhistory
    Defines functionality to collect information from and archive the
    contents of 'workhistory' folders, created by viperleed.calc when
    execution of an "earlier" segment follows execution of a later one,
    as in looped delta--search pairs, or in a run of the form 1-3 1.

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
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'
