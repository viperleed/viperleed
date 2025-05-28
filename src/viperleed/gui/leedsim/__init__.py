"""Package leedsim of viperleed.gui.

Schematic simulation of a LEED pattern. Allows exporting a list of LEED
spots (i.e., the pattern file) for properly indexing an experimental
LEED pattern.

Modules
-------
exportcsv
    Qt-independent functionality for creating a "Spot Pattern File".
mainwindow
    The front panel of the LEED Pattern Simulator. Requires PyQt5.
mainwindow_old
    An older version of the front panel that supports only a single
    structural domain an normal beam incidence. Requires PyQt5.
utils
    Qt-independent functions used in multiple spots in leedsim. Exists
    mostly to prevent circular imports.

Packages
--------
classes
    Qt-independent objects used in leedsim.
dialogs
    Dialog-like (i.e., free-standing windows) Qt widgets. Requires
    PyQt5.
widgets
    Qt QWidgets for leedsim.

Non-python data
---------------
exported
    Examples of Spot Pattern files exported via the LEED Pattern
    Simulator.
input_examples
    Sample input files for the LEED Pattern Simulator.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-03'
__license__ = 'GPLv3+'
