"""Package gui of viperleed.

======================================
  ViPErLEED Graphical User Interface
======================================

Modules
-------
base
    Qt-independent functions used in multiple spots in the GUI.
    Defines also the API use by calc for listing equivalent beams.
basewidgets
    Basic Qt widgets used in other parts of the GUI. Requires PyQt5.
cli
    The startup functionality for invoking the GUI in graphical
    or command-line mode.
constants
    Constant Qt-independent values used in multiple places in the GUI.
decorators
    Qt-independent decorators. Mostly used for debugging and profiling
    during development.
detect_graphics
    Functionality for determining whether the system has graphical
    capabilities.
helpers
    Other Qt-independent functions not defined in base. Exists mostly
    to avoid circular imports.
mathparse
    A simple parser for mathematical expressions. Qt-independent.
mpl_graphics
    Functionality for interaction with the matplotlib module.
    Particularly for what concerns the selection of the correct
    plotting engine.
pluginsbase
    Base functionality for any of the functional modules available
    in the graphical user interface. Requires PyQt5.
selecplugin
    The main window of the GUI, allowing interactive selection of
    which part of the GUI to use. Requires PyQt5.

Packages
--------
classes
    Qt-independent objects used in multiple places in the GUI
    as well as in calc.
leedsim
    The LEED Pattern Simulator plugin.
measure
    The LEED-IV measurement plugin.
widgets
    Custom widgets and widget-related functionality. Requires PyQt5.

Non-python data
---------------
fonts
    Fonts used in the GUI.
icons
    Custom icons for the GUI.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-10'
__license__ = 'GPLv3+'
