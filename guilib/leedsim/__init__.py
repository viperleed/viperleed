"""
======================================
  ViPErLEED Graphical User Interface
======================================
    *** module guilib.leedsim ***

Schematic simulation of a LEED pattern. Allow exporting list of LEED spots for
properly indexing an experimental LEED pattern
"""

from vprglobals import GLOBALS

from .classes import *
from .exportcsv import export_pattern_csv
if GLOBALS['USE_GUI']:
    from .widgets import *
    from .mainwindow import LEED_GUI
    from .NewFileDialog import NewFileDialog
    from .ExportCSVDialog import ExportCSVDialog

def print_inspect():
    import sys, inspect
    print('Here are the functions and classes available:')
    for name, obj in inspect.getmembers(sys.modules[__name__]):
        #if inspect.isclass(obj):
        print(obj)

# print('You have imported', __name__)
# print_inspect()

