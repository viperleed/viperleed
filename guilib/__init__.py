"""
======================================
  ViPErLEED Graphical User Interface
======================================
       *** module guilib ***

Created: 2020-01-10
Author: Michele Riva
       
"""

import os
import sys

cd = os.path.realpath(os.path.dirname(__file__))
guilib_path = os.path.realpath(os.path.join(cd, '..'))
for import_path in (cd, guilib_path):
    if import_path not in sys.path:
        sys.path.append(import_path)
        
from vprglobals import GLOBALS

USE_GUI = GLOBALS['USE_GUI']

if (os.name == 'posix' and 'DISPLAY' not in os.environ.keys()) or not USE_GUI:
    # The environment does not have graphics capabilities.
    BACKEND = None
    # print("Command line only")
else:
    # Import GUI modules
    # print("Environment has GUI")
    BACKEND = 'mplcairo'
    from .decorators import *
    from .widgetdecorators import *
    from .basewidgets import *
    from .widgetslib import *
    from .leedsim import *
    from .measure import *

from .base import *
from .leedsim.exportcsv import export_pattern_csv
