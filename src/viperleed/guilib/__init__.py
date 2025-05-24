"""Module viperleed.guilib.

======================================
  ViPErLEED Graphical User Interface
======================================
       *** module guilib ***

Created: 2020-01-10
Author: Michele Riva
"""

import os

from viperleed import GLOBALS
from viperleed.guilib.base import check_py_version

USE_GUI = GLOBALS['USE_GUI']

if (os.name == 'posix' and 'DISPLAY' not in os.environ.keys()) or not USE_GUI:
    # The environment does not have graphics capabilities.
    BACKEND = None
    GLOBALS['USE_GUI'] = False
else:
    BACKEND = 'mplcairo' if check_py_version('3.8', 'earlier') else 'agg'
