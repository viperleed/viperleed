"""
======================================
  ViPErLEED Graphical User Interface
======================================
       *** module guilib ***

Created: 2020-01-10
Author: Michele Riva
       
"""

import os

from viperleed.vprglobals import GLOBALS

USE_GUI = GLOBALS['USE_GUI']

from viperleed.guilib.decorators import profile_calls, exec_time, profile_lines
from viperleed.guilib.base import (get_equivalent_beams,
                                   project_to_first_domain, check_type,
                                   check_leed_params, check_multi_leed_params,
                                   catch_gui_crash, string_matrix_to_numpy,
                                   format_floats,  # probably not needed globally
                                   integer_part_length,
                                   parallel, orientation, screen_radius,
                                   BeamIndex, PlaneGroup, Lattice)

if (os.name == 'posix' and 'DISPLAY' not in os.environ.keys()) or not USE_GUI:
    # The environment does not have graphics capabilities.
    BACKEND = None
else:
    # Import GUI modules
    BACKEND = 'mplcairo'
    from viperleed.guilib.decorators import ensure_decorates_class
    from viperleed.guilib.widgetdecorators import (receive_mouse_broadcast,
                                                   broadcast_mouse)
    from viperleed.guilib.basewidgets import (Figure, FigureCanvas,
                                              MPLFigureCanvas, PainterMatrix,
                                              TextBox, TextBoxWithButtons)
    from viperleed.guilib.widgetslib import (AllGUIFonts, drawText,
                                             editStyleSheet,
                                             get_all_children_widgets)
    from viperleed.guilib.leedsim import *
    from viperleed.guilib.measure import *


from viperleed.guilib.leedsim.exportcsv import export_pattern_csv
from viperleed.guilib.leedsim.leedparameters import LEEDParameters
from viperleed.guilib.leedsim.leedparameters import LEEDParametersList
from viperleed.guilib.leedsim.classes import LEEDPattern, LEEDSymmetryDomains