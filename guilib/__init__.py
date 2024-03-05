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

USE_GUI = GLOBALS['USE_GUI']

from viperleed.guilib.base import (get_equivalent_beams,
                                   project_to_first_domain, check_type,
                                   check_leed_params, check_multi_leed_params,
                                   catch_gui_crash, check_py_version,
                                   string_matrix_to_numpy,
                                   format_floats,  # probably not needed globally
                                   integer_part_length,
                                   parallel, orientation, screen_radius)

from viperleed.guilib.mathparse import MathParser, UnsupportedMathError

if (os.name == 'posix' and 'DISPLAY' not in os.environ.keys()) or not USE_GUI:
    # The environment does not have graphics capabilities.
    BACKEND = None
    GLOBALS['USE_GUI'] = False
else:
    # Import GUI modules
    if check_py_version('3.8', 'earlier'):
        BACKEND = 'mplcairo'
    else:
        BACKEND = 'agg'
    from viperleed.guilib.basewidgets import (Figure, FigureCanvas,
                                              MPLFigureCanvas, PainterMatrix,
                                              TextBox, TextBoxWithButtons)
    from viperleed.guilib.widgetslib import (AllGUIFonts, drawText,
                                             change_control_text_color,
                                             editStyleSheet,
                                             get_all_children_widgets)
    from viperleed.guilib.pluginsbase import (ViPErLEEDPluginBase,
                                              AboutViPErLEED)
    from viperleed.guilib.leedsim import *
    from viperleed.guilib.measure import *
    from viperleed.guilib.selectplugin import ViPErLEEDSelectPlugin


from viperleed.guilib.leedsim.exportcsv import export_pattern_csv
from viperleed.guilib.leedsim.classes import (LEEDParameters,
                                              LEEDParametersList)

# Perhaps better to have the * import from up there moved down here, as
# leedsim already takes care of selecting only Qt-independent classes
# when there is no X-server by looking at GLOBALS['USE_GUI']
# I would then move also the * import of measure, making sure it also
# cares about the __all__ correctly
from viperleed.guilib.leedsim.classes import (LEEDPattern, LEEDSymmetryDomains,
                                              LEEDEquivalentBeams,
                                              LEEDStructuralDomains,
                                              LEEDParser)
