"""
======================================
  ViPErLEED Graphical User Interface
======================================
    *** module guilib.leedsim ***

Schematic simulation of a LEED pattern. Allow exporting list of LEED spots for
properly indexing an experimental LEED pattern
"""

from viperleed import GLOBALS
# from vprglobals import GLOBALS

__all__ = ['LEEDSymmetryDomains', 'LEEDsubpattern', 'Woods', 'LEEDPattern',
           'RealSpace', 'export_pattern_csv', 'LEEDParameters',
           'LEEDParametersList', 'LEEDEquivalentBeams', 'LEEDStructuralDomains']

from viperleed.guilib.leedsim.classes import (LEEDEquivalentBeams,    # not sure I need this to be exposed
                                              LEEDSymmetryDomains,    # not sure I need this to be exposed
                                              LEEDsubpattern, Woods,  # not sure I need these to be exposed, there were also mpl_colors and degrees
                                              LEEDPattern, RealSpace,
                                              LEEDStructuralDomains)  # not sure I need this to be exposed
from viperleed.guilib.leedsim.exportcsv import export_pattern_csv
from viperleed.guilib.leedsim.leedparameters import (LEEDParameters,
                                                     LEEDParametersList)
# from .classes import *
# from .exportcsv import export_pattern_csv
# from .leedparameters import *
if GLOBALS['USE_GUI']:
    from viperleed.guilib.leedsim.widgets import (DomsBlock, EnergyBlock,
                                                  RotationBlock, ToggleButton,
                                                  LEEDCanvas, RealCanvas,
                                                  HoverAnnot, MatricesPopup,
                                                  TEST)                            # This probably not...
    from viperleed.guilib.leedsim.mainwindow import LEED_GUI
    from viperleed.guilib.leedsim.dialogbulk3dsym import Bulk3DSymDialog
    from viperleed.guilib.leedsim.newfiledialog import NewFileDialog
    from viperleed.guilib.leedsim.exportcsvdialog import ExportCSVDialog
    # from .widgets import *
    # from .mainwindow import LEED_GUI
    # from .NewFileDialog import NewFileDialog
    # from .ExportCSVDialog import ExportCSVDialog
    __all__.extend(['DomsBlock', 'EnergyBlock', 'RotationBlock', 'ToggleButton',
                    'LEEDCanvas', 'RealCanvas', 'HoverAnnot', 'MatricesPopup',
                    'TEST', 'LEED_GUI', 'Bulk3DSymDialog', 'NewFileDialog',
                    'ExportCSVDialog'])


