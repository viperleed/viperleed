"""Module guilib of viperleed.

======================================
  ViPErLEED Graphical User Interface
======================================
    *** module guilib.leedsim ***

Schematic simulation of a LEED pattern. Allow exporting list of LEED spots for
properly indexing an experimental LEED pattern
"""

from viperleed.vprglobals import GLOBALS

__all__ = ['LEEDSymmetryDomains',
           # 'LEEDsubpattern',  # OLD
           'Woods', 'LEEDPattern',
           'RealSpace', 'export_pattern_csv', 'LEEDParameters',
           'LEEDParametersList', 'LEEDEquivalentBeams', 'LEEDStructuralDomains']

from viperleed.guilib.leedsim.classes import (LEEDEquivalentBeams,    # not sure I need this to be exposed
                                              LEEDSymmetryDomains,    # not sure I need this to be exposed
                                              LEEDStructuralDomains,  # not sure I need this to be exposed
                                              # LEEDsubpattern,         # OLD
                                              Woods,                  # not sure I need these to be exposed
                                              LEEDPattern, RealSpace)
from viperleed.guilib.leedsim.exportcsv import export_pattern_csv
from viperleed.guilib.leedsim.leedparameters import (LEEDParameters,
                                                     LEEDParametersList)
if GLOBALS['USE_GUI']:
    from viperleed.guilib.leedsim.widgets import (DomsBlock, EnergyBlock,
                                                  RotationBlock, ToggleButton,
                                                  LEEDCanvas, RealCanvas,
                                                  HoverAnnot, MatricesPopup,
                                                  TEST)                            # This probably not...
    from viperleed.guilib.leedsim.mainwindow import LEED_GUI
    from viperleed.guilib.leedsim.newfiledialog import NewFileDialog
    from viperleed.guilib.leedsim.exportcsvdialog import ExportCSVDialog
    __all__.extend(['DomsBlock', 'EnergyBlock', 'RotationBlock', 'ToggleButton',
                    'LEEDCanvas', 'RealCanvas', 'HoverAnnot', 'MatricesPopup',
                    'TEST', 'LEED_GUI', 'NewFileDialog', 'ExportCSVDialog'])


