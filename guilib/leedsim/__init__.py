"""Module guilib of viperleed.

======================================
  ViPErLEED Graphical User Interface
======================================
    *** module guilib.leedsim ***

Schematic simulation of a LEED pattern. Allow exporting list of LEED spots for
properly indexing an experimental LEED pattern
"""

from viperleed import GLOBALS

__all__ = ['LEEDSymmetryDomains', 'Woods', 'LEEDPattern', 'RealSpace',
           'export_pattern_csv', 'LEEDParameters', 'LEEDParametersList',
           'LEEDEquivalentBeams', 'LEEDStructuralDomains', 'LEEDParser']

from viperleed.guilib.leedsim.classes import (LEEDEquivalentBeams,    # not sure I need this to be exposed
                                              LEEDSymmetryDomains,    # not sure I need this to be exposed
                                              LEEDStructuralDomains,  # not sure I need this to be exposed
                                              LEEDParser,             # not sure I need this to be exposed
                                              Woods,                  # not sure I need this to be exposed
                                              LEEDPattern, RealSpace)
from viperleed.guilib.leedsim.exportcsv import export_pattern_csv
from viperleed.guilib.leedsim.leedparameters import (LEEDParameters,
                                                     LEEDParametersList)

if GLOBALS['USE_GUI']:
    __all__.extend(['DomsBlock', 'EnergyBlock', 'RotationBlock',
                    'ToggleButton', 'LEEDCanvas', 'RealCanvas',
                    'HoverAnnot', 'MatricesPopup', 'TEST',
                    'LatticeInput', 'EditableMatrix', 'BulkInput',
                    'SurfaceStructureInput',
                    'LEED_GUI',
                    'NewFileDialog',
                    'Bulk3DSymDialog', 'ExportCSVDialog'])
    from viperleed.guilib.leedsim.widgets import (DomsBlock, EnergyBlock,
                                                  RotationBlock, ToggleButton,
                                                  LEEDCanvas, RealCanvas,
                                                  HoverAnnot, MatricesPopup,
                                                  TEST,                         # This probably not...
                                                  LatticeInput, EditableMatrix,
                                                  BulkInput,
                                                  SurfaceStructureInput)
    from viperleed.guilib.leedsim.mainwindow import LEED_GUI
    from viperleed.guilib.leedsim.dialogs import (Bulk3DSymDialog,              # Perhaps not needed if moving all dialogs in their folder
                                                  NewFileDialog,                # Perhaps not needed if moving all dialogs in their folder
                                                  ExportCSVDialog)              # Perhaps not needed if moving all dialogs in their folder
