"""Module guilib of viperleed.

======================================
  ViPErLEED Graphical User Interface
======================================
    *** module guilib.leedsim ***

Schematic simulation of a LEED pattern. Allow exporting list of LEED spots for
properly indexing an experimental LEED pattern
"""

from viperleed import GLOBALS

__all__ = ['LEEDSymmetryDomains', 'LEEDPattern', 'RealSpace',
           'export_pattern_csv', 'LEEDParameters', 'LEEDParametersList',
           'LEEDEquivalentBeams', 'LEEDStructuralDomains', 'LEEDParser']

from viperleed.guilib.leedsim.classes import (LEEDEquivalentBeams,    # not sure I need this to be exposed
                                              LEEDSymmetryDomains,    # not sure I need this to be exposed
                                              LEEDStructuralDomains,  # not sure I need this to be exposed
                                              LEEDParser,             # not sure I need this to be exposed
                                              LEEDPattern, RealSpace)
from viperleed.guilib.leedsim.exportcsv import export_pattern_csv
from viperleed.guilib.leedsim.classes.leedparameters import (
    LEEDParameters, LEEDParametersList
    )

if GLOBALS['USE_GUI']:
    __all__.extend(['DomsBlock', 'EnergyBlock', 'RotationBlock',
                    'ToggleButton', 'LEEDCanvas', 'RealCanvas',
                    'HoverAnnot', 'MatricesPopup', 'TEST',
                    'LatticeInput', 'EditableMatrix', 'BulkInput',
                    'SurfaceStructureInput',
                    'LEEDPatternSimulator',
                    'NewFileDialog',
                    'Bulk3DSymDialog', 'ExportCSVDialog', 'ErrorBox'])
    from viperleed.guilib.leedsim.widgets import (DomsBlock, EnergyBlock,
                                                  RotationBlock, ToggleButton,
                                                  LEEDCanvas, RealCanvas,
                                                  HoverAnnot, MatricesPopup,
                                                  TEST,                         # This probably not...
                                                  LatticeInput, EditableMatrix,
                                                  BulkInput,
                                                  SurfaceStructureInput)
    from viperleed.guilib.leedsim.mainwindow import LEEDPatternSimulator
    from viperleed.guilib.leedsim.mainwindow import LEED_GUI
    from viperleed.guilib.leedsim.mainwindow import show_use_betatest_version_popup
    from viperleed.guilib.leedsim.dialogs import (Bulk3DSymDialog,              # Perhaps not needed if moving all dialogs in their folder
                                                  NewFileDialog,                # Perhaps not needed if moving all dialogs in their folder
                                                  ExportCSVDialog,              # Perhaps not needed if moving all dialogs in their folder
                                                  ErrorBox)
    # from .widgets import *
    # from .mainwindow import LEED_GUI
    # from .NewFileDialog import NewFileDialog
    # from .ExportCSVDialog import ExportCSVDialog
    __all__.extend(['DomsBlock', 'EnergyBlock', 'RotationBlock', 'ToggleButton',
                    'show_use_betatest_version_popup',
                    'LEEDCanvas', 'RealCanvas', 'HoverAnnot', 'MatricesPopup',
                    'TEST', 'LEED_GUI', 'Bulk3DSymDialog', 'NewFileDialog',
                    'ExportCSVDialog'])
