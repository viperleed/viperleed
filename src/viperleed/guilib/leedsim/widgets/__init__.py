"""Module widgets of viperleed.guilib.leedsim.

======================================
  ViPErLEED Graphical User Interface
======================================

Contains QWidgets used by within the pattern simulator
plug-in of the viperleed Graphical User Interface.

Created: 2021-06-01
Author: Michele Riva
"""

__all__ = ['LatticeInput', 'EditableMatrix', 'BulkInput',
           'SurfaceStructureInput',
           'ToggleButton', 'DomsBlock', 'EnergyBlock', 'RotationBlock',
           'HoverAnnot', 'LEEDCanvas', 'RealCanvas', 'MatricesPopup', 'TEST']

from viperleed.guilib.leedsim.widgets.latticeinput import LatticeInput
from viperleed.guilib.leedsim.widgets.editablematrix import EditableMatrix
from viperleed.guilib.leedsim.widgets.bulkinput import BulkInput
from viperleed.guilib.leedsim.widgets.surfaceinput import SurfaceStructureInput

from viperleed.guilib.leedsim.widgets.togglebutton import ToggleButton
from viperleed.guilib.leedsim.widgets.domainsblock import DomsBlock
from viperleed.guilib.leedsim.widgets.energyblock import EnergyBlock
from viperleed.guilib.leedsim.widgets.rotationblock import RotationBlock
from viperleed.guilib.leedsim.widgets.hoverannot import HoverAnnot
from viperleed.guilib.leedsim.widgets.leedcanvas import (LEEDCanvas, TEST)
from viperleed.guilib.leedsim.widgets.realcanvas import RealCanvas
from viperleed.guilib.leedsim.widgets.matricespopup import MatricesPopup