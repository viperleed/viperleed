"""Module modulebase of viperleed.guilib.

======================================
  ViPErLEED Graphical User Interface
======================================

Defines the ViPErLEEDModuleBase class from which all ViPErLEED
modules should inherit. Any concrete implementation of a module
should use super().closeEvent() if it wants to accept a closeEvent()
rather than just .accept()ing the event.

Author: Michele Riva
Created: 2021-06-29
"""

from PyQt5 import (QtWidgets as qtw,
                   QtCore as qtc,
                   QtGui as qtg)

from viperleed.gui import resources_path
from viperleed import guilib as gl

LOGO = resources_path('guilib/icons/viperleed_logo_circled_48x48.png')


class ViPErLEEDModuleBase(qtw.QMainWindow):
    """Base class for a ViPErLEED module."""

    module_closed = qtc.pyqtSignal(object)  # The class being destroyed

    def __init__(self, parent=None):
        """Initialize module."""
        super().__init__(parent)
        self.setAttribute(qtc.Qt.WA_DeleteOnClose)
        self.setWindowIcon(qtg.QIcon(LOGO))

    def closeEvent(self, event):
        """Reimplement closeEvent to emit a module_closed."""
        self.module_closed.emit(self)
        super().closeEvent(event)