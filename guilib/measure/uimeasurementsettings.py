"""Module uimeasurementsettings of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-10-29
Author: Michele Riva
Author: Florian Doerr

Defines the SettingsEditor class.
"""
import configparser
from pathlib import Path

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw
import PyQt5.QtGui as qtg

# ViPErLEED modules
from viperleed import guilib as gl

TITLE = 'Measurement Settings'

class SettingsEditor(gl.ViPErLEEDPluginBase):
    """A class that allows simulating LEED Patterns."""

    def __init__(self, parent=None):
        """Initialize window."""
        super().__init__(parent, name=TITLE)
        # Keep references to controls, dialogs, and some globals
        self._ctrls = {
            }