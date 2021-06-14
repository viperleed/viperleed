"""
=======================================
   ViPErLEED Graphical User Interface
=======================================
 *** module guilib.leedsim.widgets ***

Created: 2020-01-12
Author: Michele Riva

Blah blah TODO
"""
import PyQt5.QtWidgets as qtw


class ToggleButton(qtw.QPushButton):

    def sizeHint(self):
        self.ensurePolished()
        refbutton = qtw.QPushButton(self.text())
        refbutton.setFont(self.font())
        return refbutton.sizeHint()*1.2

    def minimumSizeHint(self):
        return self.sizeHint()
