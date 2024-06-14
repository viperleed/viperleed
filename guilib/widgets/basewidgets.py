"""Module basewidgets of viperleed.guilib.widgets.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-06-13
Author: Michele Riva
Author: Florian Dörr

This module defines basic, non-specific widgets.
"""

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw


class ButtonWithLabel(qtw.QWidget):
    """QLabel and QPushButton."""

    def __init__(self, tight=True, **kwargs):
        """Initialise widget.

        Parameters
        ----------
        tight : bool, optional
            Whether the layout should not have
            spacing around it. Default is true.

        Returns
        -------
        None
        """
        super().__init__(**kwargs)
        self.label = qtw.QLabel()
        self.button = qtw.QPushButton()
        self._compose(tight)

    def _compose(self, tight):
        """Compose widget."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.button)
        if tight:
            layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

    @qtc.pyqtSlot(str)
    def set_label_text(self, text):
        """Set the text of the label."""
        self.label.setText(text)

    @qtc.pyqtSlot(str)
    def set_button_text(self, text):
        """Set the text of the button."""
        self.button.setText(text)
