"""Module spinboxes of viperleed.guilib.measure.widgets.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-12-06
Author: Michele Riva

Defines convenience subclasses of QSpinBox and QDoubleSpinBox.
"""

from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw,
                   QtGui as qtg)


class TolerantCommaSpinBox(qtw.QDoubleSpinBox):
    """A QDoubleSpinBox that tolerates both '.' and ',' as separators.

    Each time the user presses the "comma" key, it is interpreted
    as a if a "dot" was pressed instead. This is useful to keep
    intact the user experience of pressing the "comma" key as a
    decimal separator while enforcing "dot" as a decimal separator
    everywhere.

    Reimplemented methods
    ---------------------
    keyPressEvent(event)
        Replace comma key-presses with dot key presses.
    """

    def keyPressEvent(self, event):  # pylint: disable=invalid-name
        """Replace commas with dots."""
        key = event.key()
        if key == qtc.Qt.Key_Comma:
            # Simulate that a dot was actually pressed
            event = qtg.QKeyEvent(
                event.type(), qtc.Qt.Key_Period, event.modifiers(),
                text=event.text().replace(',', '.'),
                autorep=event.isAutoRepeat(), count=event.count()
                )
        super().keyPressEvent(event)


class InfIntSpinBox(qtw.QDoubleSpinBox):
    """A spin-box that allows inputting an infinite integer number."""

    def __init__(self, parent=None):
        """Initialize instance."""
        super().__init__(parent)
        self.setRange(0, float('inf'))
        self.setDecimals(0)
