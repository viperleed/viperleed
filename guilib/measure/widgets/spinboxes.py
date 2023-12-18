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
    keyReleaseEvent(event)
        Replace comma key-releases with dot key releases.
    """

    def keyPressEvent(self, event):      # pylint: disable=invalid-name
        """Replace commas with dots."""
        super().keyPressEvent(self.__make_dot_key_event(event))
    
    def keyReleaseEvent(self, event):    # pylint: disable=invalid-name
        """Replace commas with dots."""
        super().keyReleaseEvent(self.__make_dot_key_event(event))

    @staticmethod
    def __make_dot_key_event(event):
        """Return a KeyEvent with comma replaced by dot."""
        if event.key() != qtc.Qt.Key_Comma:
            return event
        return qtg.QKeyEvent(
            event.type(), qtc.Qt.Key_Period, event.modifiers(),
            text=event.text().replace(',', '.'),
            autorep=event.isAutoRepeat(), count=event.count()
            )


class InfIntSpinBox(qtw.QDoubleSpinBox):
    """A spin-box that allows inputting an infinite integer number."""

    def __init__(self, parent=None):
        """Initialize instance."""
        super().__init__(parent)
        self.setRange(0, float('inf'))
        self.setDecimals(0)