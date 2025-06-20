"""Module spinboxes of viperleed.gui.measure.widgets.

Defines convenience subclasses of QSpinBox and QDoubleSpinBox.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-12-06'
__license__ = 'GPLv3+'

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import QtWidgets as qtw


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
