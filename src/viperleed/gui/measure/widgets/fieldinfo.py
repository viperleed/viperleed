"""Module fieldinfo of viperleed.gui.measure.widgets.

Defines the FieldInfo class: a QLabel with icon, used for displaying
extra information in a tooltip. Less conspicuous than just having a
tooltip appear on mouse hover.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-09-27'
__license__ = 'GPLv3+'


import html

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import QtWidgets as qtw


class FieldInfo(qtw.QLabel):
    """A simple widget for displaying information in a tooltip."""

    def __init__(self, tooltip='', size=16, parent=None):
        """Initialize instance."""
        super().__init__(parent=parent)
        _style = self.style()
        _icon = _style.standardIcon(_style.SP_FileDialogInfoView)
        self.setPixmap(_icon.pixmap(size))
        self.set_info_text(tooltip)

        # Timer is used to show tooltip on mouse press
        self.__tooltip_timer = qtc.QTimer(self)
        self.__tooltip_timer.setSingleShot(True)
        self.__tooltip_timer.setInterval(150)
        self.__tooltip_timer.timeout.connect(self.__show_tooltip)

    def mouseDoubleClickEvent(self, event):  # pylint: disable=invalid-name
        """Show info on mouse double click."""
        timer = self.__tooltip_timer
        if not qtw.QToolTip.isVisible() and not timer.isActive():
            timer.start()
        event.accept()

    def mouseMoveEvent(self, event):     # pylint: disable=invalid-name
        """Catch mouse move, as we also do with press."""
        event.accept()

    def mousePressEvent(self, event):    # pylint: disable=invalid-name
        """Show info on mouse press."""
        timer = self.__tooltip_timer
        if not qtw.QToolTip.isVisible() and not timer.isActive():
            timer.start()
        event.accept()

    def mouseReleaseEvent(self, event):  # pylint: disable=invalid-name
        """Catch mouse release, as we also do with press."""
        event.accept()

    def setText(self, text):             # pylint: disable=invalid-name
        """Set info text rather than displayed text."""
        self.set_info_text(text)

    def set_info_text(self, text):
        """Set the information text to be shown."""
        if text and not qtc.Qt.mightBeRichText(text):
            text = f"<qt>{html.escape(text)}<qt/>"
        self.setToolTip(text)

    def text(self):
        """Return the info text."""
        return self.toolTip()

    @qtc.pyqtSlot()
    def __show_tooltip(self):
        """Show tooltip for an appropriate amount of time."""
        if qtw.QToolTip.isVisible():
            return
        qtw.QToolTip.showText(qtg.QCursor.pos(), self.toolTip())
