"""Module uimeasurement of viperleed.guilib.measure

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-12-13
Author: Michele Riva
Author: Florian Doerr

Defines the CheckComboBox class.
"""

import PyQt5.QtCore as qtc
import PyQt5.QtWidgets as qtw
import PyQt5.QtGui as qtg


class CheckComboBox(qtw.QComboBox):
    """A ComboBox able to select multiple choices."""

    def __init__(self, parent=None):
        """Initialise class."""
        super().__init__(parent=parent)
        self.view().pressed.connect(self.ensure_change)
        self.textActivated.connect(self.change_checked)
        self.__changed = False

    def addItem(self, name):
        super().addItem(name)
        self.model().findItems(name)[0].setCheckState(qtc.Qt.Unchecked)

    def addItems(self, names):
        super().addItems(names)
        for name in names:
            self.model().findItems(name)[0].setCheckState(qtc.Qt.Unchecked)

    def insertItem(self, index, name):
        super().insertItem(index, name)
        self.model().findItems(name)[0].setCheckState(qtc.Qt.Unchecked)

    def insertItems(self, index, names):
        super().insertItem(index, names)
        for name in names:
            self.model().findItems(name)[0].setCheckState(qtc.Qt.Unchecked)

    def change_checked(self, name):
        item = self.model().findItems(name)[0]
        self.__changed = True
        if item.checkState() == qtc.Qt.Checked:
            item.setCheckState(qtc.Qt.Unchecked)
        else:
            item.setCheckState(qtc.Qt.Checked)

    def hidePopup(self):
        if not self.__changed:
            super().hidePopup()
        self.__changed = False

    def showPopup(self):
        super().showPopup()

    def ensure_change(self):
        self.__changed = True

    def is_item_checked(self, item):
        return item.checkState() == qtc.Qt.Checked

    def get_items(self):
        checkedItems = []
        for i in range(self.count()):
            item = self.model().item(i, 0)
            if self.is_item_checked(item):
                checkedItems.append(item.text())
        return checkedItems

    def paintEvent(self, event):
        painter = qtw.QStylePainter(self)
        painter.setPen(self.palette().color(qtg.QPalette.Text))
        opt = qtw.QStyleOptionComboBox()
        self.initStyleOption(opt)
        opt.currentText = ", ".join(self.get_items())
        if not opt.currentText:
            opt.currentText = '<Select option(s)>'
        painter.drawComplexControl(qtw.QStyle.CC_ComboBox, opt)
        painter.drawControl(qtw.QStyle.CE_ComboBoxLabel, opt)
