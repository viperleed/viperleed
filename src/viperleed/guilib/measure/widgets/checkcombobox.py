"""Module checkcombobox of viperleed.guilib.measure.widgets.

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2021-12-13
Author: Michele Riva
Author: Florian Doerr

Defines the CheckComboBox class.
"""

from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw)


class CheckComboBox(qtw.QComboBox):
    """A QComboBox that allows selecting multiple options."""

    check_changed = qtc.pyqtSignal()

    def __init__(self, parent=None):
        """Initialize instance.

        Parameters
        ----------
        parent : QWidget
            The parent widget of this combo box.

        Returns
        -------
        None.
        """
        super().__init__(parent=parent)
        self.view().pressed.connect(self.__ensure_changed)
        self.textActivated.connect(self.toggle_checked)
        self.__changed = False

    @staticmethod
    def is_item_checked(item):
        """Return whether an item is checked."""
        return item.checkState() == qtc.Qt.Checked

    @property
    def selected_items(self):
        """Return a list with the text of checked items."""
        checked = []
        for i in range(self.count()):
            item = self.model().item(i, 0)
            if self.is_item_checked(item):
                checked.append(item.text())
        return checked

    def addItem(self, name):            # pylint: disable=invalid-name
        """Add an item to the combo box.

        Extend base-class behavior to make the new item checkable.

        Parameters
        ----------
        name : str
            The new entry to be added.

        Returns
        -------
        None.
        """
        super().addItem(name)
        self.model().findItems(name)[0].setCheckState(qtc.Qt.Unchecked)

    def addItems(self, names):          # pylint: disable=invalid-name
        """Add multiple items to the combo box.

        Extend base-class behavior to make the new items checkable.

        Parameters
        ----------
        names : Sequence
            The new entries to be added. Each element should be
            a string.

        Returns
        -------
        None.
        """
        super().addItems(names)
        for name in names:
            self.model().findItems(name)[0].setCheckState(qtc.Qt.Unchecked)

    def hidePopup(self):                # pylint: disable=invalid-name
        """Hide the item list.

        Extend base-class behavior to hide the pop-up only if
        the user clicked on the outside of the combo box (or,
        more generally, if the click did not lead to a change
        of the checked state of options).

        Returns
        -------
        None.
        """
        if not self.__changed:
            super().hidePopup()
        self.__changed = False

    def insertItem(self, index, name):  # pylint: disable=invalid-name
        """Insert one item at a specific position.

        Extend base-class behavior to make the new item checkable.

        Parameters
        ----------
        index : int
            The position at which the new item should be inserted.
            All items already present in the combo box at indices
            >= index will be moved by one index "down". If the index
            is equal to or higher than the total number of items,
            the new item is appended to the list of existing items.
            If the index is zero or negative, the new item is
            prepended to the list of existing items.
        name : str
            The new entry to be added. Each element should be
            a string.

        Returns
        -------
        None.
        """
        super().insertItem(index, name)
        self.model().findItems(name)[0].setCheckState(qtc.Qt.Unchecked)

    # pylint: disable=invalid-name
    def insertItems(self, index, names):
        """Insert multiple items starting at a specific position.

        Extend base-class behavior to make the new items checkable.

        Parameters
        ----------
        index : int
            The position at which the new items should be inserted.
            All items already present in the combo box at indices
            >= index will be moved "down" by len(names). If the index
            is equal to or higher than the total number of items,
            the new items are appended to the list of existing items.
            If the index is zero or negative, the new items are
            prepended to the list of existing items.
        names : Sequence
            The new entries to be added. Each element should be
            a string.

        Returns
        -------
        None.
        """
        super().insertItem(index, names)
        for name in names:
            self.model().findItems(name)[0].setCheckState(qtc.Qt.Unchecked)
    # pylint: enable=invalid-name

    def toggle_checked(self, name):
        """Toggle the checked state of an entry with given name."""
        item = self.model().findItems(name)[0]
        self.__changed = True
        if item.checkState() == qtc.Qt.Checked:
            item.setCheckState(qtc.Qt.Unchecked)
        else:
            item.setCheckState(qtc.Qt.Checked)
        self.check_changed.emit()

    def __ensure_changed(self):
        """Set __changed flag to True.

        This is the slot connected to the view().pressed signal.
        This connection is necessary to ensure that the drop-down
        item list is properly shown/hidden, especially when
        multiple clicks happen in rapid succession.

        Returns
        -------
        None.
        """
        self.__changed = True

    def paintEvent(self, __event):       # pylint: disable=invalid-name
        """Render this combo box on screen.

        Reimplement base-class behavior to display the list
        of checked items in the 'header' line.

        Returns
        -------
        None.
        """
        painter = qtw.QStylePainter(self)
        painter.setPen(self.palette().color(self.palette().Text))
        opt = qtw.QStyleOptionComboBox()
        self.initStyleOption(opt)
        opt.currentText = ", ".join(self.selected_items)
        if not opt.currentText:
            opt.currentText = '-- Select option(s) --'
        painter.drawComplexControl(qtw.QStyle.CC_ComboBox, opt)
        painter.drawControl(qtw.QStyle.CE_ComboBoxLabel, opt)

    def uncheck_all(self):
        """Uncheck all checked options."""
        for i in range(self.count()):
            item = self.model().item(i, 0)
            if self.is_item_checked(item):
                item.setCheckState(qtc.Qt.Unchecked)
            item.setEnabled(True)
