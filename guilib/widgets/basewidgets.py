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
        self.button = QNoDefaultPushButton()
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


class QCheckBoxInvertedSignal(qtw.QCheckBox):
    """QCheckBox with extra unchecked signal."""

    unchecked = qtc.pyqtSignal(bool)

    def __init__(self, **kwargs):
        """Initialise widget."""
        super().__init__(**kwargs)
        self.stateChanged.connect(self.emit_inverted_signal)

    @qtc.pyqtSlot(int)
    def emit_inverted_signal(self, value):
        """Emit unchecked signal."""
        self.unchecked.emit(not bool(value))


class QNoDefaultDialogButtonBox(qtw.QDialogButtonBox):
    """QDialogButtonBox without default button."""

    def event(self, event):
        """Overwrite event to skip setting default."""
        if event.type() == qtc.QEvent.Show:
            self._unset_default_buttons()
            return qtw.QWidget().event(event)
        return super().event(event)

    def _unset_default_buttons(self):
        """Ensure no buttons is a default."""
        for button in self.buttons():
            button.setAutoDefault(False)
            button.setDefault(False)


class QNoDefaultPushButton(qtw.QPushButton):
    """QPushbutton that is not a default button."""

    def __init__(self, *args, **kwargs):
        """Initialise button."""
        super().__init__(*args, **kwargs)
        self.setAutoDefault(False)


class QNoDefaultIconButton(QNoDefaultPushButton):
    """QNoDefaultPushButton with right-aligned icons."""

    def __init__(self, *args, **kwargs):
        """Initialise button."""
        super().__init__(*args, **kwargs)
        self.setLayout(qtw.QHBoxLayout())
        self._label = qtw.QLabel()
        self._label.setAlignment(qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)
        self._label.setAttribute(qtc.Qt.WA_TransparentForMouseEvents, True)
        self.layout().addWidget(self._label)

    def setIcon(self, icon):
        """Set the desired icon on the button.

        Parameters
        ----------
        icon : QIcon
            The icon to be added to the button.
        """
        icon_size_real = icon.actualSize(self.sizeHint()*0.5)
        self._label.setPixmap(icon.pixmap(icon_size_real))


class QUncheckableButtonGroup(qtw.QButtonGroup):
    """QButtonGroup that can be unchecked."""

    def __init__(self, *args, **kwargs):
        """Initialise button group."""
        super().__init__(*args, **kwargs)
        self._pressed_button_was_checked = False
        self.buttonPressed.connect(self._remember_status_of_pressed_button)

    def uncheck_buttons(self):
        """Uncheck all buttons."""
        self.setExclusive(False)
        self.checkedButton().setChecked(False)
        self.setExclusive(True)

    @qtc.pyqtSlot(qtw.QAbstractButton)
    def uncheck_if_clicked(self, button):
        """Uncheck if checked button was clicked."""
        was_checked = self._pressed_button_was_checked
        self._pressed_button_was_checked = False
        if was_checked and button == self.checkedButton():
            self.uncheck_buttons()

    @qtc.pyqtSlot(qtw.QAbstractButton)
    def _remember_status_of_pressed_button(self, button):
        if button == self.checkedButton():
            self._pressed_button_was_checked = True
