"""Module buttons of viperleed.gui.widgets.

Collects custom buttons and button-related widgets.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-06-13'
__license__ = 'GPLv3+'

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw


class ButtonWithLabel(qtw.QWidget):
    """QLabel and QPushButton."""

    def __init__(self, tight=True, **kwargs):
        """Initialize widget.

        Parameters
        ----------
        tight : bool, optional
            Whether there should be no extra space surrounding
            this widget. Default is True.
        **kwargs : object, optional
            Other optional arguments. Passed unaltered to the
            parent class.

        Returns
        -------
        None.
        """
        super().__init__(**kwargs)
        self.label = qtw.QLabel()
        self.button = QNoDefaultPushButton()
        self._compose(tight)

    def _compose(self, tight):
        """Place children widgets."""
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


class QNoDefaultDialogButtonBox(qtw.QDialogButtonBox):
    """QDialogButtonBox without default buttons."""

    def event(self, event):
        """Overwrite event to skip setting default."""
        if event.type() == qtc.QEvent.Show:
            self._unset_default_buttons()
            return qtw.QWidget().event(event)
        return super().event(event)

    def _unset_default_buttons(self):
        """Ensure no button is a default."""
        for button in self.buttons():
            button.setAutoDefault(False)
            button.setDefault(False)


class QNoDefaultPushButton(qtw.QPushButton):
    """QPushbutton that is not a default button."""

    def __init__(self, *args, **kwargs):
        """Initialize button."""
        super().__init__(*args, **kwargs)
        self.setAutoDefault(False)
        self.setDefault(False)


class QNoDefaultIconButton(QNoDefaultPushButton):
    """QNoDefaultPushButton with right-aligned icons."""

    def __init__(self, *args, **kwargs):
        """Initialize button."""
        super().__init__(*args, **kwargs)
        self.setLayout(qtw.QHBoxLayout())
        self._label = qtw.QLabel()
        self._label.setAlignment(qtc.Qt.AlignRight | qtc.Qt.AlignVCenter)
        self._label.setAttribute(qtc.Qt.WA_TransparentForMouseEvents, True)
        self.layout().addWidget(self._label)

    def setIcon(self, icon):        # pylint: disable=invalid-name
        """Set the desired icon on the button.

        Parameters
        ----------
        icon : QIcon
            The icon to be added to the button.

        Returns
        -------
        None.
        """
        icon_size_real = icon.actualSize(self.sizeHint()*0.5)
        self._label.setPixmap(icon.pixmap(icon_size_real))


class QUncheckableButtonGroup(qtw.QButtonGroup):
    """QButtonGroup that can be unchecked."""

    def __init__(self, *args, **kwargs):
        """Initialize button group."""
        super().__init__(*args, **kwargs)
        self._pressed_button_was_checked = False
        self.buttonPressed.connect(self._remember_checked_state_of_button)

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
    def _remember_checked_state_of_button(self, button):
        if button == self.checkedButton():
            self._pressed_button_was_checked = True
