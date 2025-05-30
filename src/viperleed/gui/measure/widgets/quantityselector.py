"""Module quantityselector of viperleed.gui.measure.widgets.

Defines the QuantitySelector class.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-26'
__license__ = 'GPLv3+'

from ast import literal_eval

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure.widgets.fieldinfo import FieldInfo
from viperleed.gui.widgets.basewidgets import QUncheckableButtonGroup


class QuantitySelector(qtw.QFrame):
    """A widget that allows selection of quantities from settings."""

    # This signal is emitted when the quantity selection changes.
    settings_changed = qtc.pyqtSignal()

    def __init__(self, settings, parent=None):
        """Initialise widget.

        Parameters
        ----------
        settings : ViPErLEEDSettings
            The settings of the associated controller. Contains the
            quantities the controller can measure under 'controller',
            'measurement_devices'. Each tuple stands for one
            measurement device and for each tuple, only one quantity
            can be selected at the same time.
        parent : QObject
            The parent QObject of this widget.

        Returns
        -------
        None.
        """
        super().__init__(parent=parent)
        self.setMidLineWidth(1)
        self.setFrameStyle(self.Panel | self.Raised)
        self._selections = []
        self._quantities = literal_eval(
            settings.get('controller', 'measurement_devices', fallback=())
            )
        self._compose()

    def _compose(self):
        """Compose quantity widgets."""
        main_layout = qtw.QVBoxLayout()
        label_layout = qtw.QHBoxLayout()
        label = qtw.QLabel('Measured Quantities')
        size = label.fontMetrics().boundingRect('a').height()
        info = ('<nobr>Quantities in the same column cannot </nobr>'
                'be measured at the same time.')
        label_layout.addWidget(label)
        label_layout.addWidget(FieldInfo(info, size=size))
        label_layout.addStretch(1)
        main_layout.addLayout(label_layout)
        quantity_layout = qtw.QHBoxLayout()
        for selections in self._quantities:
            layout = qtw.QVBoxLayout()
            group = QUncheckableButtonGroup()
            group.buttonClicked.connect(group.uncheck_if_clicked)
            group.buttonClicked.connect(self.settings_changed)
            for quantity in selections:
                button = qtw.QCheckBox()
                group.addButton(button)
                button.setText(quantity)
                layout.addWidget(button)
            layout.addStretch(1)
            self._selections.append(group)
            quantity_layout.addLayout(layout)
            if selections != self._quantities[-1]:
                line = qtw.QFrame()
                line.setFrameShape(qtw.QFrame.VLine)
                quantity_layout.addWidget(line)
        main_layout.addLayout(quantity_layout)
        self.setLayout(main_layout)

    def get_selected_quantities(self):
        """Return the selected quantities.

        Returns
        -------
        quantities : tuple
            The selected quantities.
        """
        quantities = []
        for group in self._selections:
            if group.checkedButton():
                quantities.append(group.checkedButton().text())
        return tuple(quantities)

    def set_quantities(self, quantities):
        """Set quantities from settings.

        Parameters
        ----------
        quantities : tuple
            The quantities to set.

        Returns
        -------
        None.
        """
        buttons = tuple(btn for group in self._selections
                        for btn in group.buttons())
        for quantity in quantities:
            for button in buttons:
                if button.text() == quantity:
                    button.click()
                    break
            else:
                raise ValueError(f'{quantity!r} is not an '
                                 'acceptable quantity.')
