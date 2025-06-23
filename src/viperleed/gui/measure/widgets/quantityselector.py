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

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure.classes.datapoints import QuantityInfo
from viperleed.gui.measure.widgets.fieldinfo import FieldInfo
from viperleed.gui.widgets.buttons import QUncheckableButtonGroup


class QuantitySelector(qtw.QGroupBox):
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
        self.setFlat(True)
        self.setTitle('Measured Quantities')
        self._selections = {}
        self._groups = []
        quantities = self._parse_quantities(settings.getsequence(
            'controller', 'measurement_devices', fallback=()
            ))
        self._compose(quantities)

    def get_selected_quantities(self):
        """Return the selected quantities.

        Returns
        -------
        quantities : tuple
            The selected quantities.
        """
        return tuple(q.label for q, btn in self._selections.items()
                     if btn.isChecked())

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
        quantities = [QuantityInfo.from_label(name) for name in quantities]
        for quantity in quantities:
            try:
                button = self._selections[quantity]
            except KeyError:
                raise ValueError(f'{quantity.label!r} is not '
                                 'an acceptable quantity.') from None
            button.click()

    def _compose(self, quantities):
        """Compose quantity widgets.

        Parameters
        ----------
        quantities : Sequence
            Elements are squences of QuantityInfo objects. Each element
            represents one measuring device installed on the controller,
            with each QuantityInfo being one of the quantities that
            measuring device in particular can acquire.

        Returns
        -------
        None.
        """
        main_layout = qtw.QVBoxLayout()
        quantity_layout = self._make_quantity_layout(quantities)
        main_layout.addLayout(quantity_layout)
        main_layout.addItem(qtw.QSpacerItem(0, 5))
        label = qtw.QLabel('Quantities in the same column cannot '
                           'be measured at the same time.')
        label.setWordWrap(True)
        main_layout.addWidget(label)
        self.setLayout(main_layout)

    def _make_quantity_layout(self, quantities):
        """Return a layout containing the quantity selection.

        Parameters
        ----------
        quantities : Sequence
            Elements are squences of QuantityInfo objects. Each element
            represents one measuring device installed on the controller,
            with each QuantityInfo being one of the quantities that
            measuring device in particular can acquire.

        Returns
        -------
        None.
        """
        quantity_layout = qtw.QHBoxLayout()
        for selections in quantities:
            layout = qtw.QVBoxLayout()
            group = QUncheckableButtonGroup()
            group.buttonClicked.connect(group.uncheck_if_clicked)
            group.buttonClicked.connect(self.settings_changed)
            # Without a reference to the groups, the groups would be deleted.
            self._groups.append(group)
            for quantity in selections:
                btn_layout = qtw.QHBoxLayout()
                button = qtw.QCheckBox()
                group.addButton(button)
                button.setText(quantity.display_label)
                btn_layout.addWidget(button)
                tip = quantity.description
                btn_layout.addWidget(FieldInfo.for_widget(button, tooltip=tip))
                btn_layout.addStretch(1)
                layout.addLayout(btn_layout)
                self._selections[quantity] = button
            layout.addStretch(1)
            quantity_layout.addLayout(layout)
            if selections != quantities[-1]:
                line = qtw.QFrame()
                line.setFrameShape(qtw.QFrame.VLine)
                quantity_layout.addWidget(line)
        return quantity_layout

    @staticmethod
    def _parse_quantities(quantity_name_by_device):
        return [[QuantityInfo.from_label(name) for name in device]
                for device in quantity_name_by_device]
