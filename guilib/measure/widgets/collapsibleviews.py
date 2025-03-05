"""Module collapsibleviews of viperleed.guilib.measure.widgets.

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2024-07-05
Author: Michele Riva
Author: Florian Doerr

Defines various CollapsibleView classes and other classes used in them.
"""

from ast import literal_eval
from pathlib import Path

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure.classes.settings import ViPErLEEDSettings
from viperleed.guilib.measure.dialogs.settingsdialog import (
    SettingsDialogSectionBase,
    SettingsTag,
    )
from viperleed.guilib.measure.hardwarebase import make_device
from viperleed.guilib.measure.hardwarebase import safe_connect
from viperleed.guilib.measure.hardwarebase import safe_disconnect
from viperleed.guilib.measure.widgets.pathselector import PathSelector
from viperleed.guilib.widgets.basewidgets import CollapsibleView
from viperleed.guilib.widgets.basewidgets import QUncheckableButtonGroup


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
        label = qtw.QLabel('Measured Quantities')
        main_layout.addWidget(label)
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
            button_not_found = True
            for button in buttons:
                if button.text() == quantity:
                    button.click()
                    button_not_found = False
                    break
            if button_not_found:
                raise ValueError(f'{quantity!r} is not an acceptable quantity.')


class CollapsibleDeviceView(CollapsibleView):
    """A CollapsibleView for measurement hardware."""

    # This signal is emitted when device settings stored
    # in the measurement settings file change.
    # (Measured quantities and the device settings file.)
    settings_changed = qtc.pyqtSignal()

    def __init__(self, parent=None):
        """Initialise widget.

        parent : QObject
            The parent QObject of this widget.

        Returns
        -------
        None.
        """
        self._settings_folder = PathSelector(select_file=False, max_chars=25)
        self._settings_file_selector = qtw.QComboBox()
        super().__init__(parent=parent)
        self._device_cls = None
        self._device_info = None
        self._handler = None
        self._enabled = False
        self._original_settings = None

    @property
    def device_info(self):
        """Return the device info."""
        return self._device_info

    @property
    def original_settings(self):
        """Return the original settings."""
        return self._original_settings

    @original_settings.setter
    def original_settings(self, settings):
        """Store a path to settings as the original ones."""
        self._original_settings = settings
        self.set_settings_folder(settings.parent)
        self._settings_file_selector.setCurrentText(settings.stem)

    @property
    def settings_file(self):
        """Return the path to the selected settings file."""
        if self.is_dummy_device():
            return self._original_settings
        return self._settings_file_selector.currentData()

    @qtc.pyqtSlot()
    def _build_device_settings(self):
        """Populate the device view with the settings of the device."""
        if not self._handler:
            return
        self._handler.update_widgets()
        self._build_device_settings_widgets()
        self._check_if_settings_changed()

    def _build_device_settings_widgets(self):
        """Get the handler widgets for the device settings.

        Places all widgets in the SettingsHandler with the
        SettingsTag.MEASUREMENT into the collapsible QFrame.

        Returns
        -------
        None.
        """
        new_widget = qtw.QWidget()
        layout = qtw.QVBoxLayout()
        measurement_widgets = self._handler.get_widgets_with_tags(
            SettingsTag.MEASUREMENT
            )
        for widget in measurement_widgets:
            if isinstance(widget, SettingsDialogSectionBase):
                layout.addWidget(widget)
                for option in widget.options:
                    option.setVisible(True)
                continue
            widget.setVisible(True)
            _form = qtw.QFormLayout()
            _form.addRow(*widget)
            layout.addLayout(_form)
        new_widget.setLayout(layout)
        self.add_collapsible_item(new_widget)
        self._adjust_bottom_space()

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(int)
    def _check_for_handler(self):
        """Check if creating a handler is possible."""
        if not self._enabled or not self._settings_file_selector.isEnabled():
            if not self.is_dummy_device():
                return
        self._get_settings_handler()

    def _check_if_settings_changed(self):
        """Check if the current settings differ from the original settings.

        Emits
        -----
        settings_changed
            If the device settings have changed.
        """
        if self.settings_file != self._original_settings:
            self._original_settings = self.settings_file
            self.settings_changed.emit()

    def _compose_and_connect(self):
        """Compose and connect."""
        super()._compose_and_connect()
        self._settings_folder.path = Path().resolve()
        self.add_collapsible_item(self._settings_folder)
        self.add_collapsible_item(self._settings_file_selector)
        self._settings_folder.path_changed.connect(
            self._get_device_settings_files
            )
        self._settings_file_selector.setEnabled(False)
        self._settings_file_selector.currentIndexChanged.connect(
            self._check_for_handler
            )

    def _can_make_device(self):
        """Check if required settings to make a device are present.

        Returns
        -------
        can_make_device : bool
            True if making a device is possible.
        """
        return bool(self.settings_file
                    and self._device_cls
                    and self._device_info)

    @qtc.pyqtSlot()
    def _get_device_settings_files(self):
        """Get settings files of the device handled by this view."""
        self._settings_file_selector.clear()
        self._settings_file_selector.setEnabled(False)
        if not self._device_cls or self._settings_folder.path == Path():
            return

        if self.is_dummy_device():
            if self._original_settings:
                matching_settings = (self._original_settings,)
            else:
                matching_settings = ()
        else:
            matching_settings = self._device_cls.find_matching_settings_files(
                self._device_info, self._settings_folder.path, False,
                )

        for settings in matching_settings:
            self._settings_file_selector.addItem(settings.stem,
                                                 userData=settings)
        self._remove_handler()
        if any(matching_settings):
            self._settings_file_selector.setEnabled(True)
            self._check_for_handler()

    def is_dummy_device(self):
        """Returns whether the view is a dummy object or not."""
        return not self._device_info.has_hardware_interface

    def _make_handler_for_device(self, device):
        """Internally store a SettingsHandler for `device`.

        Parameters
        ----------
        device : DeviceABC
            A device that can return a SettingsHandler.

        Returns
        -------
        None.
        """
        if self._handler:
            self._remove_handler()
        device.uses_default_settings = False
        self._handler = device.get_settings_handler()

    @qtc.pyqtSlot()
    def _remove_handler(self):
        """Remove the internal reference to the handler."""
        if not self._handler:
            return
        layout = self._frame.layout()
        if layout.count() == 3:
            item = layout.itemAt(2)
            layout.removeItem(item)
            item.widget().deleteLater()
        self._handler = None

    @qtc.pyqtSlot(int)
    @qtc.pyqtSlot(bool)
    def enable_view(self, enable):
        """Collapse inner layout and remove handler when disabled.

        Parameters
        ----------
        enable : bool or int
            Converted to bool. Decides whether the collapsible
            widgets are visible/expanded. True means they become/stay
            visible.

        Returns
        -------
        None.
        """
        self._enabled = bool(enable)
        if not self._enabled:
            self._remove_handler()
        else:
            self._check_for_handler()
        super().enable_view(enable)

    def set_device(self, device_cls, device_info):
        """Set the device for which this view is meant.

        Parameters
        ----------
        device_cls : type
            The class of the device.
        device_info : SettingsInfo
            The SettingsInfo necessary to determine the settings.

        Returns
        -------
        None.
        """
        self._device_cls = device_cls
        self._device_info = device_info

    def set_settings_folder(self, settings_folder_path):
        """Set the path of the settings folder to the given path.

        Parameters
        ----------
        settings_folder_path : Path or str
            The path to the folder containing the settings.

        Returns
        -------
        None.
        """
        self._settings_folder.path = settings_folder_path

    def store_settings(self):
        """Store the edited settings."""
        if not self._handler:
            return
        try:
            self._handler.settings.update_file()
        except FileNotFoundError:
            # The file must have been moved before the settings could be saved.
            pass


class CollapsibleCameraView(CollapsibleDeviceView):
    """A CollapsibleView for cameras."""

    def _connect_camera(self, device):
        """Connect camera.

        Parameters
        ----------
        device : CameraABC
            A camera class instance.

        Returns
        -------
        None.
        """
        if self.is_dummy_device():
            return
        device.connect_()
        device.start()

    def _get_settings_handler(self):
        """Get the settings handler of the handled device."""
        if not self._can_make_device():
            return
        device = make_device(self.settings_file, self._device_cls,
                             self._device_info)
        try:
            self._connect_camera(device)
        except device.exceptions:
            pass
        self._make_handler_for_device(device)
        self._build_device_settings()
        device.disconnect_()


class CollapsibleControllerView(CollapsibleDeviceView):
    """A CollapsibleView for controllers."""

    def __init__(self, parent=None):
        """Initialise widget.

        Parameters
        ----------
        parent : QObject
            The parent QObject of this widget.

        Returns
        -------
        None.
        """
        super().__init__(parent=parent)
        self._is_primary = False
        self._trying_to_get_settings = False
        self._primary_changed = False
        self._quantity_selector = None
        self._quantities_to_set = ()

    def set_quantities(self, quantities):
        """Set the quantities to measure.

        Parameters
        ----------
        quantities : tuple
            The quantities to set.

        Returns
        -------
        None.
        """
        self._quantities_to_set = quantities
        if not self._quantity_selector:
            return
        safe_disconnect(self._quantity_selector.settings_changed,
                        self._check_if_quantities_changed)
        self._quantity_selector.set_quantities(quantities)
        safe_connect(self._quantity_selector.settings_changed,
                     self._check_if_quantities_changed,
                     type=qtc.Qt.UniqueConnection)

    @property
    def selected_quantities(self):
        """Return the selected quantities."""
        if not self._quantity_selector:
            return tuple()
        return self._quantity_selector.get_selected_quantities()

    @qtc.pyqtSlot()
    def _build_device_settings_widgets(self):
        """Get the handler and quantity widgets for the device settings."""
        if self.is_dummy_device():
            settings = ViPErLEEDSettings.from_settings(self.original_settings)
        else:
            settings = self.sender().settings
        super()._build_device_settings_widgets()
        # Get the layout in which the SettingsHandler widgets are contained.
        layout = self._frame.layout().itemAt(2).widget().layout()
        self._quantity_selector = QuantitySelector(settings)
        if self._quantities_to_set:
            self._quantity_selector.set_quantities(self._quantities_to_set)
        safe_connect(self._quantity_selector.settings_changed,
                     self._check_if_quantities_changed,
                     type=qtc.Qt.UniqueConnection)
        layout.addWidget(self._quantity_selector)

    def _check_if_primary_changed(self):
        """Check if the primary status changed while getting settings."""
        self._trying_to_get_settings = False
        if self._primary_changed:
            self._primary_changed = False
            self._check_for_handler()

    @qtc.pyqtSlot()
    def _check_if_quantities_changed(self):
        """Check if the current settings differ from the original settings.

        Emits
        -----
        settings_changed
            If the selected quantities to measure have changed.
        """
        if self._quantities_to_set != self.selected_quantities:
            self._quantities_to_set = self.selected_quantities
            self.settings_changed.emit()

    def _get_settings_handler(self):
        """Get the settings handler of the handled device."""
        if not self._can_make_device():
            return
        device = make_device(self.settings_file, self._device_cls,
                             self._device_info, sets_energy=self._is_primary)
        self._trying_to_get_settings = True
        self._make_handler_for_device(device)
        if self.is_dummy_device():
            self._trying_to_get_settings = False
            self._build_device_settings()
            return
        device.ready_to_show_settings.connect(self._build_device_settings)
        device.ready_to_show_settings.connect(device.disconnect_)
        device.ready_to_show_settings.connect(self._check_if_primary_changed)
        device.prepare_to_show_settings()

    @qtc.pyqtSlot(bool)
    def set_primary(self, is_primary):
        """Set whether the device is the primary controller.

        Parameters
        ----------
        is_primary : bool
            Whether this controller is supposed to be the primary
            controller. If true, the controller will be treated as the
            primary controller and its settings which are relevant for
            setting energies will be displayed in the collapsible view.

        Returns
        -------
        None.
        """
        was_primary = self._is_primary
        self._is_primary = is_primary
        if was_primary == is_primary:
            return
        if self._trying_to_get_settings:
            self._primary_changed = True
        else:
            self._check_for_handler()
