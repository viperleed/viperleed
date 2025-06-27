"""Module collapsibleviews of viperleed.gui.measure.widgets.

Defines CollapsibleView subclasses for displaying settings of
a single measurement device.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-07-05'
__license__ = 'GPLv3+'

from abc import abstractmethod
from pathlib import Path

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure.classes.abc import QMetaABC
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings
from viperleed.gui.measure.dialogs.settingsdialog import (
    SettingsDialogSectionBase,
    SettingsTag,
    )
from viperleed.gui.measure.hardwarebase import DEFAULTS_PATH
from viperleed.gui.measure.hardwarebase import make_device
from viperleed.gui.measure.hardwarebase import safe_connect
from viperleed.gui.measure.hardwarebase import safe_disconnect
from viperleed.gui.measure.widgets.pathselector import PathSelector
from viperleed.gui.measure.widgets.quantityselector import QuantitySelector
from viperleed.gui.widgets.collapsible import CollapsibleView


class CollapsibleDeviceView(CollapsibleView, metaclass=QMetaABC):
    """A CollapsibleView for measurement hardware."""

    # This signal is emitted when device settings stored in the
    # measurement settings file change. (These settings are the
    # quantities measured by the device and the settings file
    # the device takes its settings from.)
    settings_changed = qtc.pyqtSignal()

    def __init__(self, parent=None):
        """Initialise widget."""
        self._settings_folder = PathSelector(select_file=False, max_chars=25)
        self._settings_file_selector = qtw.QComboBox()
        self._device_cls = None
        self._device_info = None
        self._handler = None
        self._expanded = False
        self._original_settings = None
        super().__init__(parent=parent)

    @property
    def device_info(self):
        """Return the device info."""
        return self._device_info

    @property
    def has_hardware_interface(self):
        """Returns whether the view has a hardware interface or not."""
        return self.device_info and self.device_info.has_hardware_interface

    @property
    def original_settings(self):
        """Return the original settings.

        These settings are the settings the device has been loaded from.

        Returns
        -------
        Path or None
            None if the _original_settings has naver been set.
        """
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
        if not self.has_hardware_interface:
            return self._original_settings
        return self._settings_file_selector.currentData()

    @qtc.pyqtSlot(int)
    @qtc.pyqtSlot(bool)
    def set_expanded_state(self, expanded):
        """Collapse inner layout and remove handler when collapsed.

        Parameters
        ----------
        expanded : bool or int
            Converted to bool. Decides whether the collapsible
            widgets are visible/expanded. True means they become/stay
            visible.

        Returns
        -------
        None.
        """
        self._expanded = bool(expanded)
        if not self._expanded:
            self._remove_handler()
        else:
            self._try_fetching_handler()
        super().set_expanded_state(expanded)

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
        self._on_settings_folder_changed()

    def store_settings(self):
        """Store the edited settings."""
        if not self._handler:
            return
        try:
            self._handler.settings.update_file()
        except FileNotFoundError:
            # The file must have been moved before the settings could be saved.
            pass

    @property
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
    def _try_fetching_handler(self):
        """Check if creating a handler is possible."""
        enabled = self._expanded and self._settings_file_selector.isEnabled()
        if enabled:
            # Always fetch a settings handler for visible devices
            self._get_settings_handler()
        elif not self.has_hardware_interface:
            # Devices without a hardware interface are not enabled, as
            # they lack the ability to select a settings file, but they
            # are allowed to create a settings handler in order to
            # display the settings with which a previous measurement
            # was performed.
            self._get_settings_handler()
        else:
            # Disabled devices with a hardware interface were detected
            # but not selected by the user. We will fetch a handler in
            # the 'if enabled' case above whenever they are selected.
            pass

    def _check_if_settings_changed(self):
        """Check if the current settings differ from the original settings.

        Emits
        -----
        settings_changed
            If the device settings have changed.
        """
        if self.settings_file != self.original_settings:
            # Setting the original_settings is necessary here in order
            # to correctly emit the settings_changed signal if one
            # switches back and forth between two settings files.
            # Using _original_settings because we do not want to
            # trigger a _on_settings_folder_changed().
            self._original_settings = self.settings_file
            self.settings_changed.emit()

    def _compose_and_connect(self):
        """Compose and connect."""
        super()._compose_and_connect()
        self.add_collapsible_item(self._settings_folder)
        self.add_collapsible_item(self._settings_file_selector)
        self._settings_folder.path_changed.connect(
            self._on_settings_folder_changed
            )
        self._settings_file_selector.setEnabled(False)
        self._settings_file_selector.currentIndexChanged.connect(
            self._try_fetching_handler
            )

    @abstractmethod
    def _get_settings_handler(self):
        """Get the settings handler of the handled device.

        The reimplementation of this method must first check if it
        _can_make_device. If it can, then it must make_device and
        _make_handler_for_device. When the device has been connected,
        it must _update_widgets_from_device_settings and disconnect
        the device afterwards.

        Returns
        -------
        None.
        """

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
        self._handler = device.get_settings_handler()

    @qtc.pyqtSlot()
    def _on_settings_folder_changed(self):
        """Get settings files of the device handled by this view."""
        self._settings_file_selector.clear()
        self._settings_file_selector.setEnabled(False)
        if not self._device_cls or self._settings_folder.path == DEFAULTS_PATH:
            return

        if not self.has_hardware_interface:
            matching_settings = (() if not self._original_settings
                                 else (self._original_settings,))
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
            self._try_fetching_handler()

    @qtc.pyqtSlot()
    def _remove_handler(self):
        """Remove the internal reference to the handler and delete widgets."""
        if not self._handler:
            return
        layout = self._frame.layout()
        if layout.count() == 3:
            item = layout.itemAt(2)
            layout.removeItem(item)
            item.widget().deleteLater()
        self._handler = None

    @qtc.pyqtSlot()
    def _update_widgets_from_device_settings(self):
        """Populate the device view with the settings of the device."""
        if not self._handler:
            return
        self._handler.update_widgets()
        self._build_device_settings_widgets()
        self._check_if_settings_changed()


class CollapsibleCameraView(CollapsibleDeviceView):
    """A CollapsibleView for cameras."""

    def _connect_camera(self, camera):
        """Connect camera."""
        if not self.has_hardware_interface:
            return
        camera.connect_()
        camera.start()

    def _get_settings_handler(self):
        """Get the settings handler of the handled camera.

        Get and update settings handler widgets and
        populate the collapsible view with them.

        Returns
        -------
        None.
        """
        if not self._can_make_device:
            return
        camera = make_device(self.settings_file, self._device_cls,
                             self._device_info)
        try:
            self._connect_camera(camera)
        except camera.exceptions:
            pass
        self._make_handler_for_device(camera)
        self._update_widgets_from_device_settings()
        camera.disconnect_()


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
        self.loading_measurement_settings = False
        self._is_primary = False
        self._trying_to_get_settings = False
        self._primary_changed = False
        self._quantity_selector = None
        self._quantities_to_set = ()
        super().__init__(parent=parent)

    @property
    def selected_quantities(self):
        """Return the selected quantities."""
        if not self._quantity_selector:
            return tuple()
        return self._quantity_selector.get_selected_quantities()

    @qtc.pyqtSlot(bool)
    def set_primary(self, is_primary):
        """Set whether the controller is the primary controller.

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
            self._try_fetching_handler()

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
        self._quantity_selector.settings_changed.connect(
            self._check_if_quantities_changed, type=qtc.Qt.UniqueConnection
            )

    @qtc.pyqtSlot()
    def _build_device_settings_widgets(self):
        """Get the handler and quantity widgets for the device settings."""
        if not self.has_hardware_interface:
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
        if not self.loading_measurement_settings:
            self.settings_changed.emit()
        else:
            self.loading_measurement_settings = False

    @qtc.pyqtSlot()
    def _check_if_primary_changed(self):
        """Check if the primary status changed while getting settings."""
        self._trying_to_get_settings = False
        if self._primary_changed:
            self._primary_changed = False
            self._try_fetching_handler()

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
        """Get the settings handler of the handled controller.

        Get and update settings handler widgets and
        populate the collapsible view with them.

        Returns
        -------
        None.
        """
        if not self._can_make_device:
            return
        ctrl = make_device(self.settings_file, self._device_cls,
                           self._device_info, sets_energy=self._is_primary)
        self._trying_to_get_settings = True
        self._make_handler_for_device(ctrl)
        if not self.has_hardware_interface:
            self._trying_to_get_settings = False
            self._update_widgets_from_device_settings()
            return
        ctrl.ready_to_show_settings.connect(
            self._update_widgets_from_device_settings
            )
        ctrl.ready_to_show_settings.connect(ctrl.disconnect_)
        ctrl.ready_to_show_settings.connect(self._check_if_primary_changed)
        ctrl.prepare_to_show_settings()
