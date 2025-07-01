"""Module collapsiblelists of viperleed.gui.measure.widgets.

Defines CollapsibleList subclasses for displaying the settings of
measurement devices.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-14'
__license__ = 'GPLv3+'

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure.classes.abc import NO_HARDWARE_INTERFACE
from viperleed.gui.measure.classes.abc import QObjectSettingsErrors
from viperleed.gui.measure.classes.abc import SettingsInfo
from viperleed.gui.measure.classes.decorators import emit_default_faulty
from viperleed.gui.measure.classes.settings import DefaultSettingsError
from viperleed.gui.measure.classes.settings import NoSettingsError
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings
from viperleed.gui.measure.classes.settings import interpolate_config_path
from viperleed.gui.measure.hardwarebase import class_from_name
from viperleed.gui.measure.hardwarebase import disconnected_slot
from viperleed.gui.measure.hardwarebase import emit_error
from viperleed.gui.measure.hardwarebase import get_devices
from viperleed.gui.measure.widgets.collapsibleviews import (
    CollapsibleCameraView,
    CollapsibleControllerView,
    CollapsibleDeviceView,
    )
from viperleed.gui.measure.widgets.pathselector import PathSelector
from viperleed.gui.widgets.buttons import QNoDefaultPushButton
from viperleed.gui.widgets.buttons import QUncheckableButtonGroup
from viperleed.gui.widgets.collapsible import CollapsibleList
from viperleed.gui.widgets.collapsible import _PIXEL_SPACING
from viperleed.gui.widgets.lib import remove_spacing_and_margins


class CollapsibleDeviceList(CollapsibleList):
    """A widget containing an arbitrary number of CollapsibleDeviceViews."""

    _top_labels = ('Device', 'Use')
    _device_label = 'device'
    _view_type = CollapsibleDeviceView

    # This signal is emitted when the device/quantity selection changes.
    settings_changed = qtc.pyqtSignal()

    # This signal is emitted when settings, which decide
    # whether a measurement can start at all, change.
    settings_ok_changed = qtc.pyqtSignal()

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, parent=None):
        """Initialise widget.

        Parameters
        ----------
        parent : QObject, optional
            The parent QObject of this widget.

        Returns
        -------
        None.
        """
        super().__init__(parent=parent)
        self.setVerticalScrollBarPolicy(qtc.Qt.ScrollBarAlwaysOn)
        self.setHorizontalScrollBarPolicy(qtc.Qt.ScrollBarAlwaysOff)
        self.default_settings_folder = None
        # If requires_device is True, the CollapsibleDeviceList must
        # contain at least one selected valid device for the
        # are_settings_ok() method to return True.
        self.requires_device = False
        self._add_top_widget_types(qtw.QCheckBox)

    def add_new_view(self, name, cls_and_info):
        """Add a new view.

        Parameters
        ----------
        name : str
            The display name of the device that will be displayed
            on the button of the CollapsibleDeviceView.
        cls_and_info : tuple
            Tuple of the device class and the SettingsInfo required to
            determine the device settings.

        Returns
        -------
        view : CollapsibleDeviceView
            The added view.
        """
        view = self._view_type()
        view.button.setEnabled(False)
        view.button.setText(name)
        view.set_device(*cls_and_info)
        if self.default_settings_folder:
            view.set_settings_folder(self.default_settings_folder)
        self.insert_view(view)
        view.settings_changed.connect(self._emit_and_update_settings)
        return view

    def are_settings_ok(self):
        """Return whether the device selection is acceptable.

        Returns
        -------
        settings_ok : bool
            Whether the settings selected in the CollapsibleDeviceList
            are acceptable or not. If a selected device is either a
            dummy (disconnected) device, or does not have any settings,
            the return value must be False because a measurement with
            such a device will fail. If requires_device is False after
            the first check, then the return value will be True because
            the settings do not have to contain a selected device. If
            requires_device is True, then the return value will only be
            True if at least one device has been selected.
        reason : str
            A descriptive string elaborating why the settings
            are not acceptable.
        """
        for view in self.enabled_views:
            if not view.has_hardware_interface:
                # The selected device is not connected.
                reason = (f'At least one of the selected {self._device_label}s'
                          ' is not connected.')
                return False, reason
            if not view.settings_file:
                # The selected device lacks a suitable settings file.
                reason = (f'At least one of the selected {self._device_label}s'
                          ' does not have a valid settings file.')
                return False, reason
        # Check if any device has been selected, if a selected device
        # is required.
        if not self.requires_device or any(self.enabled_views):
            return True, ''
        return False , f'At least one {self._device_label} must be selected.'

    def event(self, event):
        """Extend event to match QScrollArea width to required width."""
        width = self.widget().sizeHint().width()+10 if self.widget() else 10
        if width > self.minimumWidth():
            self.setMinimumWidth(width)
        return super().event(event)

    def store_settings(self):
        """Store the settings of the selected devices."""
        for view in self.enabled_views:
            view.store_settings()

    def _checkbox(self, view):
        """Return the QCheckBox of a specific view."""
        return self.views[view][0]

    @qtc.pyqtSlot()
    def _detect_and_add_devices(self):
        """Detect devices and add them as views."""
        self.clear()
        self._make_top_items()
        try:
            detected_devices = self._detect_devices()
        except DefaultSettingsError:
            detected_devices = {}
        for name, cls_and_info in detected_devices.items():
            self.add_new_view(name, cls_and_info)

    @emit_default_faulty
    def _detect_devices(self):
        """Detect and return devices."""
        return get_devices(self._device_label)

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(int)
    @qtc.pyqtSlot(bool)
    def _emit_and_update_settings(self, *_):
        """Emit settings changed, check and update settings.

        Emits
        -----
        settings_changed
            To notify other widgets of a settings change.
        settings_ok_changed
            To trigger a settings check.
        """
        self.settings_changed.emit()
        self.settings_ok_changed.emit()
        self._update_stored_settings()

    def _get_relative_path(self, path):
        """Get a str path that is relative to the default path if possible.

        Parameters
        ----------
        path : Path
            The settings path that should be turned into a relative path.

        Returns
        -------
        new_path : str
            The relative path as a str.
        """
        try:
            path = path.relative_to(self.default_settings_folder)
        except (ValueError, AttributeError):
            new_path = str(path.as_posix())
        else:
            new_path = str(path.as_posix())
            new_path = '__CONFIG__/' + new_path
        return new_path

    def _make_top_items(self):
        """Make top labels and add them to the QScrollArea.

        Create the top refresh button, force the QScrollArea to
        take a certain size, and set the top labels.

        Returns
        -------
        None.
        """
        button = QNoDefaultPushButton()
        button.setText(f'Refresh {self._device_label}s')
        button.clicked.connect(self._detect_and_add_devices)
        # Button height*15 ensures that the frame shows at least
        # one fully expanded ViPErLEED primary controller view.
        minimum_height = button.sizeHint().height()*15
        minimum_height *= 1.3
        self.setMinimumHeight(round(minimum_height))
        self._layout.insertWidget(0, button)

        top_labels = qtw.QHBoxLayout()
        remove_spacing_and_margins(top_labels)
        top_labels.addStretch(1)
        device_label = qtw.QLabel(self._top_labels[0])
        top_labels.addWidget(device_label)
        top_labels.addStretch(1)
        for label in self._top_labels[1:]:
            q_label = qtw.QLabel(label)
            top_labels.addWidget(q_label)
            self._widths[label] = q_label.sizeHint().width()
            top_labels.addSpacing(_PIXEL_SPACING)
        self._layout.insertLayout(1, top_labels)

    def _update_stored_settings(self):
        """Update the internally stored settings."""
        # Must be implemented in subclasses if the subclass
        # stores device settings internally.


class CollapsibleCameraList(CollapsibleDeviceList):
    """A CollapsibleList for cameras."""

    _top_labels = ('Cameras', 'Use')
    _device_label = 'camera'
    _view_type = CollapsibleCameraView

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
        self._make_top_items()
        self.setMinimumWidth(self.widget().sizeHint().width() +
                             PathSelector().sizeHint().width())
        self._camera_settings = ()

    def get_camera_settings(self):
        """Return a tuple of camera settings."""
        settings_files = [self._get_relative_path(view.settings_file)
                          for view in self.enabled_views
                          if view.settings_file]
        return tuple(settings_files)

    def set_cameras_from_settings(self, meas_settings):
        """Attempt to select the cameras from settings.

        Parameters
        ----------
        meas_settings : ViPErLEEDSettings
            The settings of the loaded measurement.

        Returns
        -------
        None.
        """
        self._detect_and_add_devices()
        self._camera_settings = meas_settings.getsequence(
            'devices', 'cameras', fallback=()
            )
        self._set_camera_settings()

    def _add_top_widgets_to_view(self, view):
        """Add the top widget types to the CollapsibleView.

        Parameters
        ----------
        view : CollapsibleView
            The CollapsibleView to which the widgets will be attached.

        Returns
        -------
        None.
        """
        super()._add_top_widgets_to_view(view)
        self._checkbox(view).stateChanged.connect(view.set_expanded_state)
        self._checkbox(view).stateChanged.connect(
            self._emit_and_update_settings
            )
        view.set_top_widget_geometry(
            self._checkbox(view), width=self._widths[self._top_labels[1]]
            )

    @qtc.pyqtSlot()
    def _detect_and_add_devices(self):
        """Detect cameras, add them as views and preselect them.

        Emits
        -----
        settings_ok_changed
            At the end of setting devices from
            settings to trigger a check.
        """
        super()._detect_and_add_devices()
        self._set_camera_settings()
        self.settings_ok_changed.emit()

    def _set_camera_settings(self):
        """Set camera settings."""
        for settings in self._camera_settings:
            self._set_single_camera_settings(settings)

    def _set_single_camera_settings(self, camera_settings):
        """Set settings of one camera.

        Parameters
        ----------
        camera_settings : path-like
            The path to the settings of the camera.

        Returns
        -------
        None.

        Emits
        -----
        error_occurred
            If the settings are missing or corrupted.
        """
        try:
            settings = ViPErLEEDSettings.from_settings(camera_settings)
        except NoSettingsError:
            error_settings = [camera_settings,]
            interpolate_config_path(error_settings)
            emit_error(self,
                       QObjectSettingsErrors.SPECIFIED_SETTINGS_CORRUPTED,
                       error_settings[0],)
            return
        device_name = settings.get('camera_settings', 'device_name')
        correct_view = None
        for view in self.views:
            if device_name == view.device_info.unique_name:
                correct_view = view
                break
        else:
            name = device_name + ' (not found)'
            cls = class_from_name(
                'camera', settings.get('camera_settings', 'class_name')
                )
            info = {}
            present = False  # Because there is no hardware interface.
            settings_info = SettingsInfo(name, present, info)
            correct_view = self.add_new_view(name, (cls, settings_info))

        with disconnected_slot(self._emit_and_update_settings,
                               self._checkbox(correct_view).stateChanged,
                               type=qtc.Qt.UniqueConnection):
            correct_view.original_settings = settings.last_file
            self._checkbox(correct_view).setChecked(True)

    def _update_stored_settings(self):
        """Update the internally stored camera settings."""
        self._camera_settings = self.get_camera_settings()


class CollapsibleControllerList(CollapsibleDeviceList):
    """A CollapsibleList for controllers."""

    _top_labels = ('Controllers', 'Use', 'Sets energy')
    _device_label = 'controller'
    _view_type = CollapsibleControllerView

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
        self._add_top_widget_types(qtw.QRadioButton)
        self._radio_buttons = QUncheckableButtonGroup()
        self._make_top_items()
        self.setMinimumWidth(self.widget().sizeHint().width() +
                             PathSelector().sizeHint().width())
        self._primary_settings = ()
        self._secondary_settings = ()

    def get_primary_settings(self):
        """Return a tuple of primary controller settings."""
        for view in self.enabled_views:
            if not self._radiobutton(view).isChecked():
                continue
            return self._get_selected_controller_settings(view)
        return ()

    def get_secondary_settings(self):
        """Return a tuple of secondary controller settings."""
        settings_files = []
        for view in self.enabled_views:
            if self._radiobutton(view).isChecked() or not view.settings_file:
                continue
            settings_files.append(self._get_selected_controller_settings(view))
        return tuple(settings_files)

    def set_controllers_from_settings(self, meas_settings):
        """Attempt to select controllers from settings.

        Parameters
        ----------
        meas_settings : ViPErLEEDSettings
            The settings of the loaded measurement.

        Returns
        -------
        None.
        """
        self._detect_and_add_devices()
        self._primary_settings = meas_settings.getsequence(
            'devices', 'primary_controller', fallback=()
            )
        self._secondary_settings = meas_settings.getsequence(
            'devices', 'secondary_controllers', fallback=()
            )
        self._set_primary_from_settings()
        self._set_secondary_from_settings()

    def _add_top_widgets_to_view(self, view):
        """Add the top widget types to the CollapsibleView.

        Parameters
        ----------
        view : CollapsibleView
            The CollapsibleView to which the widgets will be attached.

        Returns
        -------
        None.
        """
        super()._add_top_widgets_to_view(view)
        selected = self._checkbox(view)
        is_primary = self._radiobutton(view)
        selected.stateChanged.connect(view.set_expanded_state)
        selected.stateChanged.connect(self._set_a_primary)
        selected.stateChanged.connect(self._on_unchecked)
        view.set_top_widget_geometry(selected,
                                     width=self._widths[self._top_labels[1]])
        self._radio_buttons.addButton(is_primary)
        is_primary.toggled.connect(view.set_primary)
        view.set_top_widget_geometry(is_primary,
                                     width=self._widths[self._top_labels[2]])
        is_primary.setEnabled(False)
        is_primary.toggled.connect(self._emit_and_update_settings)

    @qtc.pyqtSlot()
    def _detect_and_add_devices(self):
        """Detect controllers, add them as views, and preselect them.

        Emits
        -----
        settings_ok_changed
            At the end of setting devices from
            settings to trigger a check.
        """
        super()._detect_and_add_devices()
        self._set_primary_from_settings()
        self._set_secondary_from_settings()
        self.settings_ok_changed.emit()

    def _get_selected_controller_settings(self, view):
        """Return the settings for the selected controller."""
        if not view.settings_file:
            return tuple()
        rel_path = self._get_relative_path(view.settings_file)
        return (rel_path, view.selected_quantities)

    @qtc.pyqtSlot(int)
    def _on_unchecked(self, state):
        """Emit settings changed, check and update settings if unchecked.

        Parameters
        ----------
        state : int
            The state the checkbox is in now.
            0 (qtc.Qt.Unchecked) means unchecked.

        Emits
        -----
        settings_changed
            To notify other widgets of a settings change.
        settings_ok_changed
            To trigger a settings check.
        """
        if state == qtc.Qt.Unchecked:
            self._emit_and_update_settings()

    def _radiobutton(self, view):
        """Return the QRadioButton of a specific view."""
        return self.views[view][1]

    @qtc.pyqtSlot(int)
    @qtc.pyqtSlot(bool)
    def _set_a_primary(self, enable):
        """Enable/disable QRadioButton and set a controller as primary
           if necessary.

        Parameters
        ----------
        enable : bool or int
            Converted to bool. True means the associated controller
            will be selected as the primary controller if there is no
            other controller that is selected as the primary controller.
            False means another controller will be selected as the
            primary controller if the associated controller was selected
            as the primary controller.

        Returns
        -------
        None.
        """
        check_box = self.sender()
        enable = bool(enable)
        no_primary_selected = self._radio_buttons.checkedId() == -1
        try:
            radio_btn = next(r for c, r in self.views.values()
                             if c is check_box)
        except StopIteration:
            raise RuntimeError(f'Could not find checkbox {check_box}')
        was_primary = radio_btn.isEnabled()
        radio_btn.setEnabled(enable)
        if no_primary_selected and enable:
            radio_btn.setChecked(True)
        elif radio_btn.isChecked() and not enable:
            self._radio_buttons.uncheck_buttons()
        # The following code selects the first available controller that
        # can be a primary controller as the primary controller if the
        # previously selected primary controller has been disabled.
        if was_primary and not enable:
            for _, radio in self.views.values():
                if radio.isEnabled():
                    radio.setChecked(True)
                    break

    def _set_controller_settings(self, controller_settings):
        """Set settings of controller.

        Parameters
        ----------
        controller_settings : tuple
            A tuple composed of the path to the settings of the
            controller and another tuple containing the quantities
            measured by the controller for this measurement.

        Emits
        -----
        error_occurred
            If the settings are missing or corrupted.
        """
        file, quantities = controller_settings
        try:
            settings = ViPErLEEDSettings.from_settings(file)
        except NoSettingsError:
            error_settings = list(controller_settings)
            interpolate_config_path(error_settings)
            emit_error(self,
                       QObjectSettingsErrors.SPECIFIED_SETTINGS_CORRUPTED,
                       error_settings[0],)
            return
        device_name = settings.get('controller', 'device_name')
        correct_view = None
        for view in self.views:
            if device_name == view.device_info.more['name']:
                correct_view = view
                break
        else:
            name = device_name + ' (not found)'
            cls = class_from_name(
                'controller', settings.get('controller', 'controller_class')
                )
            info = {
                'address': NO_HARDWARE_INTERFACE,
                'name': device_name,
                }
            present = False  # Because there is no hardware interface.
            settings_info = SettingsInfo(name, present, info)
            correct_view = self.add_new_view(name, (cls, settings_info))

        with disconnected_slot(self._emit_and_update_settings,
                               self._checkbox(correct_view).stateChanged,
                               self.views[correct_view][1].toggled,
                               type=qtc.Qt.UniqueConnection):
            correct_view.loading_measurement_settings = True
            correct_view.original_settings = settings.last_file
            self._checkbox(correct_view).setChecked(True)
            correct_view.set_quantities(quantities)

    def _set_primary_from_settings(self):
        """Attempt to select the primary controller from settings."""
        if not self._primary_settings:
            return
        self._set_controller_settings(self._primary_settings)

    def _set_secondary_from_settings(self):
        """Attempt to select the secondary controllers from settings."""
        if not self._secondary_settings:
            return
        for settings in self._secondary_settings:
            self._set_controller_settings(settings)

    def _update_stored_settings(self):
        """Update the internally stored controller settings."""
        self._primary_settings = self.get_primary_settings()
        self._secondary_settings = self.get_secondary_settings()
