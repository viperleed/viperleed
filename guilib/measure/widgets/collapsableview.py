"""Module collapsableview of viperleed.guilib.measure.widgets.

===============================================
      ViPErLEED Graphical User Interface
===============================================

Created: 2024-07-05
Author: Michele Riva
Author: Florian Doerr

Defines the CollapsableView, CollapsableList, CollapsableCameraList, and
CollapsableControllerList classes.
"""

from ast import literal_eval
from pathlib import Path

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure.classes.abc import SettingsInfo
from viperleed.guilib.measure.classes.decorators import emit_default_faulty
from viperleed.guilib.measure.classes.settings import DefaultSettingsError
from viperleed.guilib.measure.classes.settings import ViPErLEEDSettings
from viperleed.guilib.measure.controller.abc import NO_HARDWARE_INTERFACE
from viperleed.guilib.measure.dialogs.settingsdialog import (
    SettingsDialogSectionBase
    )
from viperleed.guilib.measure.dialogs.settingsdialog import SettingsTag
from viperleed.guilib.measure.hardwarebase import class_from_name
from viperleed.guilib.measure.hardwarebase import get_devices
from viperleed.guilib.measure.hardwarebase import make_device
from viperleed.guilib.measure.hardwarebase import safe_connect
from viperleed.guilib.measure.hardwarebase import safe_disconnect
from viperleed.guilib.measure.widgets.pathselector import PathSelector
from viperleed.guilib.widgets.basewidgets import QNoDefaultPushButton
from viperleed.guilib.widgets.basewidgets import QNoDefaultIconButton
from viperleed.guilib.widgets.basewidgets import QUncheckableButtonGroup


_ALIGN_CTR = qtc.Qt.AlignHCenter
_PIXEL_SPACING = 4


def remove_spacing_and_margins(layout):
    """Remove spacing and margins from a layout."""
    layout.setSpacing(0)
    layout.setContentsMargins(0, 0, 0, 0)


class QuantitySelector(qtw.QFrame):
    """A widget that allows selection of quantities from settings."""

    # This is signal is emitted when the quantity selection changes.
    settings_changed = qtc.pyqtSignal()

    def __init__(self, settings, parent=None):
        """Initialise widget."""
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
        """Return the selected quantties."""
        quantities = []
        for group in self._selections:
            if group.checkedButton():
                quantities.append(group.checkedButton().text())
        return tuple(quantities)

    def set_quantities(self, quantities):
        """Set quantities from settings."""
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
                raise RuntimeError('Passed wrong quantities.')


class CollapsableView(qtw.QWidget):
    """A widget that can be expanded/collapsed by the press of a button."""

    def __init__(self, parent=None):
        """Initialise widget."""
        super().__init__(parent=parent)
        self._button = QNoDefaultIconButton()
        self._outer_layout = qtw.QHBoxLayout()
        self._frame = qtw.QFrame(parent=self)
        self._compose_and_connect()

    @property
    def button(self):
        """Return the main button."""
        return self._button

    def _adjust_button_icon(self, button_up):
        """Change in which direction the icon is pointing."""
        if button_up:
            icon = qtw.QStyle.SP_TitleBarShadeButton
        else:
            icon = qtw.QStyle.SP_TitleBarUnshadeButton
        self.button.setIcon(self.style().standardIcon(icon))

    def _compose_and_connect(self):
        """Compose and connect."""
        layout = qtw.QVBoxLayout()
        remove_spacing_and_margins(layout)
        layout.addWidget(self.button)
        policy = self.button.sizePolicy()
        self.button.setIcon(
            self.style().standardIcon(qtw.QStyle.SP_TitleBarUnshadeButton)
            )

        inner_layout = qtw.QVBoxLayout()
        remove_spacing_and_margins(inner_layout)
        self._frame.setFrameStyle(self._frame.StyledPanel | self._frame.Plain)
        self._frame.setLayout(inner_layout)
        self._frame.setVisible(False)

        frame_layout = qtw.QHBoxLayout()
        frame_layout.addSpacing(1)
        frame_layout.addWidget(self._frame)
        frame_layout.addSpacing(1)
        layout.addLayout(frame_layout)
        policy.setHorizontalPolicy(policy.Expanding)
        self.button.setSizePolicy(policy)

        layout.addStretch(1)
        remove_spacing_and_margins(self._outer_layout)
        self._outer_layout.addLayout(layout)
        self.setLayout(self._outer_layout)

        self.button.clicked.connect(self._change_frame_visibility)

    @qtc.pyqtSlot()
    def _change_frame_visibility(self):
        """Switch frame visibility on and off."""
        self._frame.setVisible(not self._frame.isVisible())
        self._adjust_button_icon(self._frame.isVisible())

    def add_collapsable_item(self, item):
        """Add widget to the widgets in the inner collapsable layout."""
        if isinstance(item, qtw.QWidget):
            self._frame.layout().addWidget(item)
        elif isinstance(item, qtw.QLayout):
            self._frame.layout().addLayout(item)
        else:
            raise RuntimeError('Cannot add items to the collapsable view '
                               'that are neither a layout nor a widget.')

    def add_top_widget(self, widget, width=None, align=_ALIGN_CTR):
        """Add widget to the widgets in the outer top layout."""
        layout = qtw.QVBoxLayout()
        height = int((self.button.sizeHint().height() -
                      widget.sizeHint().height())/2)
        if height > 0:
            layout.addSpacing(height)
        surrounding_widget = qtw.QWidget(parent=self)
        widget_layout = qtw.QHBoxLayout()
        remove_spacing_and_margins(widget_layout)
        widget_layout.addWidget(widget)
        layout.addWidget(surrounding_widget)
        surrounding_widget.setLayout(widget_layout)
        layout.addStretch(1)
        self._outer_layout.addLayout(layout)
        self._outer_layout.addSpacing(_PIXEL_SPACING)
        self.set_top_widget_geometry(widget, width=width, align=align)

    @qtc.pyqtSlot(int)
    @qtc.pyqtSlot(bool)
    def enable_view(self, enable):
        """Collapse inner layout."""
        enable = bool(enable)
        self.button.setEnabled(enable)
        self._frame.setVisible(enable)
        self._adjust_button_icon(enable)

    def set_top_widget_geometry(self, widget, width=None, align=_ALIGN_CTR):
        """Set top widget geometry to the given parameters."""
        for lay_i in range(1, self._outer_layout.count()):
            layout = self._outer_layout.itemAt(lay_i)
            try:
                lay_widget = layout.itemAt(1).widget()
            except AttributeError:
                # This means the item in question is a QSpacerItem.
                continue
            inner_layout = lay_widget.layout()
            inner_widget = inner_layout.itemAt(0).widget()
            if inner_widget is not widget:
                continue
            policy = lay_widget.sizePolicy()
            inner_layout.setAlignment(widget, align)
            if width:
                lay_widget.setMinimumWidth(width)
                lay_widget.setMaximumWidth(width)
            policy.setHorizontalPolicy(policy.Fixed)
            lay_widget.setSizePolicy(policy)
            return


class CollapsableDeviceView(CollapsableView):
    """A CollapsableView for measurement hardware."""

    # This is signal is emitted when device settings stored
    # in the measurement settings file change.
    # (Measured quantities and the device settings file.)
    settings_changed = qtc.pyqtSignal()

    def __init__(self, parent=None):
        """Initialise widget."""
        self._settings_folder = PathSelector(select_file=False)
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
        """Set the original settings."""
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
        """Get the handler widgets for the device settings."""
        new_widget = qtw.QWidget()
        layout = qtw.QVBoxLayout()
        for widget in self._handler.get_widgets_with_tags(
                        SettingsTag.MEASUREMENT):
            if isinstance(widget, SettingsDialogSectionBase):
                layout.addWidget(widget)
                for option in widget.options:
                    option.setVisible(True)
                continue
            else:
                widget.setVisible(True)
                _form = qtw.QFormLayout()
                _form.addRow(*widget)
                layout.addLayout(_form)
        new_widget.setLayout(layout)
        self.add_collapsable_item(new_widget)

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(int)
    def _check_for_handler(self):
        """Check if creating a handler is possible."""
        if not self._enabled or not self._settings_file_selector.isEnabled():
            if not self.is_dummy_device():
                return
        self._get_settings_handler()

    def _check_if_settings_changed(self):
        """Check if the current settings differ from the original settings."""
        if self.settings_file != self._original_settings:
            self._original_settings = self.settings_file
            self.settings_changed.emit()

    def _compose_and_connect(self):
        """Compose and connect."""
        super()._compose_and_connect()
        self._settings_folder.path = Path().resolve()
        self.add_collapsable_item(self._settings_folder)
        self.add_collapsable_item(self._settings_file_selector)
        self._settings_folder.path_changed.connect(
            self._get_device_settings_files
            )
        self._settings_file_selector.setEnabled(False)
        self._settings_file_selector.currentIndexChanged.connect(
            self._check_for_handler
            )

    def _make_device(self, **kwargs):
        """Get an instance of the handled device."""
        if not self.settings_file:
            return
        if not self._device_cls:
            return
        if not self._device_info:
            return
        device = make_device(self.settings_file, self._device_cls,
                             self._device_info, **kwargs)
        return device

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
        return not self._device_info.hardware_interface

    def _make_handler_for_device(self, device):
        """Make a SettingsHandler for the device."""
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
        """Collapse inner layout and remove handler when disabled."""
        self._enabled = bool(enable)
        if not self._enabled:
            self._remove_handler()
        else:
            self._check_for_handler()
        super().enable_view(enable)

    def set_device(self, device_cls, device_info):
        """Set the device for which this view is meant."""
        self._device_cls = device_cls
        self._device_info = device_info

    def set_settings_folder(self, settings_folder_path):
        """Set the path of the settings folder to the given path."""
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


class CollapsableCameraView(CollapsableDeviceView):
    """A CollapsableView for cameras."""

    def _connect_camera(self, device):
        """Connect camera."""
        if self.is_dummy_device():
            return
        device.connect_()
        device.start()

    def _get_settings_handler(self):
        """Get the settings handler of the handled device."""
        device = self._make_device()
        # TODO: handle cameras from settings that are not connected
        if not device:
            return
        try:
            self._connect_camera(device)
        except device.exceptions:
            pass
        self._make_handler_for_device(device)
        self._build_device_settings()
        device.disconnect_()


class CollapsableControllerView(CollapsableDeviceView):
    """A CollapsableView for controllers."""

    def __init__(self, parent=None):
        """Initialise widget."""
        super().__init__(parent=parent)
        self._is_primary = False
        self._trying_to_get_settings = False
        self._primary_changed = False
        self._quantity_selector = None
        self._quantities_to_set = ()

    def set_quantities(self, quantities):
        """Set the quantities to measure."""
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
        """Return the selected quantties."""
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
        """Check if the current settings differ from the original settings."""
        if self._quantities_to_set != self.selected_quantities:
            self._quantities_to_set = self.selected_quantities
            self.settings_changed.emit()

    def _get_settings_handler(self):
        """Get the settings handler of the handled device."""
        device = self._make_device(sets_energy=self._is_primary)
        if not device:
            return
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
        """Set whether the device is the primary controller."""
        was_primary = self._is_primary
        self._is_primary = is_primary
        if was_primary == is_primary:
            return
        if self._trying_to_get_settings:
            self._primary_changed = True
        else:
            self._check_for_handler()


class CollapsableDeviceList(qtw.QScrollArea):
    """A widget composed of an arbitrary number of CollapsableViews."""

    _top_labels = ('Device', )

    # This is signal is emitted when the device/quantity selection changes.
    settings_changed = qtc.pyqtSignal()

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, parent=None):
        """Initialise widget."""
        super().__init__(parent=parent)
        self._views = {}
        self._widths = {}
        self._top_widget_types = []
        self._device_type = 'device'
        self._layout = qtw.QVBoxLayout()
        self._layout.setSpacing(0)
        self._default_settings_folder = None
        self.setWidgetResizable(True)
        self.setVerticalScrollBarPolicy(qtc.Qt.ScrollBarAlwaysOn)
        self._make_scroll_area()

    @property
    def default_settings_folder(self):
        """Return the default settings folder."""
        return self._default_settings_folder

    @default_settings_folder.setter
    def default_settings_folder(self, settings_folder_path):
        """Set the default settings folder."""
        self._default_settings_folder = settings_folder_path

    @property
    def views(self):
        """Return views."""
        return self._views

    def _add_top_widget_types(self, *widg_types):
        """Extend the list of widgets to add."""
        self._top_widget_types.extend(widg_types)

    def _add_top_widgets_to_view(self, view):
        """Add the top widget types to the CollapsableView."""
        self.views[view] = []
        for widget_type in self._top_widget_types:
            widget = widget_type()
            view.add_top_widget(widget)
            self.views[view].append(widget)

    def _update_stored_settings(self):
        """Update the internally stored settings."""
        # Must be implemented in subclasses if the subclass
        # stores device settings internally.

    @qtc.pyqtSlot()
    def _detect_and_add_devices(self):
        """Detect devices and add them as views."""
        self._views = {}
        self._layout = qtw.QVBoxLayout()
        self._layout.setSpacing(0)
        self._make_scroll_area()
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
        return get_devices(self._device_type)

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(int)
    @qtc.pyqtSlot(bool)
    def _emit_settings_changed(self, *_):
        """Emit if any settings changes."""
        self.settings_changed.emit()
        self._update_stored_settings()

    def _get_relative_path(self, path):
        """Get a str path that is relative to the default path if possible."""
        try:
            new_path = str(path.relative_to(self.default_settings_folder))
        except ValueError:
            new_path = str(path)
        else:
            new_path = '__CONFIG__/' + new_path
        new_path.replace('\\', '/')
        return new_path

    def _make_scroll_area(self):
        """Compose QScrollArea."""
        widget = qtw.QWidget(parent=self)
        self._layout.addStretch(1)
        widget.setLayout(self._layout)
        self.setWidget(widget)

    def _make_top_items(self):
        """Make top labels and add them to the QScrollArea."""
        button = QNoDefaultPushButton()
        button.setText('Refresh ' + self._device_type + 's')
        button.clicked.connect(self._detect_and_add_devices)
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

    def add_new_view(self, view, name, cls_and_info):
        """Add a new view."""                                                   # TODO: add doc string
        view.button.setEnabled(False)
        view.button.setText(name)
        view.set_device(*cls_and_info)
        if self.default_settings_folder:
            view.set_settings_folder(self.default_settings_folder)
        self._add_top_widgets_to_view(view)
        self._layout.insertWidget(self._layout.count()-1, view)
        view.settings_changed.connect(self._emit_settings_changed)
        return view

    def store_settings(self):
        """Store the settings of the selected devices."""
        try:
            index = self._top_labels.index('Use') - 1
        except ValueError:
            return
        for view, widgets in self.views.items():
            if widgets[index].isChecked():
                view.store_settings()


class CollapsableCameraList(CollapsableDeviceList):
    """A CollapsableList for cameras."""

    _top_labels = ('Cameras', 'Use',)

    def __init__(self, parent=None):
        """Initialise widget."""
        super().__init__(parent=parent)
        self._device_type = 'camera'
        self._add_top_widget_types(qtw.QCheckBox)
        self._make_top_items()
        self.setMinimumWidth(self.widget().sizeHint().width() +
                             PathSelector().sizeHint().width())
        self._camera_settings = ()

    def _add_top_widgets_to_view(self, view):
        """Add the top widget types to the CollapsableView."""
        super()._add_top_widgets_to_view(view)
        self._views[view][0].stateChanged.connect(view.enable_view)
        self._views[view][0].stateChanged.connect(self._emit_settings_changed)
        view.set_top_widget_geometry(
            self._views[view][0], width=self._widths[self._top_labels[1]]
            )

    @qtc.pyqtSlot()
    def _detect_and_add_devices(self):
        """Detect controllers, add them as views and preselect them."""
        super()._detect_and_add_devices()
        self._set_camera_settings()

    def _set_camera_settings(self):
        """Set camera settings."""
        for settings in self._camera_settings:
            self._set_single_camera_settings(settings)

    def _set_single_camera_settings(self, camera_settings):
        """Set settings of camera."""
        settings = ViPErLEEDSettings.from_settings(camera_settings)
        device_name = settings.get('camera_settings', 'device_name')
        correct_view = None
        for view in self.views:
            if device_name == view.device_info.unique_name:
                correct_view = view
                break

        if not correct_view:
            name = device_name + ' (not found)'
            cls = class_from_name(
                'camera', settings.get('camera_settings', 'class_name')
                )
            info = {}
            present = False # Because there is no hardware interface.
            settings_info = SettingsInfo(name, present, info)
            correct_view = self.add_new_view(name, (cls, settings_info))

        safe_disconnect(self._views[correct_view][0].stateChanged,
                        self._emit_settings_changed)
        correct_view.original_settings = settings.last_file
        self._views[correct_view][0].setChecked(True)
        safe_connect(self._views[correct_view][0].stateChanged,
                     self._emit_settings_changed,
                     type=qtc.Qt.UniqueConnection)

    def add_new_view(self, name, cls_and_info):
        """Add a new CollapsableCameraView."""                                  # TODO: doc
        view = CollapsableCameraView()
        return super().add_new_view(view, name, cls_and_info)

    def get_camera_settings(self):
        """Return a tuple of camera settings."""
        settings_files = []
        for view, widgets in self.views.items():
            if not widgets[0].isChecked():
                continue
            settings_files.append(self._get_relative_path(view.settings_file))
        return tuple(settings_files)

    def set_cameras_from_settings(self, meas_settings):
        """Attempt to select the cameras from settings."""
        self._detect_and_add_devices()
        self._camera_settings = meas_settings.getsequence(
            'devices', 'cameras', fallback=()
            )
        self._set_camera_settings()


class CollapsableControllerList(CollapsableDeviceList):
    """A CollapsableList for controllers."""

    _top_labels = ('Controllers', 'Use', 'Primary',)

    def __init__(self, parent=None):
        """Initialise widget."""
        super().__init__(parent=parent)
        self._device_type = 'controller'
        self._add_top_widget_types(qtw.QCheckBox, qtw.QRadioButton)
        self._radio_buttons = QUncheckableButtonGroup()
        self._make_top_items()
        self.setMinimumWidth(self.widget().sizeHint().width() +
                             PathSelector().sizeHint().width())
        self._primary_settings = ()
        self._secondary_settings = ()

    def _add_top_widgets_to_view(self, view):
        """Add the top widget types to the CollapsableView."""
        super()._add_top_widgets_to_view(view)
        self._views[view][0].stateChanged.connect(view.enable_view)
        self._views[view][0].stateChanged.connect(self._enable_primary)
        self._views[view][0].stateChanged.connect(self._emit_settings_changed)
        view.set_top_widget_geometry(
            self._views[view][0], width=self._widths[self._top_labels[1]]
            )
        self._radio_buttons.addButton(self._views[view][1])
        self._views[view][1].toggled.connect(view.set_primary)
        view.set_top_widget_geometry(
            self._views[view][1], width=self._widths[self._top_labels[2]]
            )
        self._views[view][1].setEnabled(False)
        self._views[view][1].toggled.connect(self._emit_settings_changed)

    def _update_stored_settings(self):
        """Update the interally stored controller settings."""
        self._primary_settings = self.get_primary_settings()
        self._secondary_settings = self.get_secondary_settings()

    @qtc.pyqtSlot()
    def _detect_and_add_devices(self):
        """Detect controllers, add them as views and preselect them."""
        super()._detect_and_add_devices()
        self._set_primary_from_settings()
        self._set_secondary_from_settings()

    @qtc.pyqtSlot(int)
    @qtc.pyqtSlot(bool)
    def _enable_primary(self, enable):
        """Enable/disable QRadioButton."""
        check_box = self.sender()
        enable = bool(enable)
        for check, radio in self.views.values():
            if check == check_box:
                was_primary = radio.isEnabled()
                radio.setEnabled(enable)
                if self._radio_buttons.checkedId() == -1 and enable:
                    radio.setChecked(True)
                elif radio.isChecked() and not enable:
                    self._radio_buttons.uncheck_buttons()
                break
        if was_primary and not enable:
            for _, radio in self.views.values():
                if radio.isEnabled():
                    radio.setChecked(True)
                    break

    def _set_controller_settings(self, controller_settings):
        """Set settings of controller."""
        file, quantities = controller_settings
        settings = ViPErLEEDSettings.from_settings(file)
        device_name = settings.get('controller', 'device_name')
        correct_view = None
        for view in self.views:
            if device_name == view.device_info.more['name']:
                correct_view = view
                break

        if not correct_view:
            name = device_name + ' (not found)'
            cls = class_from_name(
                'controller', settings.get('controller', 'controller_class')
                )
            info = {}
            info['address'] = NO_HARDWARE_INTERFACE
            info['name'] = device_name
            present = False # Because there is no hardware interface.
            settings_info = SettingsInfo(name, present, info)
            correct_view = self.add_new_view(name, (cls, settings_info))

        safe_disconnect(self._views[correct_view][0].stateChanged,
                        self._emit_settings_changed)
        safe_disconnect(self._views[correct_view][1].toggled,
                        self._emit_settings_changed)
        correct_view.original_settings = settings.last_file
        self._views[correct_view][0].setChecked(True)
        correct_view.set_quantities(quantities)
        safe_connect(self._views[correct_view][0].stateChanged,
                     self._emit_settings_changed,
                     type=qtc.Qt.UniqueConnection)
        safe_connect(self._views[correct_view][1].toggled,
                     self._emit_settings_changed,
                     type=qtc.Qt.UniqueConnection)

    def _set_primary_from_settings(self):
        """Attempt to select the primary controller from settings."""
        if not self._primary_settings:
            return
        self._set_controller_settings(self._primary_settings)

    def _set_secondary_from_settings(self):
        """Attempt to select the secondary controllers from settings."""
        if not self._secondary_settings:
            return
        for controller_settings in self._secondary_settings:
            self._set_controller_settings(controller_settings)

    def add_new_view(self, name, cls_and_info):
        """Add a new CollapsableControllerView."""                              # TODO: doc
        view = CollapsableControllerView()
        return super().add_new_view(view, name, cls_and_info)

    def get_primary_settings(self):
        """Return a tuple of camera settings."""
        for view, widgets in self.views.items():
            if not widgets[0].isChecked() or not widgets[1].isChecked():
                continue
            rel_path = self._get_relative_path(view.settings_file)
            return (rel_path, view.selected_quantities)
        return ()

    def get_secondary_settings(self):
        """Return a tuple of camera settings."""
        settings_files = []
        for view, widgets in self.views.items():
            if not widgets[0].isChecked() or widgets[1].isChecked():
                continue
            settings_files.append((self._get_relative_path(view.settings_file),
                                   view.selected_quantities))
        return tuple(settings_files)

    def set_controllers_from_settings(self, meas_settings):
        """Attempt to select controllers from settings."""
        self._detect_and_add_devices()
        self._primary_settings = meas_settings.getsequence(
            'devices', 'primary_controller', fallback=()
            )
        self._secondary_settings = meas_settings.getsequence(
            'devices', 'secondary_controllers', fallback=()
            )
        self._set_primary_from_settings()
        self._set_secondary_from_settings()
