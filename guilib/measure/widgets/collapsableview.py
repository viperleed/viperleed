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

from viperleed.guilib.measure.dialogs.settingsdialog import SettingsDialog
from viperleed.guilib.measure.widgets.pathselector import PathSelector
from viperleed.guilib.widgets.basewidgets import QNoDefaultPushButton
from viperleed.guilib.measure.hardwarebase import get_devices
from viperleed.guilib.measure.hardwarebase import make_device


_ALIGN_CTR = qtc.Qt.AlignHCenter
_PIXEL_SPACING = 4


def remove_spacing_and_margins(layout):
    """Remove spacing and margins from a layout."""
    layout.setSpacing(0)
    layout.setContentsMargins(0, 0, 0, 0)


class QuantitySelector(qtw.QWidget):                                            # TODO: make title and borders (QFrame?)
    """A widget that allows selection of quantities from settings."""

    # This is signal is emitted when the quantity selection changes.
    settings_changed = qtc.pyqtSignal()

    def __init__(self, settings, parent=None):
        """Initialise widget."""
        super().__init__(parent=parent)
        self._selections = []
        quantities = settings.get('controller', 'measurement_devices',
                                  fallback=())
        self._compose(literal_eval(quantities))

    def _compose(self, all_quantities):
        """Compose quantity widgets."""
        main_layout = qtw.QHBoxLayout()
        for selections in all_quantities:
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
            main_layout.addLayout(layout)
        self.setLayout(main_layout)

    def get_selected_quantities(self):
        """Return the selected quantties."""
        quantities = []
        for group in self._selections:
            if group.checkedButton():
                quantities.append(group.checkedButton().text())
        return tuple(quantities)


class CollapsableView(qtw.QWidget):
    """A widget that can be expanded/collapsed by the press of a button."""

    def __init__(self, parent=None):
        """Initialise widget."""
        super().__init__(parent=parent)
        self._button = QNoDefaultPushButton()
        self._outer_layout = qtw.QHBoxLayout()
        self._frame = qtw.QFrame(parent=self)
        self._compose_and_connect()

    @property
    def button(self):
        """Return the main button."""
        return self._button

    def _compose_and_connect(self):
        """Compose and connect."""
        layout = qtw.QVBoxLayout()
        remove_spacing_and_margins(layout)
        layout.addWidget(self._button)
        policy = self._button.sizePolicy()

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
        self._button.setSizePolicy(policy)

        layout.addStretch(1)
        remove_spacing_and_margins(self._outer_layout)
        self._outer_layout.addLayout(layout)
        self.setLayout(self._outer_layout)

        self._button.clicked.connect(self._change_frame_visibility)

    @qtc.pyqtSlot()
    def _change_frame_visibility(self):
        """Switch frame visibility on and off."""
        self._frame.setVisible(not self._frame.isVisible())

    def add_collapsable_widget(self, widget):
        """Add widget to the widgets in the inner collapsable layout."""
        self._frame.layout().addWidget(widget)

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

    def __init__(self, parent=None):
        """Initialise widget."""
        self._settings_folder = PathSelector(select_file=False)
        self._settings_folder.path = Path().resolve()
        self._settings_file = qtw.QComboBox()
        self._settings_file.setEnabled(False)
        self._get_settings_button = QNoDefaultPushButton()
        self._get_settings_button.setEnabled(False)
        super().__init__(parent=parent)
        self._device_cls = None
        self._device_info = None

    def _compose_and_connect(self):
        """Compose and connect."""
        super()._compose_and_connect()
        self.add_collapsable_widget(self._settings_folder)
        self.add_collapsable_widget(self._settings_file)
        self._settings_folder.path_changed.connect(
            self._get_device_settings_files
            )

        self._get_settings_button.setText('Get settings')
        self._get_settings_button.clicked.connect(self._get_settings_handler)
        self.add_collapsable_widget(self._get_settings_button)

    def set_device(self, device_cls, device_info):
        """Set the device for which this view is meant."""
        self._device_cls = device_cls
        self._device_info = device_info

    def set_settings_folder(self, settings_folder_path):
        """Set the path of the settings folder to the given path."""
        self._settings_folder.path = settings_folder_path

    @qtc.pyqtSlot()
    def _get_device_settings_files(self):
        """Get settings files of the device handled by this view."""
        self._settings_file.clear()
        self._settings_file.setEnabled(False)
        self._get_settings_button.setEnabled(False)
        if not self._device_cls or self._settings_folder.path == Path():
            return

        matching_settings = self._device_cls.find_matching_settings_files(
            self._device_info, self._settings_folder.path, False, False)

        for settings in matching_settings:
            self._settings_file.addItem(settings.stem, userData=settings)
        if any(matching_settings):
            self._settings_file.setEnabled(True)
            self._get_settings_button.setEnabled(True)

    def _make_dialog_for_device(self, device):
        """Make and return a SettingsDialog for the device."""
        device.uses_default_settings = False
        dialog = SettingsDialog(device, parent=self)
        dialog.measurement_settings_only = True
        dialog.adv_button.setVisible(False)
        dialog.finished.connect(device.disconnect_)
        return dialog


class CollapsableCameraView(CollapsableDeviceView):
    """A CollapsableView for cameras."""

    def _get_settings_handler(self):
        """Get the settings handler of the handled device."""
        device = make_device(self._settings_file.currentData(),
                             self._device_cls, self._device_info)
        if not device:
            return
        device.connect_()
        device.start()
        dialog = self._make_dialog_for_device(device)
        dialog.open()


class CollapsableControllerView(CollapsableDeviceView):
    """A CollapsableView for controllers."""

    def __init__(self, parent=None):
        """Initialise widget."""
        super().__init__(parent=parent)
        self._is_primary = False

    def _get_settings_handler(self):
        """Get the settings handler of the handled device."""
        device = make_device(self._settings_file.currentData(),
                             self._device_cls, self._device_info,
                             sets_energy=self._is_primary)
        if not device:
            return
        dialog = self._make_dialog_for_device(device)
        device.ready_to_show_settings.connect(dialog.open)
        device.prepare_to_show_settings()

    @qtc.pyqtSlot(bool)
    def set_primary(self, is_primary):
        """Set whether the device is the primary controller."""
        self._is_primary = is_primary


class CollapsableDeviceList(qtw.QScrollArea):
    """A widget composed of an arbitrary number of CollapsableViews."""

    _top_labels = ('Device', )

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
        self._views[view] = []
        for widget_type in self._top_widget_types:
            widget = widget_type()
            view.add_top_widget(widget)
            self._views[view].append(widget)

    @qtc.pyqtSlot()
    def _detect_devices(self):
        """Detect devices and add them as views."""
        self._layout = qtw.QVBoxLayout()
        self._layout.setSpacing(0)
        self._make_scroll_area()
        self._make_top_items()
        detected_devices = get_devices(self._device_type)
        for name, cls_and_info in detected_devices.items():
            self.add_new_view(name, cls_and_info)

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
        button.clicked.connect(self._detect_devices)
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
        """Add a new view."""
        view.button.setEnabled(False)
        view.button.setText(name)
        view.set_device(*cls_and_info)
        if self.default_settings_folder:
            view.set_settings_folder(self.default_settings_folder)
        self._add_top_widgets_to_view(view)
        self._layout.insertWidget(self._layout.count()-1, view)


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

    def _add_top_widgets_to_view(self, view):
        """Add the top widget types to the CollapsableView."""
        super()._add_top_widgets_to_view(view)
        self._views[view][0].stateChanged.connect(view.enable_view)
        view.set_top_widget_geometry(
            self._views[view][0], width=self._widths[self._top_labels[1]]
            )

    def add_new_view(self, name, cls_and_info):
        """Add a new CollapsableCameraView."""
        view = CollapsableCameraView()
        super().add_new_view(view, name, cls_and_info)


class CollapsableControllerList(CollapsableDeviceList):
    """A CollapsableList for controllers."""

    _top_labels = ('Controllers', 'Use', 'Primary',)

    def __init__(self, parent=None):
        """Initialise widget."""
        super().__init__(parent=parent)
        self._device_type = 'controller'
        self._add_top_widget_types(qtw.QCheckBox, qtw.QRadioButton)
        self._radio_buttons = qtw.QButtonGroup()
        self._make_top_items()
        self.setMinimumWidth(self.widget().sizeHint().width() +
                             PathSelector().sizeHint().width())

    def _add_top_widgets_to_view(self, view):
        """Add the top widget types to the CollapsableView."""
        super()._add_top_widgets_to_view(view)
        self._views[view][0].stateChanged.connect(view.enable_view)
        view.set_top_widget_geometry(
            self._views[view][0], width=self._widths[self._top_labels[1]]
            )
        self._radio_buttons.addButton(self._views[view][1])
        self._views[view][1].toggled.connect(view.set_primary)
        view.set_top_widget_geometry(
            self._views[view][1], width=self._widths[self._top_labels[2]]
            )

    def add_new_view(self, name, cls_and_info):
        """Add a new CollapsableControllerView."""
        view = CollapsableControllerView()
        super().add_new_view(view, name, cls_and_info)
