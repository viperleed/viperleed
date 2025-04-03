"""Private module _meassettings of viperleed.guilib.measure.measurement.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-04-18
Author: Michele Riva
Author: Florian Dörr

This module defines widgets used in the settings dialog for
MeasurementABC and subclass instances.
"""

from ast import literal_eval

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.dialogs.settingsdialog import (
    SettingsDialogSectionBase,
    SettingsTag,
    )
from viperleed.guilib.measure.widgets.collapsiblelists import (
    CollapsibleCameraList, CollapsibleControllerList,
    )
from viperleed.guilib.measure.widgets.fieldinfo import FieldInfo
from viperleed.guilib.measure.widgets.spinboxes import CoercingDoubleSpinBox
from viperleed.guilib.measure.widgets.spinboxes import CoercingSpinBox
from viperleed.guilib.widgets.basewidgets import ButtonWithLabel
from viperleed.guilib.widgets.basewidgets import QNoDefaultDialogButtonBox


DELTA_ENERGY_NAME = 'Step height'
MAX_NUM_STEPS = 7
MAX_DELAY = 65535


class DeviceEditor(SettingsDialogSectionBase):
    """Class for selecting devices and editing their settings."""

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, settings, may_have_cameras=False, default_folder=None,
                 **kwargs):
        """Initialize instance.

        Parameters
        ----------
        settings : ViPErLEEDSettings
            The settings of the loaded measurement.
        may_have_cameras : bool
            Whether a collapsible list for cameras should be added.
        default_folder : path-like or None
            Default folder to look for settings in.
            None if none was given.
        **kwargs : dict
            Keyword arguments passed on to super().__init__
            'display_name' : Displayed name of section.
            'tags' : Tags associated with the section.
            'tooltip' : Tooltip displayed with the section.

        Returns
        -------
        None.
        """
        self._settings = settings
        kwargs.setdefault('display_name', 'Device configuration')
        kwargs.setdefault('tags', SettingsTag.REGULAR)
        kwargs.setdefault('tooltip', 'This section lists devices, allows their'
                          'selection, and the editing of their settings.')
        super().__init__(**kwargs)
        self._controllers = CollapsibleControllerList()
        self._cameras = CollapsibleCameraList()
        self._default_settings_folder = None
        self.default_settings_folder = default_folder
        self.settings_changed.connect(self._store_device_settings)
        self._compose_and_connect_collapsible_lists(may_have_cameras)

    @property
    def _device_lists(self):
        """Return the device lists."""
        return(self._controllers, self._cameras)

    @property
    def default_settings_folder(self):
        """Return the default settings folder."""
        return self._default_settings_folder

    @default_settings_folder.setter
    def default_settings_folder(self, settings_path):
        """Set the default settings folder.

        Parameters
        ----------
        settings_path : Path or str
            The path to the folder containing the settings.

        Returns
        -------
        None.
        """
        self._default_settings_folder = settings_path
        self._cameras.default_settings_folder = settings_path
        self._controllers.default_settings_folder = settings_path

    def _compose_and_connect_collapsible_lists(self, may_have_cameras):
        """Compose the collapsible lists for cameras and controllers.

        Parameters
        ----------
        may_have_cameras : bool
            Whether a collapsible list for cameras should be added.

        Returns
        -------
        None.
        """
        central_layout = qtw.QHBoxLayout()
        central_layout.addWidget(self._controllers)
        if may_have_cameras:
            central_layout.addWidget(self._cameras)
        self.central_widget.setLayout(central_layout)
        for collapsible_list in self._device_lists:
            collapsible_list.settings_changed.connect(self.settings_changed)
            collapsible_list.error_occurred.connect(self.error_occurred)
            collapsible_list.settings_ok_changed.connect(self.settings_ok_changed)
        self._controllers.requires_device = True
        meas = self._settings.get('measurement_settings', 'measurement_class')
        must_have_cameras = ('IVVideo',)
        self._cameras.requires_device = meas in must_have_cameras

    @qtc.pyqtSlot()
    def _store_device_settings(self):
        """Collect device settings from the list viewers."""
        self._settings.set('devices', 'primary_controller',
                           str(self._controllers.get_primary_settings()))
        self._settings.set('devices', 'secondary_controllers',
                           str(self._controllers.get_secondary_settings()))
        self._settings.set('devices', 'cameras',
                           str(self._cameras.get_camera_settings()))

    def are_settings_ok(self):
        """Return whether the device selection is acceptable.

        Returns
        -------
        settings_ok : bool
            Whether the settings selected in the CollapsibleDeviceLists
            are acceptable or not.
        reason : str
            A descriptive string elaborating why the settings
            are not acceptable.
        """
        ctrl_ok, reason_ctrl = self._controllers.are_settings_ok()
        cameras_ok, reason_camera = self._cameras.are_settings_ok()
        reasons = reason_ctrl, reason_camera
        reason = ' '.join(r for r in reasons if r)
        return ctrl_ok and cameras_ok, reason

    def store_lower_level_settings(self):
        """Store the settings of the selected devices."""
        for collapsible_list in self._device_lists:
            collapsible_list.store_settings()

    @qtc.pyqtSlot()
    def update_widgets(self):
        """Update collapsible lists."""
        if not self.isVisible():
            # If the widget is not visible, we do not want to create
            # device objects that may occupy the serial as the
            # showEvent will trigger the device population anyway.
            return
        self._controllers.set_controllers_from_settings(self._settings)
        self._cameras.set_cameras_from_settings(self._settings)


class StepProfileViewer(ButtonWithLabel):
    """Viewer of the current step-profile type.

    Shows the current step profile to the user and opens
    the step profile editor with the push of a button.
    """

    settings_changed = qtc.pyqtSignal()

    def __init__(self, **kwargs):
        """Initialize viewer."""
        super().__init__(**kwargs)
        self.notify_ = self.settings_changed
        self.set_button_text('Edit')
        self.profile_editor = StepProfileEditor(parent=self)
        self._connect()

    def _connect(self):
        """Connect (only once) relevant signals and slots."""
        self.button.clicked.connect(self.profile_editor.show)
        self.profile_editor.accepted.connect(self._on_settings_changed)

    @qtc.pyqtSlot()
    def _on_settings_changed(self):
        """Update step profile to selected profile."""
        self._set_label_text(self.profile_editor.profile[0])
        self.settings_changed.emit()

    def _set_label_text(self, value):
        """Set label text from selected profile."""
        if isinstance(value, str):
            self.set_label_text(value.capitalize() + ' profile')
            return
        self.set_label_text('Custom profile')

    def get_(self):
        """Return the value to be stored in the config."""
        return str(self.profile_editor.profile)

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Connect finished signal when shown."""
        base.safe_connect(self.window().finished, self.profile_editor.reject,
                          type=qtc.Qt.UniqueConnection)
        super().showEvent(event)

    def set_(self, value):
        """Set label and load profile into step profile editor."""
        value = literal_eval(value)
        if not value:
            value = ('abrupt',)
        if isinstance(value, str):
            value = (value,)
        self.profile_editor.profile = value
        self._set_label_text(value[0])


class StepProfileEditor(qtw.QDialog):                                           # TODO: visually draw profiles
    """Editor for setting step profiles.

    Provides a selection of step profiles and allows
    editing of the corresponding settings.
    """

    def __init__(self, **kwargs):
        """Initialize editor."""
        super().__init__(**kwargs)
        self.pick_profile = qtw.QComboBox()
        self._profile_editors = {
            name: cls()
            for name, cls in ProfileEditorBase().subclasses.items()
            }
        self._populate_profile_options()
        self.setWindowTitle('Edit energy-step profile')
        self.setWindowFlags(self.windowFlags()
                            & ~qtc.Qt.WindowContextHelpButtonHint)
        self._delay = qtc.QTimer()
        self._delay.setSingleShot(True)
        self._delay.setInterval(3)
        self._delay.timeout.connect(self.adjustSize)
        self._compose_and_connect()

    @property
    def profile(self):
        """Return the currently selected profile."""
        return self.pick_profile.currentData().profile

    @profile.setter
    def profile(self, profile_data):
        """Set the profile from settings."""
        first = profile_data[0]
        name = (first if isinstance(first, str) else 'custom')
        index = self.pick_profile.findText(name.capitalize() + ' profile')
        self.pick_profile.setCurrentIndex(index)
        self.pick_profile.currentData().set_profile(profile_data)

    def _compose_and_connect(self):
        """Compose StepProfileEditor and connect signals."""
        layout = qtw.QVBoxLayout()
        layout.addLayout(self._compose_editor())
        _bbox = QNoDefaultDialogButtonBox
        buttons = _bbox(_bbox.Ok | _bbox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.setLayout(layout)
        self.pick_profile.currentTextChanged.connect(self._on_profile_selected)
        self._profile_editors['custom'].step_count_reduced.connect(
            self._delay.start
            )

    def _compose_editor(self):
        """Compose editor section."""
        layout = qtw.QHBoxLayout()
        profile_with_stretch = qtw.QVBoxLayout()
        profile_with_stretch.addWidget(self.pick_profile)
        profile_with_stretch.addStretch(1)
        layout.addLayout(profile_with_stretch)
        for profile_editor in self._profile_editors.values():
            layout.addWidget(profile_editor)
            profile_editor.hide()
        return layout

    def _on_profile_selected(self):
        """Update the displayed step profile editor."""
        for profile_editor in self._profile_editors.values():
            profile_editor.hide()
        self.pick_profile.currentData().show()
        self._delay.start()

    def _populate_profile_options(self):
        """Add profile options to profile selection."""
        for profile_editor in self._profile_editors.values():
            name = profile_editor.name.capitalize()+' profile'
            self.pick_profile.addItem(name, userData=profile_editor)

    @qtc.pyqtSlot()
    def accept(self):
        """Store selected profile then accept."""
        self.pick_profile.currentData().update_profile()
        super().accept()


class ProfileEditorBase(qtw.QWidget):
    """Base class for step profiles."""

    _subclasses = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls._subclasses[cls.name] = cls

    def __init__(self):
        """Initialise object."""
        super().__init__()
        self.profile = ()

    @property
    def subclasses(self):
        return self._subclasses.copy()


class AbruptStepEditor(ProfileEditorBase):
    """Abrupt step."""

    name = 'abrupt'

    def __init__(self):
        """Initialise object."""
        super().__init__()
        self.profile = ('abrupt', )

    def set_profile(self, *_):
        """Do nothing."""

    def show(self):
        """Do nothing."""

    def update_profile(self):
        """Do nothing."""


class LinearStepEditor(ProfileEditorBase):
    """Editor for selecting settings of a linear energy profile."""

    name = 'linear'

    def __init__(self):
        """Initialise object."""
        super().__init__()
        self._controls = {
            'n_steps' : CoercingSpinBox(soft_range=(0, MAX_NUM_STEPS)),
            'duration' : CoercingSpinBox(soft_range=(0, MAX_DELAY),
                                         suffix=' ms'),
            }
        self._compose()

    def _compose(self):
        """Place children widgets."""
        layout = qtw.QVBoxLayout()
        layout.addLayout(self._compose_step_nr_selection())
        layout.addLayout(self._compose_duration_selection())
        layout.addStretch(1)
        self.setLayout(layout)

    def _compose_step_nr_selection(self):
        """Return a layout of the step-number selection."""
        layout = qtw.QHBoxLayout()
        step_number_label = qtw.QLabel('Nr. of steps:')
        layout.addWidget(step_number_label)
        size = step_number_label.fontMetrics().boundingRect('a').height()
        info = ('<nobr>The number of intermediate steps.</nobr> '
                f'Cannot be higher than {MAX_NUM_STEPS}.')
        layout.addWidget(FieldInfo(info, size=size))
        layout.addWidget(self._controls['n_steps'])
        return layout

    def _compose_duration_selection(self):
        """Return a layout of the duration selection."""
        layout = qtw.QHBoxLayout()
        duration_label = qtw.QLabel()
        duration_label.setText('Step duration:')
        layout.addWidget(duration_label)
        size = duration_label.fontMetrics().boundingRect('a').height()
        info = 'The settle time after each energy.'
        layout.addWidget(FieldInfo(info, size=size))
        layout.addWidget(self._controls['duration'])
        return layout

    def set_profile(self, profile):
        """Set linear profile."""
        self._controls['n_steps'].setValue(profile[1])
        self._controls['duration'].setValue(profile[2])

    def update_profile(self):
        """Set the profile to the selected values."""
        self.profile = ('linear', self._controls['n_steps'].value(),
                        self._controls['duration'].value())
        if self.profile == ('linear', 0, 0):
            self.profile = ('abrupt',)


class FractionalStepEditor(ProfileEditorBase):
    """Editor for the settings of an energy profile with custom steps."""

    name = 'custom'
    step_count_reduced = qtc.pyqtSignal()

    def __init__(self):
        """Initialise object."""
        super().__init__()
        self._controls = {
            'add_step' : qtw.QPushButton(),
            'remove_step' : qtw.QPushButton(),
            }
        # In order to keep track whether a step has properly been removed
        # from the editor, we have to keep track of how many widgets have
        # been deleted.
        self._n_widgets_removed = 0
        self._steps = []
        self._connect()
        self._compose()

    @qtc.pyqtSlot()
    def _add_step(self, fraction=None, duration=None):
        """Add a step to the fractional step profile."""
        layout = qtw.QHBoxLayout()
        fraction_handler = CoercingDoubleSpinBox(decimals=2, soft_range=(0, 2),
                                                 step=0.05)
        duration_handler = CoercingSpinBox(soft_range=(0, MAX_DELAY),
                                           suffix=' ms')
        for value, handler in zip((fraction, duration),
                                  (fraction_handler, duration_handler)):
            if value:
                handler.setValue(value)
            handler.destroyed.connect(self._emit_step_count_reduced)
            layout.addWidget(handler)
        self._steps.append(layout)
        self._update_button_states()
        self.layout().insertLayout(self.layout().count() - 1, layout)

    def _compose(self):
        """Place children widgets."""
        layout = qtw.QVBoxLayout()
        layout.addLayout(self._compose_buttons())
        layout.addLayout(self._compose_labels())
        layout.addStretch(1)
        self.setLayout(layout)

    def _compose_buttons(self):
        """Return a layout of the add and remove buttons."""
        layout = qtw.QHBoxLayout()
        self._controls['add_step'].setText('Add')
        layout.addWidget(self._controls['add_step'])
        self._controls['remove_step'].setText('Remove')
        layout.addWidget(self._controls['remove_step'])
        self._controls['remove_step'].setEnabled(False)
        return layout

    def _compose_labels(self):
        """Return a layout of the labels."""
        layout = qtw.QHBoxLayout()
        fraction_label = qtw.QLabel()
        fraction_label.setText('Step fraction')
        layout.addWidget(fraction_label)
        size = fraction_label.fontMetrics().boundingRect('a').height()
        info = ('<nobr>The energies to set given as a '
                f'fraction</nobr> of "{DELTA_ENERGY_NAME}". '
                f'Number of steps cannot exceed {MAX_NUM_STEPS}.')
        layout.addWidget(FieldInfo(info, size=size))
        duration_label = qtw.QLabel()
        duration_label.setText('Duration')
        layout.addWidget(duration_label)
        info = '<nobr>The settle time after each</nobr> energy given in ms.'
        layout.addWidget(FieldInfo(info, size=size))
        return layout

    def _connect(self):
        """Connect buttons to meathods."""
        self._controls['add_step'].clicked.connect(self._add_step)
        self._controls['remove_step'].clicked.connect(self._remove_step)

    def _emit_step_count_reduced(self):
        """Emit the step_count_reduced signal once both widgets are deleted."""
        self._n_widgets_removed += 1
        if self._n_widgets_removed < 2:
            return
        self.step_count_reduced.emit()
        self._n_widgets_removed = 0

    @qtc.pyqtSlot()
    def _remove_step(self):
        """Remove a step from the fractional step profile."""
        layout = self.layout()
        if layout.count() < 4:
            return
        item = layout.itemAt(layout.count() - 2)
        layout.removeItem(item)
        del(self._steps[-1])
        self._update_button_states()
        for widget_index in range(item.layout().count()):
            item.itemAt(widget_index).widget().deleteLater()

    def _update_button_states(self):
        """Enable/disable add/remove buttons."""
        self._controls['add_step'].setEnabled(len(self._steps) < MAX_NUM_STEPS)
        self._controls['remove_step'].setEnabled(len(self._steps) != 0)

    def set_profile(self, profile):
        """Set fractional profile."""
        while self.layout().count() >= 4:
            self._remove_step()
        self.profile = profile
        for fraction, duration in zip(profile[0::2], profile[1::2]):
            self._add_step(fraction, duration)
        self._update_button_states()

    def update_profile(self):
        """Set the profile to the selected values."""
        tmp = []
        layout = self.layout()
        for index in range(2, layout.count() - 1):
            item = layout.itemAt(index)
            for widget_index in range(item.layout().count()):
                tmp.append(item.itemAt(widget_index).widget().value())
        self.profile = tuple(tmp)
        if not self.profile or not any(self.profile):
            self.profile = ('abrupt',)
