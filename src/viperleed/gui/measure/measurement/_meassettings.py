"""Private module _meassettings of viperleed.gui.measure.measurement.

This module defines widgets used in the settings dialog for
MeasurementABC and subclass instances.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-04-18'
__license__ = 'GPLv3+'

from abc import abstractmethod
from ast import literal_eval

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.abc import QMetaABC
from viperleed.gui.measure.classes.settings import SystemSettings
from viperleed.gui.measure.dialogs.settingsdialog import (
    SettingsDialogSectionBase,
    SettingsTag,
    )
from viperleed.gui.measure.widgets.collapsiblelists import (
    CollapsibleCameraList,
    CollapsibleControllerList,
    )
from viperleed.gui.measure.widgets.fieldinfo import FieldInfo
from viperleed.gui.measure.widgets.spinboxes import CoercingDoubleSpinBox
from viperleed.gui.measure.widgets.spinboxes import CoercingSpinBox
from viperleed.gui.widgets.buttons import ButtonWithLabel
from viperleed.gui.widgets.buttons import QNoDefaultDialogButtonBox
from viperleed.gui.widgets.buttons import QNoDefaultPushButton


DELTA_E_NAME = '\u0394E'
START_E_NAME = 'E start'
END_E_NAME = 'E end'
MAX_NUM_STEPS = 7
MAX_DELAY = 65535
N_COLUMNS = 2
N_HEADER_ROWS = 2


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
        may_have_cameras : bool, optional
            Whether a collapsible list for cameras should be added.
            Default is False.
        default_folder : path-like or None, optional
            Default folder to look for settings in. Default is None.
            If None is given, then the base system configuration path
            will be used.
        **kwargs : object
            Keyword arguments passed on to SettingsDialogSectionBase
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
        settings_path = SystemSettings().paths['configuration']
        if default_folder:
            self.default_settings_folder = default_folder
        else:
            self.default_settings_folder = settings_path
        self.settings_changed.connect(self._store_device_settings)
        self._compose_and_connect_collapsible_lists(may_have_cameras)

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

    @property
    def _device_lists(self):
        """Return the device lists."""
        return (self._controllers, self._cameras)

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
            collapsible_list.settings_ok_changed.connect(
                self.settings_ok_changed
                )
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
        self.profile_editor = EnergyStepProfileDialog(parent=self)
        self._connect()

    def get_(self):
        """Return the value to be stored in the config."""
        return str(self.profile_editor.profile)

    def set_(self, value):
        """Set label and load profile into step profile editor."""
        value = literal_eval(value)
        if not value:
            value = AbruptEnergyStepEditor().profile
        if isinstance(value, str):
            value = (value,)
        self.profile_editor.profile = value
        self._set_label_text(value[0])

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Connect finished signal when shown."""
        base.safe_connect(self.window().finished, self.profile_editor.reject,
                          type=qtc.Qt.UniqueConnection)
        super().showEvent(event)

    def _connect(self):
        """Connect (only once) relevant signals and slots."""
        self.button.clicked.connect(self.profile_editor.show)
        self.profile_editor.accepted.connect(self._on_settings_changed)

    @qtc.pyqtSlot()
    def _on_settings_changed(self):
        """Update step profile to selected profile."""
        self._set_label_text(self.profile_editor.profile_name)
        self.settings_changed.emit()

    def _set_label_text(self, value):
        """Set label text from selected profile."""
        if isinstance(value, str):
            self.set_label_text(value.capitalize() + ' profile')
            return
        self.set_label_text('Custom profile')


class EnergyStepProfileDialog(qtw.QDialog):                                     # TODO: visually draw profiles
    """A dialog for setting step profiles.

    Provides a selection of step profiles and allows
    editing of the corresponding settings.
    """

    def __init__(self, **kwargs):
        """Initialize editor."""
        super().__init__(**kwargs)
        self.pick_profile = qtw.QComboBox()
        self._profile_editors = {
            name: cls()
            for name, cls in EnergyStepProfileShapeEditor().subclasses.items()
            }
        self._profile_description = qtw.QLabel()
        self._populate_profile_options()
        self.setWindowTitle('Edit energy-step profile')
        self.setWindowFlags(self.windowFlags()
                            & ~qtc.Qt.WindowContextHelpButtonHint)
        # The _adjust_size_timer timer is used to trigger a replot in
        # case the step count is reduced and the EnergyStepProfileDialog
        # therefore needs less space to fit all its widgets.
        self._adjust_size_timer = qtc.QTimer()
        self._adjust_size_timer.setSingleShot(True)
        self._adjust_size_timer.setInterval(3)
        self._adjust_size_timer.timeout.connect(self.adjustSize)
        self._compose_and_connect()

    @property
    def profile_name(self):
        """Return name of the currently selected profile."""
        return self.pick_profile.currentData().name

    @property
    def profile(self):
        """Return the currently selected profile."""
        return self.pick_profile.currentData().profile

    @profile.setter
    def profile(self, profile_data):
        """Set the profile from settings."""
        first = profile_data[0]
        name = (first if isinstance(first, str)
                else FractionalEnergyStepEditor.name)
        index = self.pick_profile.findText(name.capitalize() + ' profile')
        self.pick_profile.setCurrentIndex(index)
        self.pick_profile.currentData().set_profile(profile_data)
        self._profile_description.setText(
            self.pick_profile.currentData().description
            )

    @qtc.pyqtSlot()
    def accept(self):
        """Store selected profile then accept."""
        self.pick_profile.currentData().update_profile()
        super().accept()

    def _compose_and_connect(self):
        """Compose EnergyStepProfileDialog and connect signals."""
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
            self._adjust_size_timer.start
            )

    def _compose_editor(self):
        """Compose editor section."""
        layout = qtw.QHBoxLayout()
        profile_with_stretch = qtw.QVBoxLayout()
        profile_with_stretch.addWidget(self.pick_profile)
        profile_with_stretch.addWidget(self._profile_description)
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
        selected_profile = self.pick_profile.currentData()
        self._profile_description.setText(selected_profile.description)
        selected_profile.show()
        self._adjust_size_timer.start()

    def _populate_profile_options(self):
        """Add profile options to profile selection."""
        for profile_editor in self._profile_editors.values():
            name = profile_editor.name.capitalize() +' profile'
            self.pick_profile.addItem(name, userData=profile_editor)


class EnergyStepProfileShapeEditor(qtw.QWidget):
    """Base class for step profiles."""

    _subclasses = {}
    description = None
    name = None

    def __init_subclass__(cls, **kwargs):
        """Register subclasses."""
        super().__init_subclass__(**kwargs)
        if not cls.name:
            raise ValueError(f'{type(cls).__name__} must have a '
                             '"name" class attribute.')
        if not cls.description:
            raise ValueError(f'{type(cls).__name__} must have a '
                             '"description" class attribute.')
        cls._subclasses[cls.name] = cls

    def __init__(self):
        """Initialise object."""
        super().__init__()
        self.profile = ()

    @property
    def subclasses(self):
        """Return a {name: cls} dict of subclasses."""
        return self._subclasses.copy()

    def set_profile(self, profile=None):
        """Set displayed profile.

        Parameters
        ----------
        profile : object, optional
            An object that describes the profile that should be set.

        Returns
        -------
        None.
        """
        raise NotImplementedError(
            'set_profile is abstract and must be overridden in subclasses.'
            )

    def update_profile(self):
        """Set the profile to the selected values.

        Subclasses must use this method to modify self.profile.

        Returns
        -------
        None.
        """
        raise NotImplementedError(
            'update_profile is abstract and must be overridden in subclasses.'
            )


class AbruptEnergyStepEditor(EnergyStepProfileShapeEditor):
    """Abrupt step."""

    name = 'abrupt'
    description = ('An abrupt energy step that immediately goes from\n'
                   'the current energy to the next desired energy.')

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


class LinearEnergyStepEditor(EnergyStepProfileShapeEditor):
    """Editor for selecting settings of a linear energy profile."""

    name = 'linear'
    description = ('A linear energy step that goes from the \ncurrent '
                   'energy to the next desired energy\nin equidistant '
                   'intermediate steps.')

    def __init__(self):
        """Initialise object."""
        super().__init__()
        self._controls = {
            'n_steps' : CoercingSpinBox(soft_range=(0, MAX_NUM_STEPS)),
            'duration' : CoercingSpinBox(soft_range=(0, MAX_DELAY),
                                         suffix=' ms'),
            }
        self._compose()

    def set_profile(self, profile):
        """Set linear profile.

        Parameters
        ----------
        profile : Sequence
            A tuple or list of a str and two int values. The first value
            is a str, which should say 'linear' in this case. The second
            value is an integer, which will be set as the number of
            steps in the linear profile, the third value is the duration
            of the delay after each intermediate step in milliseconds.

        Returns
        -------
        None.
        """
        if not profile[1] == self.name or len(profile) != 3:
            raise ValueError('Unsuitable settings for a linear profile.')       # TODO: Catch error on the outside.
        self._controls['n_steps'].setValue(profile[1])
        self._controls['duration'].setValue(profile[2])

    def update_profile(self):
        """Set the profile to the selected values."""
        self.profile = (self.name, self._controls['n_steps'].value(),
                        self._controls['duration'].value())
        if any(value == 0 for value in self.profile):
            self.profile = AbruptEnergyStepEditor().profile

    def _compose(self):
        """Place children widgets."""
        layout = qtw.QFormLayout()
        duration_label = qtw.QLabel('Step duration:')
        step_num_label = qtw.QLabel('Nr. of steps:')
        duration_info = ('<nobr>How long to wait until </nobr>'
                         'the next intermediate step.')
        step_num_info = ('<nobr>The number of intermediate steps.</nobr> '
                         f'Cannot be more than {MAX_NUM_STEPS}.')
        layout.addRow(self._make_info_label(duration_label, duration_info),
                      self._controls['n_steps'])
        layout.addRow(self._make_info_label(step_num_label, step_num_info),
                      self._controls['duration'])
        self.setLayout(layout)

    def _make_info_label(self, label, info):
        """Return a widget containing the label and info."""
        container = qtw.QWidget()
        layout = qtw.QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(label)
        layout.addWidget(FieldInfo.for_widget(label, tooltip=info))
        container.setLayout(layout)
        return container


class FractionalEnergyStepEditor(EnergyStepProfileShapeEditor):
    """Editor for the settings of an energy profile with custom steps."""

    name = 'custom'
    description = ('A custom energy step that goes from the \ncurrent '
                   'energy to the next desired energy\nin fractions of '
                   f'{DELTA_E_NAME}.')
    step_count_reduced = qtc.pyqtSignal()

    def __init__(self):
        """Initialise object."""
        super().__init__()
        self._controls = {
            'add_step' : QNoDefaultPushButton('Add'),
            'remove_step' : QNoDefaultPushButton('Remove'),
            }
        # In order to keep track whether a step has properly been removed
        # from the editor, we have to keep track of how many widgets have
        # been deleted: each step is composed of two widgets.
        self._n_widgets_removed = 0
        self._connect()
        self._compose()

    @property
    def n_steps(self):
        """Return the number of intermediate steps."""
        return max(0, self.layout().rowCount() - N_HEADER_ROWS)

    def set_profile(self, profile):
        """Set fractional profile.

        Parameters
        ----------
        profile : Sequence of int
            Items at even indices will be set as a fraction
            of the step height, items at odd indices are the
            delays in milliseconds.

        Returns
        -------
        None.
        """
        while self.n_steps() > 0:
            self._remove_step()
        self.profile = profile
        for fraction, duration in zip(profile[0::2], profile[1::2]):
            self._add_step(fraction, duration)
        self._update_button_states()

    def update_profile(self):
        """Set the profile to the selected values."""
        profile = []
        layout = self.layout()
        # The first two elements in the layout are the add/remove
        # buttons and the labels. After that, each step is a separate
        # item with two widgets.
        for index in range(2, layout.count()):
            item = layout.itemAt(index)
            profile.append(item.itemAt(0).widget().value())
            profile.append(item.itemAt(1).widget().value())
        self.profile = tuple(profile)
        if not self.profile or not any(self.profile):
            self.profile = AbruptEnergyStepEditor().profile

    @qtc.pyqtSlot()
    def _add_step(self, fraction=None, duration=None):
        """Add a step to the fractional step profile."""
        fraction_handler = CoercingDoubleSpinBox(decimals=2, step=0.05)
        duration_handler = CoercingSpinBox(soft_range=(0, MAX_DELAY),
                                           suffix=' ms')
        for value, handler in zip((fraction, duration),
                                  (fraction_handler, duration_handler)):
            if value:
                handler.setValue(value)
            handler.destroyed.connect(self._emit_step_count_reduced)
        # The layout is required to stop the QFormLayout
        # from compressing the fraction_handler.
        layout = qtw.QHBoxLayout()
        layout.addWidget(fraction_handler)
        layout.addWidget(duration_handler)
        self.layout().addRow(layout)
        self._update_button_states()

    def _compose(self):
        """Place children widgets."""
        layout = qtw.QFormLayout()
        layout.addRow(self._compose_buttons())
        layout.addRow(self._compose_labels())
        self.setLayout(layout)

    def _compose_buttons(self):
        """Return a layout of the add and remove buttons."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(self._controls['add_step'])
        layout.addWidget(self._controls['remove_step'])
        self._controls['remove_step'].setEnabled(False)
        return layout

    def _compose_labels(self):
        """Return a layout of the labels."""
        layout = qtw.QHBoxLayout()
        fraction_label = qtw.QLabel('Step fraction')
        layout.addWidget(fraction_label)
        info = ('<nobr>The energies to set given as a fraction</nobr> '
                f'of {DELTA_E_NAME}. Number of steps cannot exceed '
                f'{MAX_NUM_STEPS}. Any value is acceptable. Zero is '
                'equivalent to the current energy and one to the next '
                'energy. A fraction of one does not have to be '
                'explicitly included at the end, as this is added '
                'automatically with the settle time.')
        layout.addWidget(FieldInfo.for_widget(fraction_label, tooltip=info))
        layout.addWidget(qtw.QLabel('Duration'))
        info = ('<nobr>How long to wait until </nobr>'
                'the next intermediate step.')
        layout.addWidget(FieldInfo.for_widget(fraction_label, tooltip=info))
        return layout

    def _connect(self):
        """Connect buttons to meathods."""
        self._controls['add_step'].clicked.connect(self._add_step)
        self._controls['remove_step'].clicked.connect(self._remove_step)

    def _emit_step_count_reduced(self):
        """Emit the step_count_reduced signal once both widgets are deleted."""
        self._n_widgets_removed += 1
        if self._n_widgets_removed < N_COLUMNS:
            return
        self.step_count_reduced.emit()
        self._n_widgets_removed = 0

    @qtc.pyqtSlot()
    def _remove_step(self):
        """Remove a step from the fractional step profile."""
        if self.n_steps == 0:
            return
        self.layout().removeRow(self.layout().rowCount()-1)
        self._update_button_states()

    def _update_button_states(self):
        """Enable/disable add/remove buttons."""
        self._controls['add_step'].setEnabled(self.n_steps < MAX_NUM_STEPS)
        self._controls['remove_step'].setEnabled(self.n_steps > 0)
