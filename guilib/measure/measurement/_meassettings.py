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

from viperleed.guilib.measure.widgets.fieldinfo import FieldInfo
from viperleed.guilib.measure.widgets.spinboxes import CoercingDoubleSpinBox
from viperleed.guilib.measure.widgets.spinboxes import CoercingSpinBox
from viperleed.guilib.widgets.basewidgets import ButtonWithLabel


class StepProfileViewer(ButtonWithLabel):
    """Viewer of the current step profile type.

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
        self._delay = qtc.QTimer()
        self._delay.setSingleShot(True)
        self._delay.setInterval(10)
        self._delay.timeout.connect(self._connect_finished)
        self._delay.start()

    def _connect_finished(self):
        """Connect finished signal after full instantiation."""
        self.window().finished.connect(self.profile_editor.reject)

    def _connect(self):
        """Connect (only once) relevant signals and slots."""
        self.button.clicked.connect(self.profile_editor.show)
        self.profile_editor.accepted.connect(self._on_settings_changed)

    @qtc.pyqtSlot()
    def _on_settings_changed(self):
        """Update step profile to selected profile."""
        if isinstance(self.profile_editor.profile[0], str):
            chosen_profile = self.profile_editor.profile[0]
        else:
            chosen_profile = 'fractional'
        self.set_label_text(chosen_profile + ' profile')
        self.settings_changed.emit()

    def get_(self):
        """Return the value to be stored in the config."""
        return str(self.profile_editor.profile)

    def set_(self, value):
        """Set label and load profile into step profile editor."""
        value = literal_eval(value)
        if not value:
            value = ('abrupt', )
        if isinstance(value, str):
            value = (value,)
        self.profile_editor.profile = value
        if isinstance(value[0], str):
            self.set_label_text(value[0] + ' profile')
            return
        self.set_label_text('fractional profile')


class StepProfileEditor(qtw.QDialog):
    """Editor for setting step profiles.

    Provides a selection of step profiles and allows
    editing of the corresponding settings.
    """

    def __init__(self, **kwargs):
        """Initialize editor."""
        super().__init__(**kwargs)
        self._controls = {
            'profile' : qtw.QComboBox(),
            'abrupt' : AbruptStep(),
            'linear_editor' : LinearStepEditor(),
            'fraction_editor' : FractionalStepEditor(),
            'accept' : qtw.QPushButton('Accept'),
            'cancel' : qtw.QPushButton('Cancel'),
            }
        self._populate_profile_options()
        self.setWindowTitle('Step profile editor')
        self.setWindowFlags(self.windowFlags()
                            & ~qtc.Qt.WindowContextHelpButtonHint)
        self._delay = qtc.QTimer()
        self._delay.setSingleShot(True)
        self._delay.setInterval(3)
        self._delay.timeout.connect(self.adjustSize)
        self._compose()
        self._connect()

    @property
    def fraction_editor(self):
        """Return the fraction editor."""
        return self._controls['fraction_editor']

    @property
    def linear_editor(self):
        """Return the linear editor."""
        return self._controls['linear_editor']

    @property
    def profile(self):
        """Return the currently selected profile."""
        return self._controls['profile'].currentData().profile

    @profile.setter
    def profile(self, profile_data):
        """Set the profile from settings."""
        if isinstance(profile_data[0], str):
            if profile_data[0] == 'linear':
                index = self._controls['profile'].findText('Linear profile')
            else:
                index = self._controls['profile'].findText('Abrupt')
        else:
            index = self._controls['profile'].findText('Fractional profile')
        self._controls['profile'].setCurrentIndex(index)
        self._controls['profile'].currentData().set_profile(profile_data)

    def _compose(self):
        """Compose StepProfileEditor."""
        layout = qtw.QVBoxLayout()
        layout.addLayout(self._compose_editor())
        layout.addLayout(self._compose_accept_cancel())
        self.setLayout(layout)

    def _compose_accept_cancel(self):
        """Compose accept and cancel buttons."""
        layout = qtw.QHBoxLayout()
        layout.addWidget(self._controls['cancel'])
        layout.addStretch(1)
        layout.addWidget(self._controls['accept'])
        return layout

    def _compose_editor(self):
        """Compose editor section."""
        layout = qtw.QHBoxLayout()
        profile_with_stretch = qtw.QVBoxLayout()
        profile_with_stretch.addWidget(self._controls['profile'])
        profile_with_stretch.addStretch(1)
        layout.addLayout(profile_with_stretch)
        layout.addWidget(self.linear_editor)
        layout.addWidget(self.fraction_editor)
        self.linear_editor.hide()
        self.fraction_editor.hide()
        return layout

    def _connect(self):
        """Connect signals."""
        self._controls['profile'].currentTextChanged.connect(
            self._on_profile_selected
            )
        self._controls['accept'].clicked.connect(self.accept)
        self._controls['cancel'].clicked.connect(self.reject)
        self.fraction_editor.step_count_reduced.connect(self._delay.start)

    def _on_profile_selected(self):
        """Update the displayed step profile editor."""
        self.linear_editor.hide()
        self.fraction_editor.hide()
        self._controls['profile'].currentData().show()
        self._delay.start()

    def _populate_profile_options(self):
        """Add profile options to profile selection."""
        control_text = ('Abrupt', 'Linear profile', 'Fractional profile')
        control_options = (self._controls['abrupt'], self.linear_editor,
                           self.fraction_editor)
        for text, option in zip(control_text, control_options):
            self._controls['profile'].addItem(text, userData=option)

    @qtc.pyqtSlot()
    def accept(self):
        """Store selected profile then accept."""
        self._controls['profile'].currentData().update_profile()
        super().accept()


class ProfileStep(qtw.QWidget):
    """Base class for step profiles."""

    def __init__(self):
        """Initialise object."""
        super().__init__()
        self.profile = ()


class AbruptStep(ProfileStep):
    """Abrupt step."""

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


class LinearStepEditor(ProfileStep):
    """Editor for selecting linear profile settings."""

    def __init__(self):
        """Initialise object."""
        super().__init__()
        self._controls = {
            'step_number' : CoercingSpinBox(range_=(0, 32767)),
            'duration' : CoercingSpinBox(range_=(0, 9999), suffix=' ms'),
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
        """Return a layout of the step number selection."""
        layout = qtw.QHBoxLayout()
        step_number_label = qtw.QLabel()
        step_number_label.setText('Nr. of steps:')
        layout.addWidget(step_number_label)
        size = step_number_label.fontMetrics().boundingRect('a').height()
        info = ('The number of intermediate steps.')
        layout.addWidget(FieldInfo(info, size=size))
        layout.addWidget(self._controls['step_number'])
        return layout

    def _compose_duration_selection(self):
        """Return a layout of the duration selection."""
        layout = qtw.QHBoxLayout()
        duration_label = qtw.QLabel()
        duration_label.setText('Step duration:')
        layout.addWidget(duration_label)
        size = duration_label.fontMetrics().boundingRect('a').height()
        info = ('The settle time after each energy.')
        layout.addWidget(FieldInfo(info, size=size))
        layout.addWidget(self._controls['duration'])
        return layout

    def set_profile(self, profile):
        """Set linear profile."""
        self._controls['step_number'].setValue(profile[1])
        self._controls['duration'].setValue(profile[2])

    def update_profile(self):
        """Set the profile to the selected values."""
        self.profile = ('linear', self._controls['step_number'].value(),
                        self._controls['duration'].value())
        if self.profile == ('linear', 0, 0):
            self.profile = ('abrupt',)


class FractionalStepEditor(ProfileStep):
    """Editor for selecting fractional profile settings."""

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
        fraction_handler = CoercingDoubleSpinBox(decimals=2, range_=(0, 2),
                                                 step=0.05)
        duration_handler = CoercingSpinBox(range_=(0, 32767), suffix=' ms')
        for value, handler in zip((fraction, duration),
                                  (fraction_handler, duration_handler)):
            if value:
                handler.setValue(value)
            handler.destroyed.connect(self._emit_step_count_reduced)
            layout.addWidget(handler)
        self._steps.append(layout)
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
        return layout

    def _compose_labels(self):
        """Return a layout of the labels."""
        layout = qtw.QHBoxLayout()
        fraction_label = qtw.QLabel()
        fraction_label.setText('Step fraction')
        layout.addWidget(fraction_label)
        size = fraction_label.fontMetrics().boundingRect('a').height()
        info = ('<nobr>The energies to set given as a '
                'fraction</nobr> of "Delta energy".')
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
        for widget_index in range(item.layout().count()):
            item.itemAt(widget_index).widget().deleteLater()

    def set_profile(self, profile):
        """Set fractional profile."""
        while self.layout().count() >= 4:
            self._remove_step()
        self.profile = profile
        for fraction, duration in zip(profile[0::2], profile[1::2]):
            self._add_step(fraction, duration)

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
