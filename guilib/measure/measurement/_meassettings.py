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

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.widgets.basewidgets import ButtonWithLabel


class StepProfileViewer(ButtonWithLabel):
    """Viewer of current step profile type.

    Shows the current step profile to the user and opens
    the step profile editor with the push of a button.
    """

    profile_selected = qtc.pyqtSignal()

    def __init__(self, settings, **kwargs):
        """Initialize viewer."""
        super().__init__(**kwargs)
        self._settings = settings
        self.notify_ = self.profile_selected
        self._update_profile_label()
        self.set_button_text('Edit')
        self.profile_editor = StepProfileEditor(self._settings, parent=self)
        self.profile_editor.profile_selected.connect(self._update_profile_label)
        self.button.clicked.connect(self.profile_editor.show)

    @qtc.pyqtSlot()
    def _update_profile_label(self):
        """Change label to the selected step profile."""
        profile_type = self._settings.get('measurement_settings',
                                          'step_profile')
        self.set_label_text(profile_type + ' profile')

    def get_(self, *args):
        """"""

    def set_(self, *args):
        """"""


class StepProfileEditor(qtw.QDialog):
    """Editor for setting step profiles.

    Provides a selection of step profiles and allows
    editing of the corresponding settings.
    """

    profile_selected = qtc.pyqtSignal()

    def __init__(self, settings, **kwargs):
        """Initialize editor."""
        super().__init__(**kwargs)
        self._controls = {
            'profile' : qtw.QComboBox(),
            'linear_editor' : LinearStepEditor(),
            'fraction_editor' : FractionalStepEditor(),
            'accept' : qtw.QPushButton('Accept'),
            'cancel' : qtw.QPushButton('Cancel'),
            }
        self._settings = settings
        self._populate_profile_options()
        self.setWindowTitle('Step profile editor')
        self.setWindowFlags(self.windowFlags()
                            & ~qtc.Qt.WindowContextHelpButtonHint)
        self._compose()
        self._connect()

    @property
    def linear_editor(self):
        """Return the linear editor."""
        return self._controls['linear_editor']

    @property
    def fraction_editor(self):
        """Return the fraction editor."""
        return self._controls['fraction_editor']

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
        layout.addWidget(self._controls['profile'])
        layout.addWidget(self.linear_editor)
        layout.addWidget(self.fraction_editor)
        self.linear_editor.hide()
        self.fraction_editor.hide()
        return layout

    def _connect(self):
        """Connect signals."""
        self._controls['profile'].currentTextChanged.connect(
            self._on_profle_selected
            )
        self._controls['accept'].clicked.connect(self.accept)
        self._controls['cancel'].clicked.connect(self.reject)

    def _on_profle_selected(self):
        """Update the displayed step profile editor."""
        self.linear_editor.hide()
        self.fraction_editor.hide()
        profile = self._controls['profile'].currentText()
        if profile == 'Linear profile':
            self.linear_editor.show()
        if profile == 'Fractional profile':
            self.fraction_editor.show()

    def _populate_profile_options(self):
        """Add profile options to profile selection."""
        control_text = ('Abrupt', 'Linear profile', 'Fractional profile')
        # control_options = (None, self.linear_editor, self.fraction_editor)
        # for text, option in zip(control_text, control_options):
            # self._controls['profile'].addItem(text, userData=option)
        for text in control_text:
            self._controls['profile'].addItem(text)

    @qtc.pyqtSlot()
    def accept(self):
        """Store selected profile then accept."""
        #                                                                       TODO: Write profile to settings
        super().accept()


class LinearStepEditor(qtw.QWidget):                                            # TODO: subclass of?
    """Editor for selecting linear profile settings."""

    def __init__(self):
        super().__init__()
        self._controls = {
            'step_number' : qtw.QSpinBox(),
            'duration' : qtw.QSpinBox(),
            }
        self._labels = {
            'step_number' : qtw.QLabel(),
            'duration' : qtw.QLabel(),
            }
        self._compose()

    def _compose(self):
        """Place children widgets."""
        layout = qtw.QVBoxLayout()
        layout.addLayout(self._compose_step_nr_selection())
        layout.addLayout(self._compose_duration_selection())
        self.setLayout(layout)

    def _compose_step_nr_selection(self):
        """Return a layout of the step number selection."""
        layout = qtw.QHBoxLayout()
        self._labels['step_number'].setText('# of intermediate steps:')
        layout.addWidget(self._labels['step_number'])
        layout.addWidget(self._controls['step_number'])
        return layout

    def _compose_duration_selection(self):
        """Return a layout of the duration selection."""
        layout = qtw.QHBoxLayout()
        self._labels['duration'].setText('Intermediate step duration:')
        layout.addWidget(self._labels['duration'])
        layout.addWidget(self._controls['duration'])
        return layout


class FractionalStepEditor(qtw.QWidget):                                        # TODO: subclass of?
    """Editor for selecting fractional profile settings."""

    def __init__(self):
        super().__init__()
        self._controls = {
            'add_step' : qtw.QPushButton(),
            'remove_step' : qtw.QPushButton(),
            }
        self._fractions = []
        self._durations = []
        self._connect()
        self._compose()

    def _add_step(self):
        """Add a step to the fractional step profile."""
        self._fractions.append(qtw.QSpinBox())
        self._durations.append(qtw.QSpinBox())
        self._compose()

    def _compose(self):
        """Place children widgets."""
        layout = qtw.QVBoxLayout()
        layout.addLayout(self._compose_buttons())
        for fraction, duration in zip(self._fractions, self._durations):
            layout.addLayout(self._compose_step(fraction, duration))
        self.setLayout(layout)

    def _compose_buttons(self):
        """Return a layout of the add and remove buttons."""
        layout = qtw.QHBoxLayout()
        self._controls['add_step'].setText('Add')
        layout.addWidget(self._controls['add_step'])
        self._controls['remove_step'].setText('Remove')
        layout.addWidget(self._controls['remove_step'])
        return layout

    def _compose_step(self, fraction, duration):
        layout = qtw.QHBoxLayout()
        layout.addWidget(fraction)
        layout.addWidget(duration)
        return layout

    def _connect(self):
        """Connect buttons to meathods."""
        self._controls['add_step'].clicked.connect(self._add_step)
        self._controls['remove_step'].clicked.connect(self._remove_step)

    def _remove_step(self):
        """Remove a step from the fractional step profile."""
        self._fractions.pop()
        self._durations.pop()
        self._compose()
