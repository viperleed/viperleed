"""Module energysetter of viperleed.gui.measure.

Defines the EnergySetter widget for setting LEED energy without
starting a measurement.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-07-31'
__license__ = 'GPLv3+'

from pathlib import Path

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.settings import NoSettingsError
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings
from viperleed.gui.measure.widgets.spinboxes import SteppingDoubleSpinBox
from viperleed.gui.widgets.lib import AllGUIFonts


class EnergySetterErrors(base.ViPErLEEDErrorEnum):
    """Class for errors occurring in the EnergySetter widget."""

    NO_PRIMARY_CONTROLLER = (2000, 'No primary controller configured. Please '
                             'select one from the "Devices" menu using '
                             '"Select Primary Controller...".')
    CONTROLLER_FILE_MISSING = (2001, 'The primary controller settings file '
                               'no longer exists:\n{}\nPlease select a new '
                               'primary controller from the "Devices" menu.')
    CONTROLLER_LOAD_FAILED = (2002, 'Could not load the last used '
                              'controller:\n{}')
    CONTROLLER_CONNECTION_FAILED = (2003, 'Could not connect to the '
                                    'controller. Please check that the '
                                    'device is available.')
    SET_ENERGY_TIMEOUT = (2004, 'No response from controller. Check '
                          'connection.')


class EnergySetter(qtw.QWidget):
    """Widget for setting LEED energies without data acquisition.

    This widget provides a checkbox and energy input field that allow
    users to set the beam energy.
    """

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, **kwargs):
        """Initialize the EnergySetter widget."""
        super().__init__(**kwargs)
        self.energy_input = SteppingDoubleSpinBox()
        self.set_energy = qtw.QCheckBox('Set energy')
        self._primary_path = None
        self._primary_controller = None
        self._pending_energy = None
        self._operation_in_progress = False

        self._timeout_timer = qtc.QTimer(self)
        self._timeout_timer.setSingleShot(True)
        self._timeout_timer.setInterval(5000)

        self._compose()
        self._connect()

    @property
    def setting_energy(self):
        """Return whether the energy setter is setting energies or not."""
        return self.set_energy.checkState() == qtc.Qt.Checked

    @property
    def primary_path(self):
        """Return the path to the primary controller settings."""
        return self._primary_path

    @primary_path.setter
    def primary_path(self, primary):
        """Set the path to the primary controller."""
        if primary:
            self._primary_path = Path(primary)
        else:
            self._primary_path = None

    def _compose(self):
        """Set up the user interface."""
        layout = qtw.QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        for widget in (self.set_energy, self.energy_input):
            widget.setFont(AllGUIFonts().buttonFont)
            widget.ensurePolished()
            widget.setEnabled(False)
        layout.addWidget(self.set_energy, 0, 0)
        self.energy_input.setDecimals(1)
        self.energy_input.setRange(0.0, 1000.0)
        self.energy_input.setSingleStep(0.5)
        self.energy_input.setValue(0.0)
        self.energy_input.setSuffix(' eV')
        layout.addWidget(self.energy_input, 0, 1)

        self.setLayout(layout)

    def _connect(self):
        """Connect internal signals and slots."""
        self.set_energy.stateChanged.connect(self._on_set_energy_toggled)
        self.energy_input.editingFinished.connect(self._on_energy_changed)
        self.energy_input.valueChanged.connect(self._on_energy_changed)
        self._timeout_timer.timeout.connect(self._on_timeout)

    @qtc.pyqtSlot(int)
    def _on_set_energy_toggled(self, state):
        """Handle checkbox state change.

        Parameters
        ----------
        state : qtc.Qt.CheckState
            The checked state of the QCheckBox.

        Emits
        -----
        EnergySetterErrors.NO_PRIMARY_CONTROLLER
            If no primary path was given.
        EnergySetterErrors.CONTROLLER_FILE_MISSING
            If the primary path does not point to a file.
        """
        self.energy_input.setEnabled(state == qtc.Qt.Checked)
        if state != qtc.Qt.Checked:
            # Checkbox unchecked, set energy to zero before cleaning up.
            self._set_energy(0.0)
            self.energy_input.setValue(0.0)

        if not self.primary_path:
            base.emit_error(self, EnergySetterErrors.NO_PRIMARY_CONTROLLER)
            self.set_energy.setChecked(False)
            return

        if not self.primary_path.is_file():
            base.emit_error(self, EnergySetterErrors.CONTROLLER_FILE_MISSING,
                            self.primary_path)
            self.set_energy.setChecked(False)
            return

        self._primary_controller = self._get_primary_controller()

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(float)
    def _on_energy_changed(self, *_):
        """Handle energy value change."""
        if self.set_energy.checkState() != qtc.Qt.Checked:
            return

        if self._operation_in_progress:
            self._pending_energy = self.energy_input.value()
            return

        energy = self.energy_input.value()
        self._set_energy(energy)

    def _get_primary_controller(self):
        """Get or create persistent primary controller instance.

        Returns
        -------
        ControllerABC or None
            The controller instance, or None if creation failed.

        Emits
        -----
        EnergySetterErrors.CONTROLLER_LOAD_FAILED
            If making the primary controller failed.
        EnergySetterErrors.CONTROLLER_CONNECTION_FAILED
            If connecting the primary controller failed.
        """
        # If we already have a controller for this path, reuse it.
        if self._primary_controller is not None:
            ctrl_settings = self._primary_controller.settings
            if (ctrl_settings.last_file and
                ctrl_settings.last_file == self.primary_path):
                # Same controller, ensure it's connected.
                if not self._primary_controller.connected:
                    self._primary_controller.connect_()
                return self._primary_controller
            # Different controller, clean up the old one.
            self.cleanup_primary_controller()

        # Create new controller instance.
        try:
            primary_ctrl = self._make_primary_controller()
        except (NoSettingsError, ValueError) as err:
            base.emit_error(self, EnergySetterErrors.CONTROLLER_LOAD_FAILED,
                            err)
            return None

        # Connect and store the controller.
        base.safe_connect(primary_ctrl.error_occurred, self._on_error,
                          type=qtc.Qt.UniqueConnection)
        base.safe_connect(primary_ctrl.serial.busy_changed,
                          self._on_ctrl_finished, type=qtc.Qt.QueuedConnection)
        primary_ctrl.connect_()
        if not primary_ctrl.connected:
            base.emit_error(self,
                            EnergySetterErrors.CONTROLLER_CONNECTION_FAILED)
            self.cleanup_primary_controller()
            return None

        return primary_ctrl

    def _make_primary_controller(self):
        """Make and return a new primary controller.

        Returns
        -------
        primary_ctrl : ControllerABC
            A controller capable of setting the energy.

        Raises
        ------
        NoSettingsError
            If the settings file cannot be read.
        ValueError
            If ctrl_cls_name was not found.
        """
        ctrl_settings = ViPErLEEDSettings()
        ctrl_settings.read(self.primary_path)
        ctrl_cls_name = ctrl_settings.get('controller', 'controller_class')
        ctrl_settings.prepare_aliases(ctrl_cls_name)
        ctrl_cls = base.class_from_name('controller', ctrl_cls_name)
        address = ctrl_settings.get('controller', 'address')
        primary_ctrl = ctrl_cls(settings=ctrl_settings, address=address,
                                sets_energy=True)
        return primary_ctrl

    def cleanup_primary_controller(self):
        """Clean up the persistent primary controller."""
        if self._primary_controller is None:
            return
        base.safe_disconnect(self._primary_controller.error_occurred,
                             self._on_error,)
        base.safe_disconnect(self._primary_controller.serial.busy_changed,
                             self._on_ctrl_finished)
        self._primary_controller.disconnect_()
        self._primary_controller.deleteLater()
        self._primary_controller = None

    @qtc.pyqtSlot(tuple)
    def _on_error(self, error_info):
        """Handle controller errors.

        Parameters
        ----------
        error_info : tuple

        Emits
        -----
        error_occurred
            With error_info on the error that occurred.
        """
        self.error_occurred.emit(error_info)
        self._flush()

    @qtc.pyqtSlot()
    def _on_timeout(self):
        """Handle timeout.

        Emits
        -----
        EnergySetterErrors.SET_ENERGY_TIMEOUT
            When a timeout happened.
        """
        base.emit_error(self, EnergySetterErrors.SET_ENERGY_TIMEOUT)
        self._flush()

    def _flush(self):
        """Reset on error."""
        self.set_energy.setChecked(False)
        self._operation_in_progress = False
        self._pending_energy = None
        self._timeout_timer.stop()
        self.cleanup_primary_controller()

    def _set_energy(self, energy):
        """Set energy on the primary controller.

        Parameters
        ----------
        energy : float
            Energy value in eV.

        Returns
        -------
        None.
        """
        self._operation_in_progress = True
        if self._primary_controller is None:
            self._operation_in_progress = False
            self.set_energy.setChecked(False)
            return
        self._timeout_timer.start()
        self._primary_controller.set_energy(energy, 0, trigger_meas=False)

    @qtc.pyqtSlot(bool)
    def _on_ctrl_finished(self, busy):
        """Clean up after energy has been set."""
        if busy:
            return
        self._operation_in_progress = False
        self._timeout_timer.stop()

        # If set energy is no longer toggled, we want to disconnect the ctrl.
        if not self.set_energy.isChecked():
            self.cleanup_primary_controller()
            return

        # If a new energy value was queued during the operation, process it.
        if self._pending_energy is not None:
            energy = self._pending_energy
            self._pending_energy = None
            self._set_energy(energy)

    def set_enabled(self, enable):
        """Switch enabled status of widgets."""
        enable &= bool(self.primary_path)
        self.set_energy.setEnabled(enable)
        self.energy_input.setEnabled(enable and self.set_energy.isChecked())
