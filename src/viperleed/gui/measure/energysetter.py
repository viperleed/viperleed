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

import configparser
from pathlib import Path

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.settings import NoSettingsError
from viperleed.gui.measure.classes.settings import SettingsError
from viperleed.gui.measure.classes.settings import ViPErLEEDSettings
from viperleed.gui.measure.widgets.spinboxes import SteppingDoubleSpinBox
from viperleed.gui.widgets.lib import AllGUIFonts


class EnergySetterErrors(base.ViPErLEEDErrorEnum):
    """Class for errors occurring in the EnergySetter widget."""

    NO_CONTROLLER = (2000, 'No controller configured. Please select one '
                     'from the "Devices" menu using "Select Controller...".')
    CONTROLLER_FILE_MISSING = (2001, 'The controller settings file no '
                               'longer exists:\n{}\nPlease select a new'
                               ' controller from the "Devices" menu.')
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
        # Path to the settings of the controller
        # that is supposed to set the energy.
        self._path = None
        # Controller object that sets the energy.
        self._controller = None
        # Tracks if the EnergySetter is busy setting an energy.
        self._operation_in_progress = False
        # Stores the last requested energy if setting
        # the energy was not possible at that time.
        self._pending_energy = None

        self._timeout_timer = qtc.QTimer(self)
        self._timeout_timer.setSingleShot(True)
        self._timeout_timer.setInterval(5000)

        self._compose()
        self._connect()

    @property
    def is_busy(self):
        """Return whether an operation is in flight."""
        return self._operation_in_progress

    @property
    def path(self):
        """Return the path to the controller settings."""
        return self._path

    @path.setter
    def path(self, path):
        """Set the path to the controller."""
        if path:
            self._path = Path(path)
        else:
            self._path = None

    @property
    def setting_energy(self):
        """Return whether the energy setter is setting energies or not."""
        return self.set_energy.checkState() == qtc.Qt.Checked

    def _compose(self):
        """Set up the user interface."""
        layout = qtw.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        for widget in (self.set_energy, self.energy_input):
            widget.setFont(AllGUIFonts().buttonFont)
            widget.ensurePolished()
            widget.setEnabled(False)
        layout.addWidget(self.set_energy)
        self.energy_input.setDecimals(1)
        self.energy_input.setRange(0.0, 1000.0)
        self.energy_input.setSingleStep(0.5)
        self.energy_input.setValue(0.0)
        self.energy_input.setSuffix(' eV')
        layout.addWidget(self.energy_input)

        self.setLayout(layout)

    def _connect(self):
        """Connect internal signals and slots."""
        self.set_energy.stateChanged.connect(self._on_set_energy_toggled)
        self.energy_input.editingFinished.connect(self._on_energy_changed)
        self.energy_input.stepped.connect(self._on_energy_changed)
        self._timeout_timer.timeout.connect(self._on_timeout)

    def _connect_controller(self, ctrl):
        """Attempt to connect the controller.

        Returns
        -------
        connected : bool
            True if the controller is connected.

        Emits
        -----
        EnergySetterErrors.CONTROLLER_CONNECTION_FAILED
            If connecting the controller failed.
        """
        ctrl.connect_()
        if not ctrl.connected:
            base.emit_error(self,
                            EnergySetterErrors.CONTROLLER_CONNECTION_FAILED)
            return False
        return True

    def _flush(self):
        """Reset on error."""
        self.set_energy.setChecked(False)
        self._operation_in_progress = False
        self._pending_energy = None
        self._timeout_timer.stop()
        self.cleanup_controller()

    def _get_controller(self):
        """Get or create persistent controller instance.

        Returns
        -------
        ctrl : ControllerABC or None
            The controller instance, or None if creation failed.

        Emits
        -----
        EnergySetterErrors.CONTROLLER_LOAD_FAILED
            If making the controller or loading settings failed.
        EnergySetterErrors.CONTROLLER_CONNECTION_FAILED
            If connecting the controller failed.
        """
        # If we already have a controller for this path, reuse it.
        if self._controller is not None:
            if (self._controller.settings.last_file and
                self._controller.settings.last_file == self.path):
                # Same controller, ensure that the settings
                # are ok and it is connected.
                if not self._reload_ctrl_settings():
                    self.cleanup_controller()
                    return None
                return self._controller
            # Different controller, clean up the old one.
            self.cleanup_controller()

        # Create new controller instance.
        try:
            ctrl = self._make_controller()
        except (NoSettingsError, ValueError) as exc:
            base.emit_error(self, EnergySetterErrors.CONTROLLER_LOAD_FAILED,
                            exc)
            return None

        # Connect and store the controller.
        base.safe_connect(ctrl.error_occurred, self._on_error,
                          type=qtc.Qt.UniqueConnection)
        base.safe_connect(ctrl.serial.busy_changed,
                          self._on_ctrl_finished, type=qtc.Qt.QueuedConnection)
        if not self._connect_controller(ctrl):
            return None

        return ctrl

    def _get_controller_address(self, ctrl_cls, ctrl_settings):
        """Return the address of the used controller.

        Parameters
        ----------
        ctrl_cls : type
            Controller class, used for device detection.
        ctrl_settings : ViPErLEEDSettings
            The controller settings.

        Returns
        -------
        address : str
            The address to use for the controller. An empty string
            makes ControllerABC fall back on the value stored in
            ``controller/address``.
        """
        devices = ctrl_cls().list_devices()
        controller_name = ctrl_settings.get('controller', 'device_name',
                                            fallback=None)
        for device in devices:
            detected_name = (device.more.get('name') or '')
            if detected_name and detected_name == controller_name:
                return device.more.get('address', '')
        return ''

    def _make_controller(self):
        """Make and return a new controller.

        Returns
        -------
        ctrl : ControllerABC
            A controller capable of setting the energy.

        Raises
        ------
        NoSettingsError
            If the settings file cannot be read.
        ValueError
            If ctrl_cls_name was not found.
        """
        ctrl_settings = ViPErLEEDSettings()
        ctrl_settings.read(self.path)
        ctrl_cls_name = ctrl_settings.get('controller', 'controller_class')
        ctrl_settings.prepare_aliases(ctrl_cls_name)
        ctrl_cls = base.class_from_name('controller', ctrl_cls_name)
        address = self._get_controller_address(ctrl_cls, ctrl_settings)
        ctrl = ctrl_cls(settings=ctrl_settings, address=address,
                        sets_energy=True)
        return ctrl

    @qtc.pyqtSlot(bool)
    def _on_ctrl_finished(self, busy):
        """Clean up after energy has been set."""
        if busy:
            return
        self._operation_in_progress = False
        self._timeout_timer.stop()

        # If the setter was switched off, the energy must be set to zero.
        if not self.set_energy.isChecked():
            if self._pending_energy == 0.0:
                # The setter was un-toggled while an energy step
                # was in flight. Now set the energy to zero.
                self._pending_energy = None
                self._set_energy(0.0)
                return
            # Energy is zero, disconnect the controller.
            if self._controller:
                self._controller.disconnect_()
            return

        # If a new energy value was queued during the operation, process it.
        if self._pending_energy is not None:
            energy = self._pending_energy
            self._pending_energy = None
            self._set_energy(energy)

    @qtc.pyqtSlot()
    def _on_energy_changed(self):
        """Handle energy value change."""
        if self.set_energy.checkState() != qtc.Qt.Checked:
            return

        if self._operation_in_progress:
            self._pending_energy = self.energy_input.value()
            return

        self._set_energy(self.energy_input.value())

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

    @qtc.pyqtSlot(int)
    def _on_set_energy_toggled(self, state):
        """Handle checkbox state change.

        Parameters
        ----------
        state : qtc.Qt.CheckState
            The checked state of the QCheckBox.

        Emits
        -----
        EnergySetterErrors.NO_CONTROLLER
            If no ctrl path was given.
        EnergySetterErrors.CONTROLLER_FILE_MISSING
            If the ctrl path does not point to a file.
        """
        self.energy_input.setEnabled(state == qtc.Qt.Checked)
        if state != qtc.Qt.Checked:
            # Checkbox unchecked. If an energy step is in flight, defer
            # setting the energy to zero until the current energy step is
            # completed. Otherwise set the energy to zero.
            if self._operation_in_progress:
                self._pending_energy = 0.0
            else:
                self._set_energy(0.0)
            self.energy_input.setValue(0.0)
            return

        if not self.path:
            base.emit_error(self, EnergySetterErrors.NO_CONTROLLER)
            self.set_energy.setChecked(False)
            return

        if not self.path.is_file():
            base.emit_error(self, EnergySetterErrors.CONTROLLER_FILE_MISSING,
                            self.path)
            self.set_energy.setChecked(False)
            return

        self._controller = self._get_controller()

    @qtc.pyqtSlot()
    def _on_timeout(self):
        """Handle timeout.

        Emits
        -----
        EnergySetterErrors.SET_ENERGY_TIMEOUT
            When a timeout happened.
        """
        self._on_error(EnergySetterErrors.SET_ENERGY_TIMEOUT)

    def _reload_ctrl_settings(self):
        """Reload the settings of the controller.

        Returns
        -------
        reload_ok : bool
            True if reloading the settings worked.

        Emits
        -----
        EnergySetterErrors.CONTROLLER_LOAD_FAILED
            If loading settings failed.
        """
        try:
            settings_ok = self._controller.settings.read_again()
        except (SettingsError, configparser.Error):
            settings_ok = False
        if not settings_ok:
            base.emit_error(self,
                EnergySetterErrors.CONTROLLER_LOAD_FAILED,
                'Controller settings corrupted.')
            return False

        ctrl_cls = self._controller.__class__
        address = self._get_controller_address(ctrl_cls,
                                               self._controller.settings)
        if address:
            self._controller.address = address
        if not self._controller.set_settings(self._controller.settings):
            base.emit_error(self,
                EnergySetterErrors.CONTROLLER_LOAD_FAILED,
                'Could not set controller settings.')
            return False
        if not self._connect_controller(self._controller):
            return False
        return True

    def _set_energy(self, energy):
        """Set energy on the controller.

        Parameters
        ----------
        energy : float
            Energy value in eV.

        Returns
        -------
        None.
        """
        self._operation_in_progress = True
        if self._controller is None:
            self._operation_in_progress = False
            self.set_energy.setChecked(False)
            return
        self._timeout_timer.start()
        self._controller.set_energy(energy, 0, trigger_meas=False)

    def cleanup_controller(self):
        """Clean up the persistent controller."""
        if self._controller is None:
            return
        base.safe_disconnect(self._controller.error_occurred,
                             self._on_error)
        base.safe_disconnect(self._controller.serial.busy_changed,
                             self._on_ctrl_finished)
        self._controller.disconnect_()
        self._controller.deleteLater()
        self._controller = None

    def set_enabled(self, enable):
        """Switch enabled status of widgets."""
        enable &= bool(self.path)
        self.set_energy.setEnabled(enable)
        self.energy_input.setEnabled(enable and self.set_energy.isChecked())
