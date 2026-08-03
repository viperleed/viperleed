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


class EnergySetter(qtw.QWidget):
    """Widget for setting LEED energy without data acquisition.

    This widget provides a checkbox and energy input field that allow
    users to set the beam energy on the primary controller without
    starting a measurement.

    """
    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, **kwargs):
        """Initialize the EnergySetter widget."""
        super().__init__(**kwargs)
        self.energy_input = SteppingDoubleSpinBox()
        self.set_energy = qtw.QCheckBox('Set energy')
        self.primary_path = ''
        self._primary_controller = None

        self._compose()
        self._connect()

    @property
    def setting_energy(self):
        """Return whether the energy setter is setting energies or not."""
        return self.set_energy.checkState() == qtc.Qt.Checked

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

    @qtc.pyqtSlot(int)
    def _on_set_energy_toggled(self, state):
        """Handle checkbox state change."""
        if state != qtc.Qt.Checked:
            # Checkbox unchecked, set energy to zero before cleaning up.
            if self.primary_path:
                primary_path = Path(self.primary_path)
                if primary_path.is_file():
                    self._set_energy(primary_path, 0.0)
        self.cleanup_primary_controller()
        # When checked, controller will be created on first energy change.

    @qtc.pyqtSlot()
    @qtc.pyqtSlot(float)
    def _on_energy_changed(self, *_):
        """Handle energy value change."""
        if self.set_energy.checkState() != qtc.Qt.Checked:
            return

        if not self.primary_path:
            qtw.QMessageBox.warning(
                self, 'No Controller Available',
                'No primary controller configured. Please select one from '
                'the "Devices" menu using "Select Primary Controller...".'
            )
            self.set_energy.setChecked(False)
            return

        primary_path = Path(self.primary_path)
        if not primary_path.is_file():
            qtw.QMessageBox.warning(
                self, 'Controller File Missing',
                f'The primary controller settings file no longer exists:\n{primary_path}\n'
                'Please select a new primary controller from the "Devices" menu.'
            )
            self.set_energy.setChecked(False)
            return

        energy = self.energy_input.value()
        self._set_energy(primary_path, energy)

    def _get_primary_controller(self, primary_path):
        """Get or create persistent primary controller instance.

        Parameters
        ----------
        primary_path : Path
            Path to the controller settings file.

        Returns
        -------
        ControllerABC or None
            The controller instance, or None if creation failed.
        """
        # If we already have a controller for this path, reuse it.
        if self._primary_controller is not None:
            ctrl_settings = self._primary_controller.settings
            if (ctrl_settings.last_file and
                ctrl_settings.last_file == primary_path):
                # Same controller, ensure it's connected.
                if not self._primary_controller.connected:
                    self._primary_controller.connect_()
                return self._primary_controller
            else:
                # Different controller, clean up the old one.
                self.cleanup_primary_controller()

        # Create new controller instance.
        try:
            ctrl_settings = ViPErLEEDSettings()
            ctrl_settings.read(primary_path)
            ctrl_cls_name = ctrl_settings.get('controller', 'controller_class')
            ctrl_settings.prepare_aliases(ctrl_cls_name)
            ctrl_cls = base.class_from_name('controller', ctrl_cls_name)
            address = ctrl_settings.get('controller', 'address')
            primary_ctrl = ctrl_cls(settings=ctrl_settings, address=address,
                                    sets_energy=True)
        except (NoSettingsError, ValueError, KeyError) as err:
            qtw.QMessageBox.warning(
                self, 'Failed to Load Controller',
                f'Could not load the last used controller:\n{err}'
            )
            return None

        # Connect and store the controller.
        primary_ctrl.connect_()
        if not primary_ctrl.connected:
            qtw.QMessageBox.warning(
                self, 'Connection Failed',
                'Could not connect to the controller. '
                'Please check that the device is available.'
            )
            primary_ctrl.deleteLater()
            return None

        self._primary_controller = primary_ctrl
        return primary_ctrl

    def cleanup_primary_controller(self):
        """Clean up the persistent primary controller."""
        if self._primary_controller is not None:
            self._primary_controller.disconnect_()
            self._primary_controller.deleteLater()
            self._primary_controller = None

    def _set_energy(self, ctrl_path, energy):
        """Set energy on the primary controller.

        Parameters
        ----------
        ctrl_path : Path
            Path to the controller settings file.
        energy : float
            Energy value in eV.

        Returns
        -------
        None.
        """
        primary_ctrl = self._get_primary_controller(ctrl_path)
        if primary_ctrl is None:
            self.set_energy.setChecked(False)
            return

        # Create event loop to wait for completion
        loop = qtc.QEventLoop()

        def on_finished():
            loop.quit()

        def on_error(error_info):
            self.set_energy.setChecked(False)
            loop.quit()
            self.error_occurred.emit(error_info)
            qtw.QMessageBox.warning(
                self, 'Error',
                f'Failed to set energy:\n{error_info}'
            )

        # Connect to serial's connection_changed to know when done
        primary_ctrl.error_occurred.connect(on_error)

        # Set up timeout in case no response arrives
        timeout_timer = qtc.QTimer()
        timeout_timer.setSingleShot(True)
        timeout_timer.setInterval(5000)  # 5 second timeout
        timeout_timer.timeout.connect(lambda: (
            self.set_energy.setChecked(False),
            qtw.QMessageBox.warning(
                self, 'Timeout',
                'No response from controller. Check connection.'
            ),
            loop.quit()
        ))

        # Connect to busy_changed - when it goes False, operation is complete
        primary_ctrl.serial.busy_changed.connect(
            lambda busy: (
                timeout_timer.stop(),
                on_finished()
            ) if not busy else None,
            type=qtc.Qt.QueuedConnection
        )

        timeout_timer.start()
        primary_ctrl.set_energy(energy, 0, trigger_meas=False)

        loop.exec()

    def set_enabled(self, enable):
        """Switch enabled status of buttons."""
        self.set_energy.setEnabled(enable and bool(self.primary_path))
        self.energy_input.setEnabled(enable and bool(self.primary_path))
