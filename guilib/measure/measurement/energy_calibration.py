"""Module energy_setpoint of viperleed.guilib.measure.measurement.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the MeasureEnergyCalibration class
which acquires a measurement for calibrating the (linear) relationship
between "energy the user asks" and "values that have to be set in the
primary controller".
"""
from configparser import NoSectionError, NoOptionError

from numpy.polynomial.polynomial import Polynomial
from PyQt5 import QtCore as qtc

from viperleed.guilib.measure.classes.abc import QObjectSettingsErrors
from viperleed.guilib.measure.classes.datapoints import QuantityInfo
from viperleed.guilib.measure.classes.settings import NotASequenceError
from viperleed.guilib.measure.measurement.abc import MeasurementABC
from viperleed.guilib.measure.measurement.abc import MeasurementErrors


_MEASURED_EGY = QuantityInfo.HV


class MeasureEnergyCalibration(MeasurementABC):
    """Energy calibration class."""

    display_name = 'Energy calibration'

    def __init__(self, measurement_settings):
        """Initialise instance from settings."""
        super().__init__(measurement_settings)
        self._old_coefficients = ""

    @property
    def start_energy(self):
        """Return the first energy for the energy ramp.

        The returned value is limited below by a minimum energy
        (as found in 'measurement_settings/min_energy' if present,
        5 eV otherwise). This is useful to avoid calibrating for
        the non-linearity of LEED electronics in the low-energy
        regime.

        Returns
        -------
        start_energy : float
            The first energy of the energy ramp.
        """
        start_e = super().start_energy

        # pylint: disable=redefined-variable-type
        # Seems a pylint bug
        try:
            min_e = self.settings.getfloat('measurement_settings',
                                           'min_energy', fallback=5.0)
        except (TypeError, ValueError):
            # Not a float
            min_e = 5.0
        return max(min_e, start_e)

    @property
    def _delta_energy(self):
        """Return the amplitude of an energy step in eV."""
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug

        # Eventually, this will be an attribute of an energy generator,
        # and it is unclear whether we will actually need it.
        fallback = 5.0
        if not self.settings:
            return fallback
        try:
            delta = self.settings.getfloat('measurement_settings',
                                           'delta_energy')
        except (TypeError, ValueError, NoSectionError, NoOptionError):
            # Not a float or not present
            delta = fallback
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
                '', 'measurement_settings/delta_energy', fallback
                )
        return delta

    @property
    def _end_energy(self):
        """Return the energy (in eV) at which the energy ramp ends."""
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug

        # Eventually, this will be an attribute of an energy generator,
        # and it is unclear whether we will actually need it.
        fallback = 1000
        if not self.settings:
            return fallback
        try:
            egy = self.settings.getfloat('measurement_settings', 'end_energy')
        except (TypeError, ValueError, NoSectionError, NoOptionError):
            # Not a float or not present
            egy = fallback
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
                '', 'measurement_settings/end_energy', fallback
                )
        return egy                                                              # TODO: warn if end == 1000

    @qtc.pyqtSlot()
    def abort(self):
        """Abort all current actions."""
        if (hasattr(self, "_old_coefficients")
                and self._old_coefficients):
            self.primary_controller.settings.set('energy_calibration',
                                                 'coefficients',
                                                 self._old_coefficients)
        super().abort()

    def are_runtime_settings_ok(self):
        """Return whether runtime settings are ok.

        Checks whether at least one controller measures the beam energy
        and whether the energy range and sample size is large enough to
        perform any kind of energy calibration.

        Returns
        -------
        settings_ok : bool
            True if the runtime settings are
            sufficient to start a measurement.

        Emits
        -----
        error_occurred
            If the primary controller is missing or if the settings
            are insufficient to perform an energy calibration.
        """
        if not super().are_runtime_settings_ok():
            return False
        self._old_coefficients = self.primary_controller.settings.get(
            'energy_calibration', 'coefficients', fallback=''
            )

        self.primary_controller.settings.set('energy_calibration',
                                             'coefficients', '(0, 1)')

        if not any(c.measures(_MEASURED_EGY) for c in self.controllers):
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTINGS,
                'devices/primary_controller or devices/secondary_controllers',
                '\nCannot run an energy calibration if no '
                'controller measures the beam energy.'
                )
            return False

        egy_range = self._end_energy - self.start_energy
        n_steps = 1 + round(egy_range/self._delta_energy)

        if egy_range < 10:
            # Require at least 10 eV for a reasonable calibration
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTINGS,
                'measurement_settings/start_energy and /end_energy',
                f"\nToo small energy range ({abs(egy_range)} eV) for "
                "calibration. It should be at least 10 eV."
                )
            return False

        if n_steps < 10:
            # Require at least 10 data points for a decent fit
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTINGS,
                'measurement_settings/start_energy, /end_energy, '
                'and /delta_energy',
                f"\nToo few energies ({n_steps}) for a reasonable fit "
                "of the calibration curve. Expected at least 10 energies."
                )
            return False

        return True

    def begin_next_energy_step(self):
        """Set energy and measure.

        Set energy via the primary controller. Once this is done
        the returned about_to_trigger signal will start the
        measurement on all secondary controllers.

        Returns
        -------
        None.
        """
        super().begin_next_energy_step()
        for device in self.devices:
            # Make all controllers busy, so we do not risk going to the
            # next energy step too early: the secondary controllers may
            # be not yet busy when the primary becomes not busy
            device.busy = True
        self.set_leed_energy(*self.step_profile,
                             self.current_energy, self.hv_settle_time)

    def calibrate_energy_setpoint(self):                                        # TODO: move this to DataPoints?
        """Calibrate the energy setpoint of the LEED electronics.

        The offset is measured in the measure_energy_setpoint()
        function which returns the measured energies and the
        nominal energies The measured energies are then put into
        relation to the nominal energies using
        numpy.polynomial.polynomial.Polynomial. A polynomial of
        first degree is most likely accurate enough for the
        calibration, but the degree can be adjusted by changing
        the integer value in the Polynomial.fit function.

        The measured energies are used as the x-coordinates
        and the nominal energies are used as the y-coordinates.
        The resulting polynomial is written into the config file
        and used to calibrate the nominal energy via the
        true_energy_to_setpoint() function to get the desired
        output.

        Returns
        -------
        None.
        """
        data, nominal_energies = (
            self.data_points.get_energy_resolved_data(_MEASURED_EGY)
            )
        measured_energies = []
        for ctrl, measurements in data.items():
            # TODO: loop over all controllers, join data. Then
            # calculate polynomial and plot in the future.
            measured_energies = measurements[_MEASURED_EGY]
            break

        primary = self.primary_controller
        try:
            domain = primary.settings.getsequence('energy_calibration',
                                                  'domain',
                                                  fallback=(-10, 1100))
        except NotASequenceError:
            domain = (-10, 1100)
            self.emit_error(
                QObjectSettingsErrors.INVALID_SETTING_WITH_FALLBACK,
                '', 'energy_calibration/domain', domain
                )

        fit_polynomial, (residuals, *_) = (
            Polynomial.fit(measured_energies, nominal_energies, deg=1,
                           domain=domain, window=domain, full=True)
            )

        # Check quality of fit
        offs, gain = fit_polynomial.coef
        if (abs(offs) > 50               # Can't be more than 50 eV off
            or gain < 1e-3 or gain > 50  # Gain should be >0 and small
                or residuals > 100):     # Most of the time < 1
            self.emit_error(
                MeasurementErrors.RUNTIME_ERROR,
                "Energy calibration fit failed, or fit coefficients "
                f"are poor:\n\nCoefficients=({offs:.2f}, {gain:.2f}) "
                f"[expected ~(0, 1)]\nSum of residuals={residuals[0]:.2f}"
                "\n\nCalibration curve will not be saved."
                )
            return

        primary.settings.set('energy_calibration', 'coefficients',
                             str(list(fit_polynomial.coef)))
        primary.settings.set('energy_calibration', 'domain', str(domain))

        # TODO: may not want to save this to file yet and ask user.
        try:
            primary.settings.update_file()
        except RuntimeError:
            pass
        else:
            return

        try:
            file_name, _ = self.settings.getsequence('devices',
                                                     'primary_controller')
        except (NotASequenceError, ValueError):
            # Something wrong with the settings?
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                            'devices/primary_controller', '')
            return

        with open(file_name, 'w', encoding='utf-8') as configfile:
            primary.settings.write(configfile)

    def get_settings_handler(self):
        """Return a SettingsHandler object for displaying settings."""
        handler = super().get_settings_handler()
        min_energy = self.settings.getfloat('measurement_settings',
                                            'min_energy', fallback=5.)
        option = handler['measurement_settings']['start_energy']
        option.handler_widget.soft_minimum = min_energy
        option.set_info_text(
            '<nobr>The energy at which the measurement starts.</nobr> The '
            f'minimum start energy is {min_energy} eV.'
            )
        return handler

    def set_settings(self, new_settings):
        """Change settings of the measurement."""
        settings_ok = super().set_settings(new_settings)
        self.data_points.time_resolved = False
        return settings_ok

    def _is_finished(self):
        """Check if the full measurement cycle is done.

        If the energy is above the _end_energy the cycle is
        completed. If not, then the delta energy is added
        and the next measurement is started.

        Returns
        -------
        bool
        """
        super()._is_finished()
        if self.current_energy + self._delta_energy > self._end_energy:
            self.calibrate_energy_setpoint()
            return True
        self.current_energy += self._delta_energy
        return False

    def _make_cameras(self):
        """Make sure we have no camera."""
        self.settings.set('devices', 'cameras', '()')
