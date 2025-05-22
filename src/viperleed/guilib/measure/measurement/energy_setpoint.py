"""Module energy_setpoint of viperleed.guilib.measure.measurement.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the MeasureEnergySetpoint class
which gives commands to the controller classes.
"""
from configparser import NoSectionError, NoOptionError

from numpy.polynomial.polynomial import Polynomial

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.measurement.abc import (MeasurementABC,
                                                      MeasurementErrors)
from viperleed.guilib.measure.datapoints import QuantityInfo
from viperleed.guilib.measure.classes.settings import NotASequenceError


class MeasureEnergySetpoint(MeasurementABC):
    """Energy calibration class."""

    display_name = 'Energy calibration'

    def __init__(self, measurement_settings):
        """Initialise measurement class"""
        super().__init__(measurement_settings)
        self.__old_coefficients = ""
        self.data_points.time_resolved = False

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
    def __delta_energy(self):
        """Return the amplitude of an energy step in eV."""
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug

        # Eventually, this will be an attribute of an energy generator,
        # and it is unclear whether we will actualy need it.
        fallback = 5.0
        if not self.settings:
            return fallback
        try:
            delta = self.settings.getfloat('measurement_settings',
                                           'delta_energy')
        except (TypeError, ValueError, NoSectionError, NoOptionError):
            # Not a float or not present
            delta = fallback
            base.emit_error(self,
                            MeasurementErrors.INVALID_SETTING_WITH_FALLBACK,
                            '', 'measurement_settings/delta_energy', fallback)
        return delta

    @property
    def __end_energy(self):
        """Return the energy (in eV) at which the energy ramp ends."""
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug

        # Eventually, this will be an attribute of an energy generator,
        # and it is unclear whether we will actualy need it.
        fallback = 1000
        if not self.settings:
            return fallback
        try:
            egy = self.settings.getfloat('measurement_settings', 'end_energy')
        except (TypeError, ValueError, NoSectionError, NoOptionError):
            # Not a float or not present
            egy = fallback
            base.emit_error(self,
                            MeasurementErrors.INVALID_SETTING_WITH_FALLBACK,
                            '', 'measurement_settings/end_energy', fallback)
        return egy                                                              # TODO: warn if end == 1000

    def begin_preparation(self):
        """Start preparation for measurements.

        Prepare the controllers for a measurement which starts
        the measurement cycle. This function only performs the
        first part. (Everything that is done before the starting
        energy is set.)

        Emits
        -----
        __preparation_started
            Starts the measurement preparation and carries
            a tuple of energies and times with it.
        """
        if not any(c.measures(QuantityInfo.HV) for c in self.controllers):
            base.emit_error(
                self, MeasurementErrors.INVALID_SETTINGS,
                'devices/primary_controller or devices/secondary_controllers',
                '\nCannot run an energy calibration if no '
                'controller measures the beam energy.'
                )
            return

        self.__old_coefficients = self.primary_controller.settings.get(
            'energy_calibration', 'coefficients', fallback=''
            )

        self.primary_controller.settings.set('energy_calibration',
                                             'coefficients', '(0, 1)')
        super().begin_preparation()

    def start_next_measurement(self):
        """Set energy and measure.

        Set energy via the primary controller. Once this is done
        the returned about_to_trigger signal will start the
        measurement on all secondary controllers.

        Returns
        -------
        None.
        """
        super().start_next_measurement()
        for controller in self.controllers:
            # Necessary to force secondaries into busy,
            # before the primary returns not busy anymore.
            controller.busy = True
        self.set_leed_energy(*self.step_profile,
                             self.current_energy, self.hv_settle_time)

    def _is_finished(self):
        """Check if the full measurement cycle is done.

        If the energy is above the __end_energy the cycle is
        completed. If not, then the delta energy is added
        and the next measurement is started.

        Returns
        -------
        bool
        """
        super()._is_finished()
        if self.current_energy >= self.__end_energy:
            self.calibrate_energy_setpoint()
            return True
        self.current_energy += self.__delta_energy
        return False

    def calibrate_energy_setpoint(self):                                        # TODO: move this to DataPoints?
        """Calibrate the energy setpoint of the LEED electronics

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
            self.data_points.get_energy_resolved_data(QuantityInfo.HV)
            )
        measured_energies = []
        for ctrl, measurements in data.items():
            # TODO: loop over all controllers, join data. Then
            # calculate polynomial and plot in the future.
            measured_energies = measurements[QuantityInfo.HV]
            break

        primary = self.primary_controller
        try:
            domain = primary.settings.getsequence('energy_calibration',
                                                  'domain',
                                                  fallback=(-10, 1100))
        except NotASequenceError:
            domain = (-10, 1100)
            base.emit_error(self,
                            MeasurementErrors.INVALID_SETTING_WITH_FALLBACK,
                            "", 'energy_calibration/domain', domain)

        fit_polynomial = Polynomial.fit(measured_energies, nominal_energies,
                                        1, domain=domain, window=domain)

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
            base.emit_error(self, MeasurementErrors.INVALID_SETTINGS,
                            'devices/primary_controller', '')
            return

        with open(file_name, 'w', encoding='utf-8') as configfile:
            primary.settings.write(configfile)

    def abort(self):
        """Abort all current actions."""
        if (hasattr(self, "__old_coefficients")
                and self.__old_coefficients):
            self.primary_controller.settings.set('energy_calibration',
                                                 'coefficients',
                                                 self.__old_coefficients)
        super().abort()
