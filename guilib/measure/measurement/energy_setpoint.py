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
        """Initialise measurement class.

        This is an upgraded version of its parent class.
        __end_energy, __delta_energy and __hv_settle_time are
        read from the settings and made private properties.
        """
        super().__init__(measurement_settings)
        self.__end_energy = 0                                                   # TODO: all attributes to properties
        self.__delta_energy = 1
        self.__hv_settle_time = 0
        if self.settings:
            self.__delta_energy = self.settings.getfloat(                       # TODO: float exceptions
                'measurement_settings', 'delta_energy', fallback=1
                )
            self.__end_energy = self.settings.getfloat(                         # TODO: float exceptions
                'measurement_settings', 'end_energy', fallback=0
                )

        if self.primary_controller:
            self.__hv_settle_time = self.primary_controller.hv_settle_time

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
        if not self.settings:
            return 0.0
        try:
            start_e = self.settings.getfloat('measurement_settings',
                                             'start_energy', fallback=0.0)
        except (TypeError, ValueError):
            # Not a float
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                            'measurement_settings/start_energy')
            start_e = 0.0

        try:
            min_e = self.settings.getfloat('measurement_settings',
                                           'min_energy', fallback=5.0)
        except (TypeError, ValueError):
            # Not a float
            min_e = 5.0
        return max(min_e, start_e)

    def begin_measurement_preparation(self):
        """Start preparation for measurements.

        Prepare the controllers for a measurement which starts
        the measurement cycle. This function only performs the
        first part. (Everything that is done before the starting
        energy is set.)

        Returns
        -------
        None.

        Emits
        -----
        begin_preparation
            Starts the measurement preparation and carries
            a tuple of energies and times with it.
        """
        coefficients = '(0, 1)'
        self.primary_controller.settings.set(
            'energy_calibration', 'coefficients', coefficients)
        super().begin_measurement_preparation()

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
        self.set_LEED_energy(self.current_energy, self.__hv_settle_time)

    def is_finished(self):
        """Check if the full measurement cycle is done.

        If the energy is above the __end_energy the cycle is
        completed. If not, then the delta energy is added
        and the next measurement is started.

        Returns
        -------
        bool
        """
        super().is_finished()
        if self.current_energy >= self.__end_energy:
            self.calibrate_energy_setpoint()
            return True
        self.current_energy += self.__delta_energy
        return False

    def calibrate_energy_setpoint(self):
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
        None
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

        try:
            file_name, _ = self.settings.getsequence('devices',
                                                     'primary_controller')
        except (NotASequenceError, ValueError):
            # Something wrong with the settings?
            base.emit_error(self, MeasurementErrors.INVALID_MEAS_SETTINGS,
                            'devices/primary_controller')
            return

        with open(file_name, 'w') as configfile:
            primary.settings.write(configfile)

    def abort(self):
        """Abort all current actions.

        Abort and reset all variables.

        Returns
        -------
        None.
        """
        coefficients = '(0, 1)'
        if self.primary_controller.settings:
            self.primary_controller.settings.set(
                'energy_calibration', 'coefficients', coefficients)
        super().abort()
