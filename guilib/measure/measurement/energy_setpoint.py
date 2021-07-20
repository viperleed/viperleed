"""Module energy_setpoint of viperleed.?????.
========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the measure_energy_setpoint class
which gives commands to the controller classes.
"""

import ast
from numpy.polynomial.polynomial import Polynomial





    def measure_energy_setpoint(self):
        """Measure the energy offset of the LEED electronics.

        This method must be reimplemented in subclasses. The
        reimplementation is supposed to take measurements of
        the voltage applied to the filament across the full
        uncalibrated energy spectrum.

        Returns
        -------
        None.
        """
        # Edit docstring
        # TODO: Write this function

        self.current_energy += self.delta_energy
        if self.current_energy > self.end_energy:
            self.calibrate_energy_setpoint()
            self.cycle_completed.emit()
        return

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
        nominal_energies = self.data_points['nominal_energies']
        measured_energies = self.data_points['measured_energies']
        domain = ast.literal_eval(self.settings['energy_calibration']['domain'])
        fit_polynomial = Polynomial.fit(measured_energies, nominal_energies, 1,
                                        domain=domain, window=domain)
        coefficients = str(list(fit_polynomial.coef))
        self.settings.set('coefficients', coefficients)