"""Module settle_time of viperleed.?????.
========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-19
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the determine_settle_time class
which gives commands to the controller classes.
"""

from measurementabc import MeasurementABC

class DetermineSettletime(MeasurementABC):
    """Settle time determination class."""

    def __init__(self, measurement_settings):
        """Class instanstiation."""
        super().__init__(measurement_settings)
        self.__end_energy = self.settings.getfloat('measurement_settings',
                                                   'end_energy')
        self.__delta_energy = self.settings.getfloat('measurement_settings',
                                                     'delta_energy')
        self.__settle_time = 0
        # Settle time has to be 0 for calibration
        self.__counter = 0
        self.__points_to_take = self.settings.getint(
            'measurement_settings', 'continuous_measurement_points')

    def start_next_measurement(self):
        """Set energy and measure.

        Set energy via the primary controller. Once this is done
        the returned about_to_trigger signal will start the
        measurement on all secondary controllers.

        Returns
        -------
        None.
        """
        self.set_LEED_energy(self.current_energy, self.__settling_time)

    def is_finished(self):
        """Check if the full measurement cycle is done.

        If the number of measurement points has been reached the
        cycle is completed. If not, one is added to the counter
        and the next measurement is started.

        Returns
        -------
        bool
        """
        if self.__points_to_take <= self.__counter:
            self.on_finished()
            return True
        self.__counter += 1
        return False

    def on_finished(self):
        """Calculate settling time.

        Takes measured energies and calculates settling time
        from them. Writes settling time to primary controller
        configuration file afterwards."""
        # TODO: Need to do this for multiple steps and take the highest value
        data_points = 0
        step_height = 0.5
        # This step height will be the aimed for height used later in the LEED I(V) video
        int_update_rate = self.primary_controller.settings.getint(
                            'controller', 'update_rate')
        update_rate = self.primary_controller.settings.getint(
                            'adc_update_rate', update_rate)
        measured_energies = self.data_poinrs['HV']

        length = len(measured_energies) - 1
        for i, energy in enumerate(measured_energies):
            if abs(measured_energies[length - i] - self.current_energy) < step_height:
                data_poinrs += 1
            else:
                break
        self.__settle_time = 10^3*data_points/update_rate
        self.primary_controller.settings.set(
            'measurement_settings', 'settle_time', self.__settle_time)

        file_name = ast.literal_eval(
                        self.settings.get('devices', 'primary_controller')
                        )[0]
        with open(file_name, 'w') as configfile:
            self.primary_controller.settings.write(configfile)

# May want to run this class for each step and overwrite settle_time only if it is higher than the previous one
# Need to use continuous mode --> one measurement ever 20 ms at update_rate 4

