"""ViPErino Controller

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-08
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the ViPErinoController class
which gives commands to the ViPErinoSerialWorker class.
"""

from collections import defaultdict

# ViPErLEED modules
from viperleed.guilib.measure.controller.\
     measurecontrollerabc import MeasureController



class ViPErinoController(MeasureController):
    """Controller class for the ViPErLEED Arduino Micro."""

    def __init__(self, settings=None, port_name='', sets_energy=False):
        """Initialise ViPErino controller object.
        
        Initialise prepare_to_dos dictionaries. The key is a
        string that is used to call a function, the value is a
        boolean that is used to determine if the connected
        function has already been called.
        """

        # Initialise dictionaries for the measurement preparation.
        self.begin_prepare_todos['get_hardware'] = True
        self.begin_prepare_todos['calibrate_adcs'] = True
        self.begin_prepare_todos['set_up_adcs'] = True
        if sets_energy == True:
            self.begin_prepare_todos[
                'set_energy'
                # 'set_energy(self.__energies[0], self.__settle_times[0])'
                ] = True
        self.continue_prepare_todos['start_autogain'] = True
        
        self.__adc0_measurement_type = ''
        self.__adc1_measurement_type = ''
        
        self.__hardware = defaultdict()

        super().__init__(settings=settings, port_name=port_name,
                     sets_energy=sets_energy)

    def set_sets_energy(self, energy_setter):
        """Set the serial to controls energy True/False.
        
        Set the boolean and update the begin_prepare_todos
        dictionary accordingly.

        Parameters
        ----------
        energy_setter : bool
            True if the controller sets the energy.
        """
        super().set_sets_energy(energy_setter)
        # TODO: change key
        key = 'set_energy(self.__energies[0], self.__settle_times[0])'
        if self.__sets_energy == True:
            self.begin_prepare_todos[key] = True
        else:
            if key in self.begin_prepare_todos:
                del self.begin_prepare_todos[key]

    def set_energy(self, energy, time, *more_steps):
        """Convert data from gui to usable values for the DAC.

        Take the energy (or energies), get setpoint energy (or
        energies) and convert it to an integer value for the DAC.
        Afterwards send energy and time to the hardware. The
        controller will automatically trigger and start measuring
        after setting the voltage.

        Parameters
        ----------
        energy: float
            True electron energy in electronvolts (i.e., electrons
            coming out of the gun will have this energy)
        time: integer
            Interval in milliseconds that the controller will wait
            before deeming the LEED optics stable at energy.
        *more_steps: float and integer (alternating)
            The first element in each pair is again one energy, the
            second element the time interval to wait.

        Returns
        -------
        energies_and_times: array of integers
        """
        # TODO: multiple steps
        v_ref_dac = self.settings.getfloat('measurement_settings',
                                           'v_ref_dac')

        dac_out_vs_nominal_energy = 10/1000  # 10V / 1000 eV
        output_gain = 4  # Gain of the output stage on board
        conversion_factor = dac_out_vs_nominal_energy * 65536 / (v_ref_dac *
                                                                 output_gain)

        energies_and_times = [energy, time, *more_steps]
        if len(more_steps) % 2 != 0:
            raise TypeError(f"{self.__class__.__name__}.set_energy: "
                            "Number of energy and time steps do not match. "
                            "Expected an even number of arguments, found "
                            f"{len(more_steps) + 2} arguments.")
        number_of_steps = int(len(energies_and_times)/2)

        for i in range(number_of_steps):
            tmp_energy = self.true_energy_to_setpoint(energies_and_times[2*i])
            tmp_energy = int(round(tmp_energy * conversion_factor))
            if tmp_energy >= 65536:
                tmp_energy = 65535
            if tmp_energy <= 0:
                tmp_energy = 0
            energies_and_times[2*i] = tmp_energy

    def start_autogain(self):
        """Determine starting gain.

        Determine gain after setting the starting energy.

        Returns
        -------
        None.
        """
        # TODO: Had issues in the beginning back then on omicron
        # (Gain too high)
        pc_autogain = self.__settings.get('available_commands', 'PC_AUTOGAIN')
        self.__serial.send_message([pc_autogain])

    def set_up_adcs(self):
        """Set up ADCs.

        Set number of measurements to be taken and channels to
        be used. The used channels have to be calibrated beforehand,
        otherwise measurements will not be taken and an error will
        occur.

        Returns
        -------
        None.
        """
        pc_set_up_adcs = self.__settings.get('available_commands',
                                             'PC_SET_UP_ADCS')
        arduino_config = self.__settings['iv_movie_configuration']
        # TODO: How to get proper section? maybe as a string in settings?
        message = []
        measurement_n = int(arduino_config[
                            'measurement_counter']).to_bytes(2, sign)
        message.append(measurement_n[0])
        message.append(measurement_n[1])       
        for key in ('adc0_channel', 'adc1_channel'):
            message.append(int(arduino_config[key]))
        send_to_arduino([pc_set_up_adcs], message)
        
    def get_hardware(self):
        """Get hardware connected to micro controller.

        This function determines the hardware controlled by
        the micro controller. It has to be run at least once
        after starting the controller up. It will tell both
        the micro controller and this controller class which
        hardware is present in the setup. Trying to do
        anything but getting the hardware or resetting the
        controller before executing this function will result
        in an error.

        Upon receiving this command the micro controller will
        first return the firmware version and then the hardware
        configuration. The conversion is already handled in the
        viperleed serial.

        Returns
        -------
        None.
        """
        pc_configuration = self.__settings.get('available_commands',
                                               'PC_CONFIGURATION')
        self.__serial.send_message([pc_configuration])

    def calibrate_adcs(self):
        """Calibrate ADCs.

        Set update rate, select channels for calibration and
        calibrate them. This function has to be run at least once
        after starting the micro controller. It will store all of
        the necessary calibration data on the micro controller
        itself. It is only done for the channels selected in the
        settings. If the controller is trying to measure channels
        that have not been calibrated yet an error will occur.

        Returns
        -------
        None.
        """
        pc_calibration = self.__settings.get('available_commands',
                                             'PC_CALIBRATION')
        arduino_config = self.__settings['iv_movie_configuration']
        # TODO: How to get proper section? maybe as a string in settings?
        message = []
        for key in ('update_rate', 'adc0_channel', 'adc1_channel'):
            message.append(int(arduino_config[key]))
        send_to_arduino([pc_calibration], message)

    def receive_measurements(self, receive):
        """Receive measurements from the serial.

        Append measurements to the according section. Done via the
        settings. The chosen ADC channels will determine which
        value was measured.

        Parameters
        ----------
        receive : list
            Data received from the serial.

        Returns
        -------
        None.
        """
        self.__measurements['nominal_energy'] = self.__energies[-1:]
        self.__measurements[self.__adc0_measurement_type] = receive[0]
        self.__measurements[self.__adc1_measurement_type] = receive[1]
        self.__measurements['temperature'] = receive[2]
        self.data_ready.emit(self.__measurements)
        self.__measurements.clear()

    def measure_now(self):
        """Take a measurement.

        Measure without setting the energy. This function is
        supposed to be used in time resolved measurements and
        by secondary controllers which will not set the enery.

        If the controller already automatically takes a measurement
        after setting an energy it can be a no op.

        Returns
        -------
        None.
        """
        pc_measure_only = self.__settings.get('available_commands',
                                               'PC_MEASURE_ONLY')
        self.__serial.send_message([pc_measure_only])

    def abort(self):
        """Abort current task.

        Abort what the controller is doing right now and
        return to waiting for further instructions.

        Returns
        -------
        None.
        """
        return

