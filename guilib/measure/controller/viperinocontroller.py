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
# ViPErLEED modules
from viperleed.guilib.measure.controller.controllerabc import ControllerABC


class ViPErinoController(ControllerABC):
    """Controller class for the ViPErLEED Arduino Micro."""

    def __init__(self, settings=None, port_name=''):
        """Initialise controller object."""

        super().init(settings=settings, port_name=port_name)

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
        time: integer
        *more_steps: float and integer (alternating)

        Returns
        -------
        energies_and_times: array of integers
        """
        v_ref_dac = self.__settings.getfloat('measurement_settings', 
                                             'v_ref_dac')
        pc_set_voltage = self.__settings.get('available_commands', 
                                             'PC_SET_VOLTAGE')

        dac_out_vs_nominal_energy = 10/1000  # 10V / 1000 eV
        ouput_gain = 4  # Gain of the output stage on board
        conversion_factor = dac_out_vs_nominal_energy * 65536 / (v_ref_dac *
                                                                 ouput_gain)
        
        energies_and_times = [energy, time, *more_steps]
        number_of_steps = len(energies_and_times)/2

        for i in range number_of_steps:
            tmp_energy = self.true_energy_to_setpoint(energies_and_times[2*i])
            tmp_energy = int(round(tmp_energy * conversion_factor))
            if tmp_energy >= 65536:
                tmp_energy = 65535
            if tmp_energy <= 0:
                tmp_energy = 0
            energies_and_times[2*i] = tmp_energy