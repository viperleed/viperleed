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
from viperleed.guilib.measure.controller.\
     measurecontrollerabc import MeasureController, MeasureControllerError



class ViPErinoController(MeasureController):
    """Controller class for the ViPErLEED Arduino Micro."""
    
    plot_info = MeasureController.plot_info.copy()
    plot_info['temperature'] = ('°C', 'lin')
    plot_info['cold_junction'] = ('°C', 'lin')


    def __init__(self, settings=None, port_name='', controls_camera=False):
        """Initialise controller object."""
        
        # TODO: Add stuff here if needed, otherwise remove init.
        
        super().__init__(settings=settings, port_name=port_name,
                     controls_camera=controls_camera)

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
        v_ref_dac = self.settings.getfloat('measurement_settings',
                                           'v_ref_dac')

        dac_out_vs_nominal_energy = 10/1000  # 10V / 1000 eV
        output_gain = 4  # Gain of the output stage on board
        conversion_factor = dac_out_vs_nominal_energy * 65536 / (v_ref_dac *
                                                                 ouput_gain)

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
