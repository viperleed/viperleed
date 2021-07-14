"""Module measurecontrollerabc of viperleed.?????.
========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-13
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the ViPErinoController class
which gives commands to the ViPErinoSerialWorker class.
"""

from numpy.polynomial.polynomial import Polynomial
# ViPErLEED modules
from viperleed.guilib.measure.controller.controllerabc import ControllerABC


class MeasureController(ControllerABC):
    """Controller class for the ViPErLEED Arduino Micro."""

    def __init__(self, settings=None, port_name='', controls_camera=False):
        """Initialise controller object."""

        super().init(settings=settings, port_name=port_name,
                     controls_camera=controls_camera)

    @abstractmethod
    def set_energy(self, energy, *other_data, **kwargs):
        """Set electron energy on LEED controller.

        This method must be reimplemented in subclasses. The
        reimplementation should take the energy value in eV
        and other optional data needed by the serial interface
        and turn them into a message that can be sent via
        self.serial.send_message(message, *other_messages).

        Conversion from the desired, true electron energy (i.e.,
        the energy that the electrons will have when exiting the
        gun) to the setpoint value to be sent to the controller
        can be done inside this method by calling
        self.true_energy_to_setpoint(energy).

        Parameters
        ----------
        energy : float
            Nominal electron energy in electronvolts
        *other_data : object, optional
            Other information that the controller may need
            to be able to set the energy. This data should
            not contain information about recalibration of
            the energy itself, as this should be done via
            self.true_energy_to_setpoint(energy). This
            list of extra positional arguments will NOT
            be passed from the GUI during normal operation.
            Hence, it can only be used during self-calibration
            of the controller.
        **kwargs
            Other keyword arguments to set the energy. These
            keyword arguments will NOT be passed from the GUI
            during normal operation. Hence, they should only
            be only be used during self-calibration of the
            controller.

        Returns
        -------
        None.
        """
        return

    @abstractmethod
    def prepare_for_measurement(self):
        """Prepare the controller for a measurement.

        This method must be reimplemented in subclasses. The
        reimplementation should take the settings property and
        use it to do all required tasks before a measurement.
        (i.e. calibrating the electronics, selecting channels,
        determining the gain, ...)
        """
        return

    @abstractmethod
    def measure_iv_video(self):
        """Measure a ramp of energies.

        This method must be reimplemented in subclasses. This
        function should take measurements while doing an energy
        ramp. It should measure I0 (the current emitted by the
        filament used to create the electron beam) and take
        pictures of the LEED screen for every energy step.
        """
        return

    @abstractmethod
    def measure_energy_setpoint(self):
        """Measure the energy offset of the LEED electronics.

        This method must be reimplemented in subclasses. The
        reimplementation is supposed to take measurements of
        the voltage applied to the filament accross the full
        uncalibrated energy spectrum. Both the measured output
        energy and the nominal energy must be returned.

        Returns
        -------
        output_energies : array of floats
            Measured energies
        nominal_energies : array of floats
            Energies the controller set
        """
        return

    def calibrate_energy_setpoint(self):
        """Calibrate the energy setpoint of the LEED electronics

        The offset is measured in the measure_energy_setpoint()
        function which returns the output energies and the
        nominal_energies The output energies are then put into 
        relation to the nominal energies using 
        numpy.polynomial.polynomial.Polynomial. The output 
        energies are used as the x-coordinates and the nominal
        energies are used as the y-coordinates. The resulting
        polynome is written into the config file and used to
        calibrate the nominal energy via the 
        true_energy_to_setpoint() function to get the desired
        output.
        """
        domain = ast.literal_eval(settings['energy_calibration']['domain'])
        output_energies, nominal_energies = self.measure_energy_setpoint()
        fit_polynomial = Polynomial.fit(output_energies, nominal_energies, 1,
                                        domain=domain, window=domain)
        # print("[", ", ".join(f"{ci:.20f}" for ci in fit_polynomial.coef), "]",
              # sep="")
        
        