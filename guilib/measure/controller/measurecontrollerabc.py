"""Module measurecontrollerabc of viperleed.?????.
========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-13
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the MeasureController class
which gives commands to the SerialABC class.
"""
from collections import defaultdict

from numpy.polynomial.polynomial import Polynomial
# ViPErLEED modules
from viperleed.guilib.measure.controller.controllerabc import ControllerABC

class MeasureControllerError(hardwarebase.ViPErLEEDErrorEnum):
    SECTIONS_DO_NOT_MATCH = (151,
                             "The sections of the labels library and the "
                             "data_points do not match. Check if "
                             "reimplementation contains all implemented "
                             "sections in both libraries.")


class MeasureController(ControllerABC):
    """Controller class for the ViPErLEED Arduino Micro.
    
    The labels dictionary in this class may be reimplemented
    in subclasses. Each section (label) needs to contain 
    the unit and the scaling ('lin' for linear and 'log' for
    logarithmic scaling.
    
    The data_points dictionary may be reimplemented aswell
    and needs to contain the same sections as the labels
    dictionary.
    """
    
    # Reemplementations may introduce more/other keys.
    self.labels = defaultdict(list)
    self.labels['nominal_energies'] = ('eV', 'lin')
    self.labels['I0'] = ('uA', 'lin')
    self.labels['measured_energies'] = ('eV', 'lin')
    self.labels['elapsed_time'] = ('ms', 'lin')
    self.labels['aux0'] = ('', 'log')
    
    def __init__(self, settings=None, port_name='', controls_camera=False):
        """Initialise controller object."""
        
        # Reemplementations may introduce more/other keys.
        self.data_points = defaultdict(list)
        self.data_points['nominal_energies'] = []
        self.data_points['I0'] = []
        self.data_points['measured_energies'] = []
        self.data_points['elapsed_time'] = []
        self.data_points['aux0'] = []
        
        if data_points.keys() != labels.keys():
            self.error_occurred.emit(
                MeasureControllerError.SECTIONS_DO_NOT_MATCH
                )
                
        
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
        reimplementation should take the settings property
        derived from the configuration file and use it to do 
        all required tasks before a measurement.
        (i.e. calibrating the electronics, selecting channels,
        determining the gain, ...) 
        
        All communication with the hardware should be done 
        via the send_command signal,which takes instructions 
        as a touple and forwards them to the send_message 
        function of the serialabc class.
        
        It should be able to select the update rate of the
        measurement electronics and change channels if there 
        are more than one.
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
        
        The energies at the beginning and at the end of the ramp
        have to be specified in the settings property derived from
        the configuration file. Furthermore, the delta energy between
        two energy steps (step height) has to be given aswell.
        """
        # TODO: Write this function
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
        #TODO: Write this function
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

    def determine_settletime(self)
        """Determine settle time of the electronics.
        """
        # TODO
        
    def timer_resovled_measurement(self)
        # TODO
        
    @abstractmethod
    def measure_now(self):
        """Take a measurement.

        This method must be reimplemented in subclasses. It is
        supposed to be called after preparing the controller for
        a measurement and after an energy has been set on the 
        controller. It should only emit a send_message signal 
        which triggers a measurement and sends additional data 
        if needed.
            
        It should take all required data for this operation from
        the settings property derived from the configuration file.

        Returns
        -------
        None.
        """
        return

    @abstractmethod
    def receive_measurements(self)
        """Receives measurements from the serial.
        
        This function has to be reimplemented in the subclasses.
        Upon receiving the data_received signal this function is 
        supposed to process and return measurements. All of the
        settings required for processing different measurements
        should be derived from the configuration file.
        (i.e.: Different measurements may require other conversion
        factors.)
        
        Returns
        -------
        measurements : array of floats
            Received measurements
        """
        return