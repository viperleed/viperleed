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

import ast
from collections import defaultdict
from abc import abstractmethod

from numpy.polynomial.polynomial import Polynomial
from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure.controller.controllerabc import ControllerABC
from viperleed.guilib.measure import hardwarebase

class MeasureControllerError(hardwarebase.ViPErLEEDErrorEnum):
    UNKNOWN_MEASUREMENT = (150, "The measurement to do is of an unkown type.")


class MeasureController(ControllerABC):
    """Controller class for the ViPErLEED Arduino Micro.

    The plot_info dictionary in this class may be reimplemented
    in subclasses. Each section (label) needs to contain
    the unit and the scaling ('lin' for linear and 'log' for
    logarithmic scaling.
    """
    # Is emitted if a measurement cycle has been completed.
    cycle_completed = qtc.pyqtSignal()
    do_next_measurement = qtc.pyqtSignal()

    # The reimplementation may introduce more/other keys.
    # See ViPErinoController for an example on how to do this.
    plot_info = defaultdict(list)
    plot_info['nominal_energies'] = ('eV', 'lin')
    plot_info['I0'] = ('uA', 'lin')
    plot_info['measured_energies'] = ('eV', 'lin')
    plot_info['elapsed_time'] = ('ms', 'lin')
    plot_info['aux0'] = ('V', 'log')

    def __init__(self, settings=None, port_name='', controls_camera=False):
        """Initialise controller class object."""

        super().init(settings=settings, port_name=port_name,
                     controls_camera=controls_camera)

        # Connect data_received signal from the serial to
        # the receive_measurements function in this class.
        self.serial.data_received.connect(self.receive_measurements)

        # The reimplementation may introduce more/other keys.
        self.data_points = defaultdict(list)
        for key in self.plot_info:
            self.data_points[key] = []

        # Flags needed for loop operation.
        self.flags = defaultdict(list)
        self.flags['image_received'] = False
        self.flags['measurements_received'] = False

        # Attributes needed for loop operation
        self.current_energy = 0
        self.cycle_type = ''
        self.end_energy = float(
            self.settings['measurement_settings']['end_energy']
            )


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

        If this function does not already trigger a measurement
        it should emit an about_to_trigger_signal which in turn
        calls the measure_now function.

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
        via the send_command signal, which takes instructions
        as a tuple and forwards them to the send_message
        function of the serialabc class.

        It should be able to select the update rate of the
        measurement electronics and change channels if there
        are more than one. It will also need to set the
        self.current_energy attribute to the starting energy
        of the cycle.

        Returns
        -------
        None.
        """
        return

    def measure_iv_video(self):
        """Measure a ramp of energies.

        This method must be reimplemented in subclasses. This
        function should take measurements while doing an energy
        ramp. It should measure I0 (the current emitted by the
        filament used to create the electron beam) and take
        pictures of the LEED screen for every energy step.

        The energies at the beginning and at the end of the ramp
        have to be specified in the settings property derived
        from the configuration file. Furthermore, the delta
        energy between two energy steps (step height) has to be
        given as well.

        Parameters
        ----------
        current_energy

        Returns
        -------
        None.
        """
        # Edit docstring
        # Do stuff here


        self.current_energy += self.end_energy
        if self.current_energy > self.end_energy:
            self.cycle_completed.emit()
        return

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

        self.current_energy += self.end_energy
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
        return

    def determine_settle_time(self):
        """Determine settle time of the electronics.

        Returns
        -------
        None.
        """
        # Edit docstring
        # TODO
        return

    def time_resolved_measurement(self):
        """Measure same energy over time.

        This function should take measurements at the same
        energy over a specified amount of time. All required
        data should be taken from the settings property derived
        from the configuration file.

        If continuous mode is implemented, the function should
        be able to determine if using continuous mode or triggered
        mode is the more appropriate choice.

        Returns
        -------
        None.
        """
        # Edit docstring
        # TODO
        return

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
    def receive_image(self, object):
        """Need to wait for camera class.
        
        Process image.
        
        After the image has been processed, the function has to 
        set the image_received flag to True and check if the 
        next measurement can start via the 
        ready_for_next_measurement function.
        
        Returns
        -------
        None.
        """
        # Edit docstring
        
        self.flags['image_received'] = True
        self.ready_for_next_measurement()
        
        return

    @abstractmethod
    def receive_measurements(self, object):
        """Receive measurements from the serial.

        This function has to be reimplemented in subclasses.
        Upon receiving the data_received signal this function
        is supposed to process and append measurements to the
        appropriate attribute of the class.

        All of the settings required for processing different
        measurements should be derived from the configuration
        file. (i.e.: Different measurements may require other
        conversion factors.) The channels selected in the
        settings can be used to determine to which section of
        the data_points library the measurement should be
        appended to.

        After the measurements have been processed, the
        function has to set the measurements_received flag to
        True and check if the next measurement can start via 
        the ready_for_next_measurement function.
        
        Returns
        -------
        None.
        """
        
        self.flags['measurements_received'] = True
        self.ready_for_next_measurement()
        
        return

    def ready_for_next_measurement(self):
        """Check if all measurements have been received.

        The next measurement will be started as soon as both
        the image from the camera and the data from the controller
        have been received. If the controller class does not control
        a camera, it will not wait for an image to be returned. If
        the camera is in live mode, the controller class will not
        wait for measurement data from the controller.

        Returns
        -------
        None.
        """

        if not self.__controls_camera:
            self.flags['image_received'] = True
        if self.settings['camera_settings']['mode'] == 'live':
            self.flags['measurements_received'] = True
        if self.image_received and self.measurements_received:
            # Rather call a function here
            # self.do_next_measurement.emit()
            self.flags['image_received'] = False
            self.flags['measurements_received'] = False
            self.do_next_measurement()
        
        return
        
    def do_next_measurement(self):
        """Determine measurement function.
        
        Decide which measurement to do, by checking the orders
        that were given by the gui. The only possible orders
        are "iv_video" for a LEED I(V) measurement, 
        "time_resolved" for a time resolved measurement,
        "energy_calibration" for determining the energy 
        calibration setpoint and "settle_time" for determining
        the settle time of the LEED electronics.
        
        Returns
        -------
        None.
        """
        if self.cycle_type == "iv_video":
            self.measure_iv_video()
        if self.cycle_type == "time_resolved":
            self.measure_energy_setpoint()
        if self.cycle_type == "energy_calibration":
            self.measure_energy_setpoint()
        if self.cycle_type == "settle_time":
            self.measure_energy_setpoint()
        else:
            self.error_occurred.emit(MeasureControllerError.UNKNOWN_MEASUREMENT)
         
        return
