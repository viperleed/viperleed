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
from abc import abstractmethod

from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure.controller.controllerabc import ControllerABC
from viperleed.guilib.measure import hardwarebase


class MeasureController(ControllerABC):
    """Controller class for measurement controllers."""

    # Signal which is used to forward data and let the measurementabc
    # class know that the controller is done measuring.
    data_ready = qtc.pyqtSignal(object)

    # This signal is only used by the primary controller which
    # sets the energy.
    about_to_trigger = qtc.pyqtSignal()

    def __init__(self, settings=None, port_name='', sets_energy=False):
        """Initialise controller class object."""

        super().init(settings=settings, port_name=port_name,
                     sets_energy=sets_energy)

        # This dictionary must be reimplemented in subclasses.
        # It must contain all possible measurement types the controller
        # can receive using the same keys used in the measurement
        # class responsible for receiving data from this controller.
        self.__measurements = defaultdict(list)

        # These dictionaries must be reimplemented in subclasses.
        # They must contain all functions the MeasureController has to call
        # in the order they have to be called to bring the controller into
        # a state ready for measurements.
        # prepare_to_dos_1 contains everything that has to be done before
        # the starting energy has been set.
        # prepare_to_dos_2 contains everything that has to be done after
        # the starting energy has been set.
        self.prepare_to_dos_1 = defaultdict(bool)
        self.prepare_to_dos_2 = defaultdict(bool)

        # Connect data_received signal from the serial to
        # the receive_measurements function in this class.
        self.serial.data_received.connect(self.receive_measurements)

        # Connect serial about_to_trigger signal to controller
        # about_to_trigger signal.
        self.serial.about_to_trigger.connect(self.about_to_trigger.emit)

        # tuple used to store the energies sent
        # by the measurementABC class.
        self.__energies = []
        # tuple used to store the times sent
        # by the measurementABC class.
        self.__settle_times = []

    def handle_do_measurement(self, energy, *other_data, **kwargs):
        """Handle the do measurement command.

        Decide what to do depending on the settings of
        the controller.

        Returns
        -------
        None.
        """
        # TODO: move to subclass?
        self.busy = True
        # TODO: other parameters not documented and not used
        if self.sets_energy:
            energies_and_times = [energy, *other_data]
            self.set_energy(self.true_energy_to_setpoint(energy))
        # TODO: fix energy setpoint function
        else:
            # TODO: delay by time from configuration
            self.measure_now()

    @abstractmethod
    def abort(self):
        """Abort current task.

        This method must be reimplemented in subclasses.
        Abort what the controller is doing right now and
        return to waiting for further instructions.

        Returns
        -------
        None.
        """
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

        If the controller already automatically takes a measurement
        after setting an energy it can be a no op.

        It should take all required data for this operation from
        the settings property derived from the configuration file.

        Returns
        -------
        None.
        """
        return

    def begin_preparation(self, busy):
        """Prepare the controller for a measurement.

        The prepare_to_dos_1 dictionary used in this method
        must be reimplemented in subclasses. The
        reimplementation should call functions that take the
        settings property derived from the configuration file
        and use it to do all required tasks before a measurement.
        (i.e. calibrating the electronics, selecting channels,
        determining the gain, ...)

        It should be able to select the update rate of the
        measurement electronics and change channels if there
        are more than one.

        Parameters
        ----------
        busy : boolean
            Busy state of the serial.
            If not busy, send next command.

        Returns
        -------
        None.
        """
        if not busy:
            next_to_do = None
            for key, to_do in self.prepare_to_dos_1.items():
                if not to_do:
                    continue
                next_to_do = key
                break

            if next_to_do:
                self.prepare_to_dos_1[next_to_do] = False
                if next_to_do == 'self.set_energy':
                    next_to_do(self.true_energy_to_setpoint(
                                self.__energies))
                    # TODO: fix this ^
                else:
                    next_to_do()
                return
            self.serial.serial_busy.disconnect()
            self.busy = False

    # @qtc.pyqtSlot(bool)
    def continue_preparation(self, busy):
        """Prepare the controller for a measurement.

        The prepare_to_dos_2 dictionary used in this method
        must be reimplemented in subclasses. The
        reimplementation should call functions that take the
        settings property derived from the configuration file
        and use it to do all required tasks before a measurement.
        (i.e. calibrating the electronics, selecting channels,
        determining the gain, ...)

        It should be able to select the update rate of the
        measurement electronics and change channels if there
        are more than one.

        Parameters
        ----------
        busy : boolean
            Busy state of the serial.
            If not busy, send next command.

        Returns
        -------
        None.
        """
        if not busy:
            next_to_do = None
            for key, to_do in self.prepare_to_dos_2.items():
                if not to_do:
                    continue
                next_to_do = key
                break

            if next_to_do:
                self.prepare_to_dos_2[next_to_do] = False
                next_to_do()
                return
            self.serial.serial_busy.disconnect()
            self.busy = False
            self.data_ready.emit({})

    @abstractmethod
    def receive_measurements(self, receive):
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

        Parameters
        ----------
        receive : object
            Data received from the serial, most
            likely an array/a list of floats.

        Returns
        -------
        None.
        """

        self.ready()
        return

    def ready(self):
        """Check if all measurements have been received.

        Emit a signal which contains all of the measurements as soon
        as the data from the controller has been received.

        The busy attribute will let the measurement class know if
        it can continue with the next step in´the measurement cycle.
        Once all of the controllers and cameras are not busy anymore,
        the signal for the next step will be sent.

        Returns
        -------
        None.
        """
        self.busy = False
        self.data_ready.emit(self.__measurements)
        for key in self.__measurements:
            self.__measurements[key] = []

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
        it should call the measure_now function.

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

    def trigger_begin_preparation(self, energy, time):
        """Trigger the first step in the preparation for measurements.

        Set self.busy to true, reset all prepare_to_dos_1 and
        start first step of the preparation.

        Parameters
        ----------
        energy : float
            Starting energy the controller will set if
            sets_energy is true.
        time : int
            Settle time after the setting the starting
            energy

        Returns
        -------
        None.
        """
        self.busy = True
        self.__energies = energy
        self.__settle_times = time
        for key in self.prepare_to_dos_1:
            self.prepare_to_dos_1[key] = True
        self.serial.serial_busy.connect(self.begin_preparation,
                                        type=qtc.Qt.UniqueConnection)
        self.begin_preparation(False)

    def trigger_continue_preparation(self):
        """Trigger the second step in the preparation for measurements.

        Set self.busy to true, reset all prepare_to_dos_2 and
        start second step of the preparation.

        Returns
        -------
        None.
        """
        self.busy = True
        for key in self.prepare_to_dos_2:
            self.prepare_to_dos_2[key] = True
        self.serial.serial_busy.connect(self.continue_preparation,
                                        type=qtc.Qt.UniqueConnection)
        self.continue_preparation(False)
