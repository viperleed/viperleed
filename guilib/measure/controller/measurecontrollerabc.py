"""Module measurecontrollerabc of viperleed.
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


class MeasureController(ControllerABC):
    """Controller class for measurement controllers."""

    # Signal which is used to forward data and let the MeasurementABC
    # class know that the controller is done measuring.
    data_ready = qtc.pyqtSignal(object)

    # This signal is only used by the primary controller which
    # sets the energy.
    about_to_trigger = qtc.pyqtSignal()

    _mandatory_settings = [*ControllerABC._mandatory_settings,
                           ('controller', 'measurement_devices')]

    def __init__(self, settings, port_name='', sets_energy=False):
        """Initialise controller class object.

        This is an upgraded version of its parent class as it
        instantiates multiple measurement related properties.

        Parameters
        ----------
        settings : ConfigParser
            The controller settings
        port_name : str, optional
            Name of the serial port to be used to communicate with
            the controller. This parameter is optional only in case
            settings contains a 'controller'/'port_name' option. If
            this is given, it will also be stored in the settings
            file, overriding the value that may be there. Default is
            an empty string.
        sets_energy : bool, optional
            Used to determine whether this controller is responsible
            for setting the electron energy by communicating with the
            LEED optics. Only one controller may be setting the energy.
            Default is False.

        Raises
        ------
        TypeError
            If no port_name is given, and none was present in the
            settings file.
        """

        super().__init__(settings, port_name=port_name, sets_energy=sets_energy)

        # This dictionary must be reimplemented in subclasses.
        # It must contain all possible measurement types the controller
        # can receive using the same keys used in the measurement
        # class responsible for receiving data from this controller.
        self.measurements = defaultdict(list)

        # These dictionaries must be reimplemented in subclasses.
        # They must contain all functions the MeasureController has
        # to call in the order to bring the controller into a state
        # ready for measurements. begin_prepare_todos contains
        # everything that has to be done before the starting energy
        # has been set. continue_prepare_todos contains everything
        # that has to be done after the starting energy has been set.
        self.begin_prepare_todos = defaultdict(bool)
        self.continue_prepare_todos = defaultdict(bool)

        # tuple used to store the energies and times sent
        # by the MeasurementABC class in alternating order.
        self.__energies_and_times = []

    @abstractmethod
    def abort_and_reset(self):
        """Abort current task and reset the controller.

        This method must be reimplemented in subclasses.
        Abort what the controller is doing right now, reset
        it and return to waiting for further instructions.

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

        It should take all required data for this operation from
        the settings property derived from the configuration file.

        Returns
        -------
        None.
        """
        return

    # @qtc.pyqtSlot(bool)
    def begin_preparation(self, serial_busy):
        """Prepare the controller for a measurement.

        The begin_prepare_todos dictionary used in this method
        must be reimplemented in subclasses. The
        reimplementation should call functions that take the
        settings property and use it to do all required tasks
        before a measurement. (i.e. calibrating the electronics,
        selecting channels, determining the gain, ...)

        It should be able to select the update rate of the
        measurement electronics and change channels if there
        are more than one.

        Parameters
        ----------
        serial_busy : boolean
            Busy state of the serial.
            If not busy, send next command.

        Returns
        -------
        None.
        """
        if serial_busy:
            return

        next_to_do = None
        for key, to_do in self.begin_prepare_todos.items():
            if not to_do:
                continue
            next_to_do = getattr(self, key)
            break
        if next_to_do:
            self.begin_prepare_todos[next_to_do.__name__] = False
            if next_to_do == self.set_energy:
                next_to_do(*self.__energies_and_times)
            else:
                next_to_do()
            return
        self.serial.serial_busy.disconnect()
        self.busy = False

    # @qtc.pyqtSlot(bool)
    def continue_preparation(self, serial_busy):
        """Prepare the controller for a measurement.

        The continue_prepare_todos dictionary used in this
        method must be reimplemented in subclasses. The
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
        serial_busy : boolean
            Busy state of the serial.
            If not busy, send next command.

        Returns
        -------
        None.
        """
        if serial_busy:
            return
        next_to_do = None
        for key, to_do in self.continue_prepare_todos.items():
            if not to_do:
                continue
            next_to_do = getattr(self, key)
            break

        if next_to_do:
            self.continue_prepare_todos[next_to_do.__name__] = False
            next_to_do()
            return

        self.serial.serial_busy.disconnect()
        self.busy = False

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

        Parameters
        ----------
        receive : object
            Data received from the serial, most
            likely an array/a list of floats.

        Returns
        -------
        None.
        """

        self.measurements_done()
        return

    def measurements_done(self):
        """Emit measurements and change busy mode.

        The busy attribute will let the measurement class know if
        it can continue with the next step in´the measurement cycle.
        Once all of the controllers and cameras are not busy anymore,
        the signal for the next step will be sent.

        Emits
        -----
        data_ready
            A signal containing the collected data, which
            triggers a check if all controllers and cameras
            are already done.

        Returns
        -------
        None.
        """
        self.busy = False
        self.data_ready.emit(self.measurements)
        for key in self.measurements:
            self.measurements[key] = []

    @abstractmethod
    def set_energy(self, energy, *other_data):
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
            self.true_energy_to_setpoint(energy).

        Returns
        -------
        None.
        """
        return

    def trigger_begin_preparation(self, energies_and_times):
        """Trigger the first step in the preparation for measurements.

        Set self.busy to true, reset all begin_prepare_todos
        and start first step of the preparation.
        energies_and_times is a tuple containing the energies
        and times to set during the preparation. First the energy
        should be set and afterwards the gain should be determined.

        Parameters
        ----------
        energies_and_times : tuple
            Starting energies and times the controller will
            use if sets_energy is true.

        Returns
        -------
        None.
        """
        self.busy = True
        self.__energies_and_times = energies_and_times
        for key in self.begin_prepare_todos:
            self.begin_prepare_todos[key] = True
        self.serial.serial_busy.connect(self.begin_preparation,
                                        type=qtc.Qt.UniqueConnection)
        self.begin_preparation(serial_busy=False)

    def trigger_continue_preparation(self):
        """Trigger the second step in the preparation for measurements.

        Set self.busy to true, reset all continue_prepare_todos
        and start second step of the preparation.

        Returns
        -------
        None.
        """
        self.busy = True
        for key in self.continue_prepare_todos:
            self.continue_prepare_todos[key] = True
        self.serial.serial_busy.connect(self.continue_preparation,
                                        type=qtc.Qt.UniqueConnection)
        self.continue_preparation(serial_busy=False)

    @abstractmethod
    def what_to_measure(self, requested):
        """Decide what to measure.

        This method must be reimplemented in subclasses. It
        should take requested measurement types as strings
        from the MeasurementABC class (i.e.: I0, Isample, ...),
        check if those types are available and not conflicting
        with each other and decide which channels to use.

        Additionally, class attributes should be implemented
        in subclasses which remember which measurements were
        requested in order to use them afterwards when creating
        dictionaries to return data.

        Parameters
        ----------
        requested : list of strings
            Contains all of the requested
            measurement types.

        Returns
        -------
        None.
        """
        return

    def set_settings(self, new_settings):
        """Set new settings for this controller.

        Settings are accepted and loaded only if they are valid.
        This is an upgraded version of the set_settings in the
        parent class as this function also connects the
        about_to_trigger and the data_received signals.

        Parameters
        ----------
        new_settings : dict or ConfigParser
            The new settings. Will be checked for the following
            mandatory sections/options:
                'controller'/'serial_port_class'

        Emits
        -----
        error_occurred
            If the new_settings are None, or if they do not match
            up with the mandatory settings required by the controller
            or if the serial class specified in the settings could not
            be instantiated.
        """
        super().set_settings(new_settings)

        if self.serial is not None:
            # Connect data_received signal from the serial to
            # the receive_measurements function in this class.
            self.serial.data_received.connect(self.receive_measurements)

            # Connect serial about_to_trigger signal to controller
            # about_to_trigger signal.
            self.serial.about_to_trigger.connect(self.about_to_trigger.emit)

    @abstractmethod
    def set_continuous_mode(self, continuous):
        """Set continuous mode.

        Has to be reimplemented in subclasses. If true the
        controller has to continue measuring and return data
        without receiving further instructions. The serial_busy
        has to be hooked up to the busy state of the controller.
        Call super() to enable switching of busy state of controller
        once the continuous mode has been set on the hardware
        controller.

        Parameters
        ----------
        continuous : bool
            True if continuous mode is on.

        Returns
        -------
        None.
        """
        if continuous:
            self.serial.serial_busy.connect(
                self.set_busy,
                type=qtc.Qt.UniqueConnection
                )
        else:
            try:
                self.serial.serial_busy.disconnect()
            except TypeError:
                # serial_busy has not been connected to anything
                pass
