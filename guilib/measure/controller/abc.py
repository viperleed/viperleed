"""Module abc of viperleed.guilib.measure.controller.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-08
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the ControllerABC and the
MeasureController class abstract base classes and their associated
ViPErLEEDErrorEnum class ControllerErrors used for giving basic
commands to the LEED electronics.
"""

# Python standard modules
import ast
from abc import abstractmethod
from collections import defaultdict
from configparser import ConfigParser

from numpy.polynomial.polynomial import Polynomial
from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure.hardwarebase import (
    config_has_sections_and_options, class_from_name, emit_error, QMetaABC,
    ViPErLEEDErrorEnum
    )


class ControllerErrors(ViPErLEEDErrorEnum):
    # The following two are fatal errors, and should make the GUI
    # essentially unusable, apart from loading appropriate settings.
    INVALID_CONTROLLER_SETTINGS = (100,
                                   "Invalid controller settings: Required "
                                   "settings {!r} missing or values "
                                   "inappropriate. Check configuration file.")
    MISSING_SETTINGS = (101,
                        "Controller cannot operate without settings. "
                        "Load an appropriate settings file before "
                        "proceeding.")


class ControllerABC(qtc.QObject, metaclass=QMetaABC):
    """Base class for giving orders to the LEED electronics."""

    error_occurred = qtc.pyqtSignal(tuple)
    
    # This signal is only used by the primary controller which
    # sets the energy. If the primary controller does not take
    # measurements then this signal needs to be emitted after
    # the primary controller has set the energy as well as an
    # empty data_ready signal.
    about_to_trigger = qtc.pyqtSignal()
    
    # Signal which is used to forward data and let the MeasurementABC
    # class know that the controller is done measuring. If the
    # primary controller does not take measurements it should emit an
    # empty dictionary.
    data_ready = qtc.pyqtSignal(object, dict)

    _mandatory_settings = [
        ('controller', 'serial_port_class'),
        ('available_commands',)
        ]
    controller_busy = qtc.pyqtSignal(bool)

    def __init__(self, settings, port_name='', sets_energy=False):
        """Initialise the controller instance.

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

        super().__init__()
        self.__sets_energy = sets_energy
        self.__settings = None
        self.__serial = None
        self.__hash = -1
        self.measured_quantities = []

        if not port_name:
            if not settings.has_option('controller', 'port_name'):
                raise TypeError("No port name given, and none found in the "
                                "configuration file. Cannot instantiate "
                                "a controller without a valid port.")
            port_name = settings.get('controller', 'port_name',
                                     fallback='None')
        else:
            settings.set('controller', 'port_name', port_name)
        self.__port_name = port_name

        self.__init_errors = []  # Report these with a little delay
        self.__init_err_timer = qtc.QTimer(self)
        self.__init_err_timer.setSingleShot(True)

        self.error_occurred.connect(self.__on_init_errors)
        self.__init_err_timer.timeout.connect(self.__report_init_errors)

        self.set_settings(settings)

        # Is used to determine if the next step
        # in the measurement cycle can be done.
        self.__busy = False

        # __unsent_messages is a list of messages that have been
        # stored by the controller because it was not yet possible
        # to send them to the hardware controller. New messages
        # will automatically be appended to __unsent_messages if
        # the serial is still waiting for a response from the
        # hardware. All messages in __unsent_messages are processed
        # in the order they arrive at. When the serial receives
        # a valid answer from the hardware it will automatically
        # trigger an attempt to send the next message in
        # __unsent_messages. Each element of unsent_messages is
        # a tuple containing the command and its associated data
        # and the timeout parameter.
        self.__unsent_messages = []
        self.serial.serial_busy.connect(self.send_unsent_messages,
                                        type=qtc.Qt.UniqueConnection)
        if self.__init_errors:
            self.__init_err_timer.start(20)
        self.error_occurred.disconnect(self.__on_init_errors)

    def __hash__(self):
        """Return modified hash of self."""
        if self.__hash == -1:
            self.__hash = hash((id(self), self.__port_name))
        return self.__hash

    def __get_busy(self):
        """Return whether the controller is busy."""
        if self.serial.busy:
            return True
        return self.__busy

    def set_busy(self, is_busy):
        """Set the controller to busy True/False.

        Parameters
        ----------
        is_busy : bool
            True if the controller is busy

        Emits
        -----
        controller_busy
            If the busy state of the controller changed.
        """
        if self.__unsent_messages:
            return
        # was_busy = self.__busy
        was_busy = self.busy
        self.__busy = bool(is_busy) if not self.serial.busy else True
        is_busy = self.busy
        if was_busy is not is_busy:
            # self.__busy = is_busy if not self.serial.busy else True
            self.controller_busy.emit(self.busy)

    busy = property(__get_busy, set_busy)

    @property
    def initial_delay(self):
        """Return the initial time delay of a measurement in seconds."""
        return self.settings.getfloat(
            'controller', 'initial_delay', fallback=0
            )

    @property
    def serial(self):
        """Return the serial port instance used."""
        return self.__serial

    @property
    def sets_energy(self):
        """Return whether the controller sets the energy."""
        return self.__sets_energy

    @sets_energy.setter
    def sets_energy(self, energy_setter):
        """Set the serial to controls energy True/False.

        Parameters
        ----------
        energy_setter : bool
            True if the controller sets the energy.
        """
        self.__sets_energy = bool(energy_setter)

    def __get_settings(self):
        """Return the current settings used as a ConfigParser."""
        return self.__settings

    def set_settings(self, new_settings):
        """Set new settings for this controller.

        Settings are accepted and loaded only if they are valid.

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
        if new_settings is None:
            emit_error(self, ControllerErrors.MISSING_SETTINGS)
            return

        # The next extra setting is mandatory only for a controller
        # that sets the LEED energy on the optics
        extra_mandatory_list = []
        if self.sets_energy:
            extra_mandatory_list = [('measurement_settings', 'i0_settle_time'),
                                    ('measurement_settings', 'hv_settle_time')]

        new_settings, invalid = config_has_sections_and_options(
            self, new_settings,
            (*self._mandatory_settings, *extra_mandatory_list)
            )

        if invalid:
            error_msg = ', '.join(invalid)
            emit_error(self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                       error_msg)
            return

        serial_cls_name = new_settings.get('controller', 'serial_port_class')
        if self.serial.__class__.__name__ != serial_cls_name:
            try:
                serial_class = class_from_name('serial', serial_cls_name)
            except ValueError:
                emit_error(self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                           'controller/serial_port_class')
                return
            self.__serial = serial_class(new_settings,
                                         port_name=self.__port_name,
                                         parent=self)
            self.serial.error_occurred.connect(self.error_occurred)
        else:
            # The next line will also check that new_settings contains
            # appropriate settings for the serial port class used.
            self.serial.port_settings = new_settings
            self.serial.port_name = self.__port_name

        # Notice that the .connect() will run anyway, even if the
        # settings are invalid (i.e., missing mandatory fields)!
        self.serial.serial_connect()
        self.__settings = self.serial.port_settings
        self.__hash = -1

    settings = property(__get_settings, set_settings)

    def measures(self, quantity):
        """Return whether this controller measures quantity."""
        return quantity in self.measured_quantities

    @abstractmethod
    def set_energy(self, energy, *other_data, trigger_meas=True, **kwargs):
        """Set electron energy on LEED controller.

        This method must be reimplemented in subclasses. The
        reimplementation should take the energy value in eV
        and other optional data needed by the serial interface
        and turn them into a message that can be sent via
        self.__serial.send_message(message, *other_messages).

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
            self.true_energy_to_setpoint(energy).
        trigger_meas : bool, optional
            True if the controllers are supposed to take
            measurements after the energy has been set.
            Default is True. The base ControllerABC cannot
            take measurements, therefore this parameter
            is only used to emit an about to trigger
            if True.
        **kwargs : object
            Unused keyword arguments.

        Returns
        -------
        None.
        """
        return

    def true_energy_to_setpoint(self, energy):
        """Take requested energy and convert it to the energy to set.

        The conversion is done by reading a polynomial from
        the configuration files which is a function of the true
        energy and yields the energy to set. This polynomial is
        determined in the energy calibration and written to the
        config. This function should always be called when
        setting an energy.

        Parameters
        ----------
        energy : float
            Requested energy in eV.

        Returns
        -------
        new_energy : float
            Energy to set in eV in order
            to get requested energy.
        """
        calibration_coef = ast.literal_eval(
            self.settings['energy_calibration']['coefficients']
            )
        calibration_domain = ast.literal_eval(
            self.settings['energy_calibration']['domain']
            )
        calibration = Polynomial(calibration_coef, domain=calibration_domain,
                                 window=calibration_domain)
        new_energy = calibration(energy)

        return new_energy

    def make_serial(self):
        """Create serial port object and connect to it."""
        self.serial.port = self.settings.get('controller', 'port_name',
                                             fallback='None')

    def send_message(self, *data):
        """Use serial to send message.

        Send message via serial. Save it if serial is busy
        or if there are other messages that have been stored.

        Parameters
        ----------
        data : tuple
            Any data the controller class may need to send.

        Returns
        -------
        None.
        """
        if self.__unsent_messages or self.serial.busy:
            self.__unsent_messages.append(data)
            return
        self.serial.send_message(*data)

    def send_unsent_messages(self, serial_busy):
        """Send messages that have been stored.

        Parameters
        ----------
        serial_busy : boolean
            Busy state of the serial.
            If not busy, send next unsent message.

        Returns
        -------
        None.
        """
        if serial_busy:
            return
        if self.__unsent_messages:
            data = self.__unsent_messages.pop(0)
            self.serial.send_message(*data)

    def flush(self):
        """Clear unsent, set busy false and send abort"""
        self.__unsent_messages = []
        self.busy = False

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
    def stop(self):
        """Stop.

        Stop whatever the controller is doing right now
        and return to idle state.

        Returns
        -------
        None.
        """
        try:
            self.serial.serial_busy.disconnect()
        except TypeError:
            pass
        self.serial.serial_busy.connect(self.send_unsent_messages,
                                        type=qtc.Qt.UniqueConnection)
        self.serial.serial_busy.connect(self.set_busy,
                                        type=qtc.Qt.UniqueConnection)

    def __on_init_errors(self, err):
        """Collect initialization errors to report later."""
        self.__init_errors.append(err)

    def __report_init_errors(self):
        """Emit error_occurred for each initialization error."""
        for error in self.__init_errors:
            self.error_occurred.emit(error)
        self.__init_errors = []

    def what_to_measure(self, *args, **kwargs):
        """Set measured_quantities property."""
        self.measured_quantities = []


class MeasureControllerABC(ControllerABC):
    """Controller class for measurement controllers."""

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
        controller. It should only trigger a measurement and
        send additional data if needed.

        It should take all required data for this operation from
        the settings property derived from the configuration file.

        Returns
        -------
        None.
        """
        return

    @property
    @abstractmethod
    def measurement_interval(self):
        """Return the time interval between measurements in seconds."""
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
        self.serial.serial_busy.connect(self.send_unsent_messages,
                                        type=qtc.Qt.UniqueConnection)
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
        self.data_ready.emit(self, self.measurements.copy())
        for key in self.measurements:
            self.measurements[key] = []

    @abstractmethod
    def set_energy(self, energy, *other_data, trigger_meas=True, **kwargs):
        """Set electron energy on LEED controller.

        This method must be reimplemented in subclasses. The
        reimplementation should take the energy value in eV
        and other optional data needed by the serial interface
        and turn them into a message that can be sent via
        self.send_message(message, *other_messages).

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
        trigger_meas : bool, optional
            True if the controllers are supposed to take
            measurements after the energy has been set.
            Default is True.
        **kwargs : object
            Other unused keyword arguments.

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
        from the MeasurementABC class (i.e.: I0, I_Sample, ...),
        check if those types are available and not conflicting
        with each other and decide which channels to use.
            
        super().what_to_measure(requested) must be called in
        subclasses at the end of the reimplementation.

        Parameters
        ----------
        requested : list of strings
            Contains all of the requested
            measurement types.

        Returns
        -------
        None.
        """
        self.measured_quantities = requested

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

        Has to be reimplemented in subclasses. If continuous is
        true the controller has to continue measuring and return
        data without receiving further instructions. The serial_busy
        has to be hooked up to the busy state of the controller.
        Call super() to enable switching of busy state of controller
        once the continuous mode has been set on the hardware
        controller.

        Parameters
        ----------
        continuous : bool
            Wether continuous mode should be on.
            Used in subclass.

        Returns
        -------
        None.
        """
        try:
            self.serial.serial_busy.connect(self.set_busy,
                                            type=qtc.Qt.UniqueConnection)
        except TypeError:
            pass
