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

from abc import abstractmethod
from collections import defaultdict
from copy import deepcopy

from numpy.polynomial.polynomial import Polynomial
from PyQt5 import QtCore as qtc

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.datapoints import QuantityInfo
from viperleed.guilib.measure.classes.settings import (
    ViPErLEEDSettings, NoDefaultSettingsError,
    NoSettingsError, NotASequenceError,
    )


class ControllerErrors(base.ViPErLEEDErrorEnum):
    """Class for controller errors."""
    # The following three are fatal errors, and should make the GUI
    # essentially unusable, apart from allowing to load appropriate
    # settings.
    INVALID_SETTINGS = (100,
                        "Invalid controller settings: Required "
                        "settings {} missing or values "
                        "inappropriate. Check configuration file.\n{}")
    MISSING_SETTINGS = (101,
                        "Controller cannot operate without settings. "
                        "Load an appropriate settings file before "
                        "proceeding.")
    DEFAULT_SETTINGS_CORRUPTED = (102,
                                  "No or multiple default settings "
                                  "found for controller class {!r}.")
    CANNOT_MEASURE = (103,
                      "A subclass of ControllerABC is not supposed to "
                      "measure any quantities. Subclass MeasureControllerABC "
                      "instead.")
    INVALID_SETTING_WITH_FALLBACK = (
        104,
        "Invalid/unreadable controller settings value {} for setting {!r}. "
        "Using {} instead. Consider fixing your configuration file."
        )


# too-many-public-methods
# too-many-instance-attributes
class ControllerABC(qtc.QObject, metaclass=base.QMetaABC):
    """Base class for giving orders to the LEED electronics."""

    error_occurred = qtc.pyqtSignal(tuple)

    # This signal is only used by the primary controller which
    # sets the energy. If the primary controller does not take
    # measurements then this signal needs to be emitted after
    # the primary controller has set the energy. See also the
    # data_ready signal.
    about_to_trigger = qtc.pyqtSignal()

    # Signal which is used to forward data and let the MeasurementABC
    # class know that the controller is done measuring. If the
    # primary controller does not take measurements it should emit an
    # empty dictionary.
    data_ready = qtc.pyqtSignal(dict)

    controller_busy = qtc.pyqtSignal(bool)

    _mandatory_settings = [
        ('controller', 'serial_port_class'),
        ]

    def __init__(self, settings=None, port_name='', sets_energy=False):
        """Initialise the controller instance.

        Parameters
        ----------
        settings : ViPErLEEDSettings or str or path-like, optional
            The controller settings. If not given, it should be set
            via the .settings property before the controller can
            be used.
        port_name : str, optional
            Name of the serial port to be used to communicate with
            the controller. If this is given, it will also be stored
            in the settings file, overriding the value that may be
            there. If not given and no value is present in the
            "controller/port_name" field, a port name should be set
            explicitly via the .serial.port_name property. Default
            is an empty string.
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
        self.__settings = ViPErLEEDSettings()
        self.__serial = None
        self.__hash = -1

        self.__init_errors = []  # Report these with a little delay
        self.__init_err_timer = qtc.QTimer(self)
        self.__init_err_timer.setSingleShot(True)
        self.__init_err_timer.timeout.connect(self.__report_init_errors)

        self.__port_name = port_name
        self.__energy_calibration = None

        # Set in self.set_energy to the sum of the waiting times
        self._time_to_trigger = 0

        # __can_continue_preparation decides whether the controller
        # can start the second part of the preparation, i.e., the
        # one that contains steps to be done after the LEED energy
        # has be already set.
        self.__can_continue_preparation = False

        self.error_occurred.connect(self.__on_init_errors)

        self.set_settings(settings)

        # Is used to determine if the next step
        # in the measurement cycle can be done.
        self.__busy = False

        # These dictionaries must be reimplemented in subclasses.
        # They must contain all functions the MeasureController has
        # to call in the order to bring the controller into a state
        # ready for setting the energy/taking measurements.
        # begin_prepare_todos contains everything that has to be
        # done before the starting energy has been set.
        # continue_prepare_todos contains everything that has to be
        # done after the starting energy has been set.
        self.begin_prepare_todos = defaultdict(bool)
        self.continue_prepare_todos = defaultdict(bool)

        # Sequence used to store the energies and times (in
        # alternating order) to be set during the preparation.
        self.first_energies_and_times = []

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
        # a tuple whose first element is the command and its associated
        # data, and with second element the timeout parameter.
        self.__unsent_messages = []
        if self.serial:
            self.serial.serial_busy.connect(self.send_unsent_messages,
                                            type=qtc.Qt.UniqueConnection)
        if self.__init_errors:
            self.__init_err_timer.start(20)
        self.error_occurred.disconnect(self.__on_init_errors)

    def __hash__(self):
        """Return modified hash of self."""
        if self.__hash == -1:
            self.__hash = hash((id(self), self.name))
        return self.__hash

    def __get_busy(self):  # pylint: disable=unused-private-member
        """Return whether the controller is busy."""
        if self.serial and self.serial.busy:
            return True
        return self.__busy

    def set_busy(self, is_busy):
        """Set the controller to busy True/False.

        No change of the busy state can occur if the controller
        has a queue of messages that have still to be sent to
        the hardware. Also, the value given is discarded if the
        serial port is currently busy (i.e., waiting for a response).
        In this latter case, the controller will always be busy.

        Parameters
        ----------
        is_busy : bool
            True if the controller is busy.

        Emits
        -----
        controller_busy
            If the busy state of the controller changed.
        """
        if self.__unsent_messages:
            return
        was_busy = self.busy
        is_busy = bool(is_busy)
        if self.serial and self.serial.busy:
            is_busy = True
        self.__busy = is_busy
        if was_busy is not is_busy:
            self.controller_busy.emit(self.busy)

    busy = property(__get_busy, set_busy)

    @property
    def energy_calibration_curve(self):
        """Return a callable for converting energies.

        Returns
        -------
        energy_calibration_curve : callable
            When called with an array-like argument, interprets
            the argument as nominal electron energies and returns
            the corresponding setpoint energies for the LEED optics.
        """
        if self.__energy_calibration is None:
            try:
                coef = self.settings.getsequence('energy_calibration',
                                                 'coefficients',
                                                 fallback=(0, 1))
            except NotASequenceError:
                coef = (0, 1)
                base.emit_error(
                    self, ControllerErrors.INVALID_SETTING_WITH_FALLBACK,
                    '', 'energy_calibration/coefficients', coef
                    )
            try:
                domain = self.settings.getsequence('energy_calibration',
                                                   'domain',
                                                   fallback=(-10, 1100))
            except NotASequenceError:
                # pylint: disable=redefined-variable-type
                # Likely pylint bug? domain is Sequence also before.
                # No reasonable need to bother the user with this. We
                # anyway construct the polynomial from scratch.
                domain = (-10, 1100)

            self.__energy_calibration = Polynomial(coef, domain=domain,
                                                   window=domain)
        return self.__energy_calibration

    @property
    def hv_settle_time(self):
        """Return the time it takes to settle the energy (msec)."""
        if not self.sets_energy:
            raise RuntimeError("Should never ask the hv_settle_time for "
                               "a controller that does not .sets_energy")
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug
        fallback = 2000
        if not self.settings:
            return fallback
        try:
            # Mandatory for a controller that .sets_energy.
            settle_t = self.settings.getint('measurement_settings',
                                            'hv_settle_time')
        except (TypeError, ValueError):
            # Not an int
            settle_t = fallback
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            'measurement_settings/hv_settle_time')
        return settle_t

    @property
    def i0_settle_time(self):
        """Return the time it takes to settle the beam current (msec)."""
        if not self.sets_energy:
            raise RuntimeError("Should never ask the i0_settle_time for "
                               "a controller that does not .sets_energy")
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug
        fallback = 2000
        if not self.settings:
            return fallback
        try:
            # Mandatory for a controller that .sets_energy
            settle_t = self.settings.getint('measurement_settings',
                                            'i0_settle_time')
        except (TypeError, ValueError):
            # Not an int
            settle_t = fallback
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            'measurement_settings/i0_settle_time')
        return settle_t

    @property
    def initial_delay(self):
        """Return the initial time delay of a measurement (msec).

        Returns
        -------
        initial_delay : float
            The time interval between when a measurement was triggered
            and when the measurement was actually acquired. If mutliple
            measurements are averaged over, the "time when measurements
            are actually acqiuired" should be the middle time between
            the beginning and the end of the measurement.
        """
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug
        fallback = 0.0
        if not self.settings:
            return fallback
        try:
            delay = self.settings.getfloat('controller', 'initial_delay',
                                           fallback=fallback)
        except (TypeError, ValueError):
            # Not a float
            delay = fallback
            base.emit_error(self,
                            ControllerErrors.INVALID_SETTING_WITH_FALLBACK,
                            '', 'controller/initial_delay', delay)
        return delay

    @property
    def long_settle_time(self):
        """Return the first settle time (msec)."""                              # TODO: doc
        if not self.sets_energy:
            raise RuntimeError("Should never ask the long_settle_time for "
                               "a controller that does not .sets_energy")
        # pylint: disable=redefined-variable-type
        # Seems a pylint bug
        fallback = 5000
        if not self.settings:
            return fallback
        try:
            # Mandatory for a controller that .sets_energy
            settle_t = self.settings.getint('measurement_settings',
                                            'first_settle_time')
        except (TypeError, ValueError):
            # Not an int
            settle_t = fallback
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            'measurement_settings/i0_settle_time')
        return settle_t

    @property
    @abstractmethod
    def name(self):
        """Return a unique name for this controller."""
        return ""

    @property
    def measured_quantities(self):
        """Return measured quantities.

        This property can only be set by calling .set_measurements.

        Returns
        -------
        None.
        """
        # In ControllerABC, this method always returns an
        # empty list, as __measured_quantities cannot be set
        # otherwise.
        return tuple()

    @property
    def port_name(self):
        """Return the name of the serial port for this controller."""
        name = ""
        try:
            name = self.serial.port_name
        except AttributeError:
            pass
        return name if name else 'UNKNOWN PORT'

    @port_name.setter
    def port_name(self, port_name):
        """Set the port name for this controller."""
        if not isinstance(port_name, str):
            raise TypeError("port_name must be a string")
        self.__port_name = port_name
        if self.settings:
            self.settings.set('controller', 'port_name', port_name)
        if not self.serial:
            return
        self.serial.port_name = port_name
        if port_name:
            self.serial.serial_connect()

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

    @property
    def settings(self):
        """Return the current settings used as a ConfigParser."""
        return self.__settings

    @settings.setter
    def settings(self, new_settings):
        """Set new settings for this controller."""
        self.set_settings(new_settings)

    @property
    def time_to_trigger(self):
        """Return the time required for triggering measurements.

        Returns
        -------
        time_to_trigger : int
            Time interval in milliseconds between setting the
            last energy and emitting the .about_to_trigger
            signal. If this controller is not responsible for
            setting energies, time_to_trigger is zero.
        """
        if not self.sets_energy:
            return 0
        return self._time_to_trigger

    def set_settings(self, new_settings):  # too-complex
        """Set new settings for this controller.

        Settings are accepted and loaded only if they are valid.
        Notice that the serial port name in the new settings will
        be used only the first time a valid settings is loaded if
        no valid port_name was given before (either at instantiation
        or using the .port_name property). Explicitly use the
        .port_name property to set a new communication port. If
        you want the port name to be taken from new settings, set
        .port_name to an empty string, then load the new settings.

        Parameters
        ----------
        new_settings : dict or ConfigParser or None
            The new settings. If False-y, an attempt to retrieve
            settings from a default file will be made. new_settings
            (or the default) will be checked for the following
            mandatory sections/options:
                'controller'/'serial_port_class'
            (if self.sets_energy:
                'measurement_settings'/'i0_settle_time'
                'measurement_settings'/'hv_settle_time')

        Emits
        -----
        error_occurred
            If the new_settings are None, or if they do not match
            up with the mandatory settings required by the controller
            or if the serial class specified in the settings could not
            be instantiated.
        """
        # Look for a default only if no settings are given
        _name = self.__class__.__name__ if not new_settings else None
        try:
            new_settings = ViPErLEEDSettings.from_settings(new_settings,
                                                           find_from=_name)
        except NoDefaultSettingsError:
            base.emit_error(self, ControllerErrors.DEFAULT_SETTINGS_CORRUPTED,
                            _name)
            return
        except NoSettingsError:
            base.emit_error(self, ControllerErrors.MISSING_SETTINGS)
            return

        if not self.are_settings_ok(new_settings):
            return

        if not self.__port_name:
            self.__port_name = new_settings.get('controller', 'port_name',
                                                fallback=self.__port_name)
        else:
            new_settings['controller']['port_name'] = self.__port_name

        serial_cls_name = new_settings.get('controller', 'serial_port_class')
        if self.serial.__class__.__name__ != serial_cls_name:
            try:
                serial_class = base.class_from_name('serial', serial_cls_name)
            except ValueError:
                base.emit_error(
                    self, ControllerErrors.INVALID_SETTINGS,
                    'controller/serial_port_class', ''
                    )
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
        # settings are invalid (i.e., missing mandatory fields for
        # the serial)!
        if self.__port_name:
            self.serial.serial_connect()
        self.__settings = self.serial.port_settings
        self._time_to_trigger = 0
        self.__hash = -1

    # pylint: disable=no-self-use
    # Method is here for consistency with MeasureControllerABC
    def measures(self, _):
        """Return whether this controller measures quantity."""
        return False
    # pylint: enable=no-self-use

    @abstractmethod
    def are_settings_ok(self, settings):
        """Return whether a ViPErLEEDSettings is compatible with self.

        This method should be extended in subclasses, i.e., it
        should return False if super().are_settings_ok(settings)
        returns False. This method is guaranteed to be called
        once each time new settings are loaded, either via
        .set_settings() or via the .settings property. It is
        always passed a ViPErLEEDSettings instance. Settings will
        be loaded successfully only if this method returns True.
        The base implementation automatically checks for the
        presence of all _mandatory_settings. Thus, reimplementations
        may simply expand _mandatory_settings, then call super().

        Reimplementations can emit self.error_occured to inform
        users of invalid settings.

        Parameters
        ----------
        settings : ViPErLEEDSettings
            The settings to be checked.

        Returns
        -------
        settings_ok : bool
            True if settings are ok.
        """
        # The next extra settings are mandatory only for a
        # controller that sets the LEED energy on the optics
        extra_mandatory = ()
        if self.sets_energy:
            extra_mandatory = (('measurement_settings', 'i0_settle_time'),
                               ('measurement_settings', 'hv_settle_time'),
                               ('measurement_settings', 'first_settle_time'))

        invalid = settings.has_settings(*self._mandatory_settings,
                                        *extra_mandatory)

        if invalid:
            base.emit_error(self, ControllerErrors.INVALID_SETTINGS,
                            ', '.join(invalid), '')
            return False
        return True

    @abstractmethod
    def set_energy(self, energy, settle_time, *more_steps, trigger_meas=True):
        """Set electron energy on LEED controller.

        This method must be extended in subclasses by calling super()
        ant some point. The reimplementation should take the energy
        value(s) in eV, the wating times (in msec) and the keyword
        argument trigger_meas, together with other optional data to be
        read from self.settings, and turn them into messages that can
        be sent via self.send_message(*messages).

        Conversion from the desired, true electron energy (i.e.,
        the energy that the electrons will have when exiting the
        gun) to the setpoint value to be sent to the controller
        can be done inside this method by calling
        self.true_energy_to_setpoint(energy).

        Parameters
        ----------
        energy : float
            Nominal electron energy in electronvolts.
        settle_time : integer
            Interval in milliseconds that the controller will wait
            before deeming the LEED optics stable at set energy.
        *more_steps : Number
            If given, it should be an even number of elements.
            Odd elements are energies, even ones settle-time intervals.
            Multiple steps can be executed quickly after each other.
            The last step will be the final energy that is set and
            should ensure stabilization of the electronics.
        trigger_meas : bool, optional
            True if the controllers are supposed to take measurements
            after the energy has been set. Default is True. The base
            ControllerABC cannot take measurements, therefore this
            parameter is only used to emit an about to trigger if True.

        Returns
        -------
        None.
        """
        self._time_to_trigger = settle_time + sum(more_steps[1::2])

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
            Energy to set in eV in order to get requested energy.
        """
        to_setpoint = self.energy_calibration_curve
        return to_setpoint(energy)

    def send_message(self, *data, **kwargs):
        """Use serial to send message.

        Send message via serial. Save it if serial is busy
        or if there are other messages that have been stored.

        Parameters
        ----------
        *data : object
            Any data the controller class may need to send.
            Should have the same types and number of elements
            as self.serial.send_message.
        **kwargs : dict
            Keyword arguments passed on to self.serial.send_message.

        Returns
        -------
        None.
        """
        if self.__unsent_messages or self.serial.busy:
            self.__unsent_messages.append((data, kwargs))
            return
        self.serial.send_message(*data, **kwargs)

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
            data, kwargs = self.__unsent_messages.pop(0)
            self.serial.send_message(*data, **kwargs)

    def flush(self):
        """Clear unsent messages and set not busy."""
        self.__unsent_messages = []
        self.busy = False

    @abstractmethod
    def abort_and_reset(self):  # TODO: Confusing. Perhaps remove altogether?
        """Abort current task and reset the controller.

        This method must be extended in subclasses by calling
        super(). Abort what the controller is doing right now,
        reset it and return to waiting for further instructions.

        Returns
        -------
        None.
        """
        self._time_to_trigger = 0
        # TODO: check other stuff that has to be reset.

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

    def set_measurements(self, quantities):
        """Set measured_quantities property."""
        if quantities:
            base.emit_error(self, ControllerErrors.CANNOT_MEASURE)

    @abstractmethod
    def list_devices(self):
        """List all devices of this class."""
        return

    def connect_(self):
        """Connect serial port."""
        # TODO: should we complain if .__port_name is False-y?
        if not self.serial or self.serial.is_open:
            # Invalid or already connected
            return
        self.serial.serial_connect()

    def disconnect_(self):
        """Diconnect serial port."""
        try:
            self.serial.serial_disconnect()
        except (TypeError, AttributeError):
            pass

    def __do_preparation_step(self, serial_busy=False):
        """Prepare the controller for a measurement cycle.

        Runs, in key-insertion order, over the self.begin_prepare_todos
        and self.continue_prepare_todos dictionaries. The two dicts
        should be populated by subclasses. The former contains tasks
        used during the first 'step' of the preparation: all those
        that can be done with the LEED energy set to zero. The latter
        all the tasks that require a nonzero LEED energy.

        Parameters
        ----------
        serial_busy : bool
            Busy state of the serial. If not busy, send next command.

        Returns
        -------
        None.
        """
        if serial_busy:
            return

        todos = self.begin_prepare_todos
        if self.__can_continue_preparation:
            todos = self.continue_prepare_todos

        next_to_do = None
        for method_name, to_be_done in todos.items():
            if not to_be_done:
                continue
            next_to_do = getattr(self, method_name)
            break
        if next_to_do:
            todos[next_to_do.__name__] = False
            if next_to_do.__name__ == 'set_energy':
                # Never trigger measurements during preparation
                next_to_do(*self.first_energies_and_times, trigger_meas=False)
            else:
                next_to_do()
            return

        # Disconnect signal: will be reconnected
        # during the call to continue_preparation
        self.serial.serial_busy.disconnect(self.__do_preparation_step)
        if self.__can_continue_preparation:
            # The whole preparation is now over.
            # Guarantee that any unsent message that may have come
            # while the preparation was running is now sent.
            base.safe_connect(self.serial.serial_busy,
                              self.send_unsent_messages,
                              type=qtc.Qt.UniqueConnection)
        self.busy = False

    def begin_preparation(self, energies_and_times):
        """Trigger the first step in the preparation for measurements.

        Set self.busy to True, reset all begin_prepare_todos and start
        first step of the preparation. The .controller_busy() signal
        will be emitted carrying False once all steps are complete.

        Parameters
        ----------
        energies_and_times : tuple
            Starting sequence of energies and times the controller
            will use to set the energy during preparation. These
            quantities are used only if self.sets_energy is True.

        Returns
        -------
        None.
        """
        self.reset_preparation_todos()
        self.first_energies_and_times = energies_and_times

        # Clear any unsent messages. Any message that comes during
        # preparation will not be sent till the end of the whole
        # preparation. This excludes a call to .stop(), which is
        # always done as soon as possible.
        self.__unsent_messages = []
        self.__can_continue_preparation = False

        self.busy = True
        base.safe_disconnect(self.serial.serial_busy,
                             self.send_unsent_messages)
        base.safe_connect(self.serial.serial_busy, self.__do_preparation_step,
                          type=qtc.Qt.UniqueConnection)
        self.__do_preparation_step()

    def continue_preparation(self):
        """Trigger the second step in the preparation for measurements.

        Set self.busy to true, reset all continue_prepare_todos
        and start second step of the preparation.

        Returns
        -------
        None.
        """
        self.__can_continue_preparation = True

        self.busy = True
        base.safe_disconnect(self.serial.serial_busy,
                             self.send_unsent_messages)
        base.safe_connect(self.serial.serial_busy, self.__do_preparation_step,
                          type=qtc.Qt.UniqueConnection)
        self.__do_preparation_step()

    def reset_preparation_todos(self):
        """Reset all the segments of preparation.

        Segments are in the form [key: bool}. Those segments for
        which the value is True will be done one at a time, in
        the order in which keys are inserted. Segments are stored
        in self.begin_prepare_todos and self.continue_prepare_todos.

        This method can be reimplemented in subclasses if some segments
        are to be skipped (e.g., depending on the firmware version of
        the hardware).

        The base-class implementation sets all segments to True.

        Returns
        -------
        None.
        """
        for key in self.begin_prepare_todos:
            self.begin_prepare_todos[key] = True
        for key in self.continue_prepare_todos:
            self.continue_prepare_todos[key] = True


class MeasureControllerABC(ControllerABC):
    """Controller class for measurement controllers."""

    def __init__(self, settings=None, port_name='', sets_energy=False):
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

        super().__init__(settings=settings, port_name=port_name,
                         sets_energy=sets_energy)

        # This dictionary must be reimplemented in subclasses.
        # It must contain all possible measurement types the controller
        # can receive using the same keys used in the measurement
        # class responsible for receiving data from this controller.
        self.measurements = defaultdict(list)

        # This tuple contains the QuantityInfo Enums associated with
        # the selected quantities to measure.
        self.__measured_quantities = tuple()

    @abstractmethod
    def measure_now(self):
        """Take a measurement.

        This method must be extended in subclasses by calling
        super(). It is supposed to be called after preparing
        the controller for a measurement and after an energy
        has been set on the controller. It should only trigger
        measurements (and send additional data if needed).

        It should take all required data for this operation from
        self.settings (and other attributes of self).

        Returns
        -------
        None.
        """
        self._time_to_trigger = 0

    @property
    def measured_quantities(self):
        """Return a list of measured quantities."""
        return self.__measured_quantities

    @property
    @abstractmethod
    def measurement_interval(self):
        """Return the time interval between measurements (msec)."""
        return

    def measures(self, quantity):
        """Return whether this controller measures quantity."""
        return quantity in self.measured_quantities

    @abstractmethod
    def on_data_ready(self, data):
        """Receive and store data from the serial.

        This method is the slot connected to the data_received
        signal of self.serial and has to be reimplemented in
        subclasses.

        Measurements received from the serial should be processed,
        appropriately converted to the correct physical units
        and stored into self.measurements, a dictionary whose
        keys are QuantityInfo objects.

        When all the measurements expected has been received,
        this method should call self.measurements_done().

        Notice that data may also be non-measurement information.
        This should also be processed and stored, if needed. In
        this case, do not call .measurements_done().

        All of the settings required for processing data should
        be derived from self.settings. This includes, for example,
        conversion factors/curves.

        Parameters
        ----------
        data : object
            Data received from self.serial.

        Returns
        -------
        None.
        """
        self.measurements_done()

    def measurements_done(self):
        """Emit measurements and change busy mode.

        The busy attribute will let the measurement class know if
        it can continue with the next step in the measurement cycle.
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
        self.data_ready.emit(deepcopy(self.measurements))
        for key in self.measurements:
            self.measurements[key] = []

    @abstractmethod
    def set_energy(self, energy, settle_time, *more_steps, trigger_meas=True):
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
        settle_time : integer
            Interval in milliseconds that the controller will wait
            before deeming the LEED optics stable at set energy.
        *more_steps : Number
            If given, it should be an even number of elements.
            Odd elements are energies, even ones settle-time intervals.
            Multiple steps can be executed quickly after each other.
            The last step will be the final energy that is set and
            should ensure stabilization of the electronics.
        trigger_meas : bool, optional
            True if the controllers are supposed to take
            measurements after the energy has been set.
            Default is True.

        Returns
        -------
        None.
        """
        super().set_energy(energy, settle_time, *more_steps,
                           trigger_meas=trigger_meas)

    @abstractmethod
    def set_measurements(self, quantities):
        """Decide what to measure.

        This method must be reimplemented in subclasses. It
        should take requested measurement types as strings
        from the MeasurementABC class (i.e.: I0, I_Sample, ...),
        check if those types are available and not conflicting
        with each other and decide which channels to use.

        super().set_measurements(requested) must be called in
        subclasses at the end of the reimplementation to actually
        store the measured quantities.

        Parameters
        ----------
        quantities : list of strings
            Contains all of the requested
            measurement types.

        Returns
        -------
        None.
        """
        self.__measured_quantities = tuple(QuantityInfo.from_label(q)
                                           for q in quantities)

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
            # the on_data_ready function in this class.
            self.serial.data_received.connect(self.on_data_ready)

            # Connect serial about_to_trigger signal to controller
            # about_to_trigger signal.
            self.serial.about_to_trigger.connect(self.about_to_trigger)

    @abstractmethod
    def set_continuous_mode(self, continuous=True):
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
        continuous : bool, optional
            Wether continuous mode should be on.
            Used in subclass. Default is True.

        Returns
        -------
        None.
        """
        base.safe_connect(self.serial.serial_busy, self.set_busy,
                          type=qtc.Qt.UniqueConnection)
