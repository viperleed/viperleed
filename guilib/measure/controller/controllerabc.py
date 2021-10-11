"""Module controllerabc of viperleed.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-08
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the ControllerABC abstract
base class and its associated ViPErLEEDErrorEnum class ControllerErrors
used for giving basic commands to the LEED electronics.
"""

# Python standard modules
import ast
from abc import abstractmethod
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

    _mandatory_settings = [
        ('controller', 'serial_port_class'),
        ('available_commands',)
        ]
    controller_busy = qtc.pyqtSignal(bool)

    def __init__(self, settings, port_name='', sets_energy=False):                # NOTE: we will need a settings file that remembers which devices were connected to which port ------ Do we?
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

        if not port_name:
            if not settings.has_option('controller', 'port_name'):
                raise TypeError("No port name given, and none found in the "
                                "configuration file. Cannot instantiate "
                                "a controller without a valid port.")
            port_name = settings.get('controller', 'port_name')
        else:
            settings.set('controller', 'port_name', port_name)
        self.__port_name = port_name

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
        is_busy = bool(is_busy)
        if was_busy is not is_busy:
            self.__busy = is_busy if not self.serial.busy else True
            self.controller_busy.emit(self.busy)

    busy = property(__get_busy, set_busy)

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
        extra_mandatory = ('measurement_settings', 'settle_time')
        if (self.sets_energy
                and extra_mandatory not in self._mandatory_settings):
            self._mandatory_settings.append(extra_mandatory)

        new_settings, invalid = config_has_sections_and_options(
            self, new_settings,
            self._mandatory_settings
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
        else:
            # The next line will also check that new_settings contains
            # appropriate settings for the serial port class used.
            self.serial.port_settings = new_settings
            self.serial.port_name = self.__port_name

        # Notice that the .connect() will run anyway, even if the
        # settings are invalid (i.e., missing mandatory fields)!
        self.serial.serial_connect()
        self.__settings = self.serial.port_settings

    settings = property(__get_settings, set_settings)

    @abstractmethod
    def set_energy(self, energy, *other_data):
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
        self.serial.port = self.settings.get('controller', 'port_name')

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
        print(self.serial.port_name, *data)
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
            data = self.__unsent_messages.pop()
            print(self.serial.port_name, *data)
            self.serial.send_message(*data)

    def flush(self):
        """Clear unsent, set busy false and send abort"""
        self.__unsent_messages = []
        self.busy = False
        self.abort_and_reset()

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
