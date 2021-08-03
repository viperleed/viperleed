"""Module controllerabc of viperleed.?????.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-08
Author: Michele Riva
Author: Florian Doerr

This module contains the definition of the ControllerABC abstract
base class, used for giving basic commands to the LEED electronics.
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
        ]
    controller_busy = qtc.pyqtSignal(bool)

    def __init__(self, settings=None, port_name='', sets_energy=False):
        """Initialize the controller instance."""

        # Switches the controller from set_energy() to measure_now().
        self.__sets_energy = sets_energy
        self.__port_name = port_name
        self.__settings = None
        self.__serial = None

        self.set_settings(settings)

        # Is used to determine if the next step
        # in the measurement cycle can be done.
        self.__busy = False

    @property
    def busy(self):
        """Return whether the controller is busy."""
        return self.__busy

    @busy.setter
    def busy(self, is_busy):
        """Set the controller to busy True/False.

        Parameters
        ----------
        is_busy : bool
            True if the controller is busy
        """
        was_busy = self.__busy
        is_busy = bool(is_busy)
        if was_busy is not is_busy:
            self.__busy = is_busy
            self.controller_busy.emit(self.busy)

    @property
    def serial(self):
        """Return the serial port instance used."""
        return self.__serial

    def __get_sets_energy(self):
        """Return whether the controller sets the energy."""
        return self.__sets_energy

    def set_sets_energy(self, energy_setter):
        """Set the serial to controls energy True/False.

        Parameters
        ----------
        energy_setter : bool
            True if the controller sets the energy.
        """
        self.__sets_energy = bool(energy_setter)

    sets_energy = property(__get_sets_energy, set_sets_energy)

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
        """
        if new_settings is None:
            emit_error(self, ControllerErrors.MISSING_SETTINGS)
            return

        new_settings, invalid = config_has_sections_and_options(
            self, new_settings,
            self._mandatory_settings
            )
        for setting in invalid:
            emit_error(self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                       setting)

        if invalid:
            return

        serial_cls_name = new_settings.get('controller', 'serial_port_class')
        if self.serial.__class__.__name__ != serial_cls_name:
            try:
                serial_class = class_from_name('serial', serial_cls_name)
            except ValueError:
                emit_error(self, ControllerErrors.INVALID_CONTROLLER_SETTINGS,
                           setting)
                return
            self.__serial = serial_class(settings=new_settings,
                                         port_name=self.__port_name)

        # The next line will also check that new_settings contains
        # appropriate settings for the serial port class used
        self.serial.port_settings = new_settings
        self.__settings = self.serial.port_settings

    settings = property(__get_settings, set_settings)

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

    def true_energy_to_setpoint(self, energy):
        """Take requested energy and convert it to the energy to set.

        The conversion is done by reading a polynomial from
        the config files which is a function of the true energy
        and yields the energy to set.

        Parameters
        ----------
        energy : float
            Requested energy in eV

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
