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
from abc import ABCMeta, abstractmethod
from configparser import ConfigParser

from numpy.polynomial.polynomial import Polynomial
from PyQt5 import QtCore as qtc

# ViPErLEED modules
from viperleed.guilib.measure import get_serial
from viperleed.guilib.measure.hardwarebase import (
    config_has_sections_and_options,
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


class ControllerABC(metaclass=ABCMeta):
    """Base class for giving orders to the LEED electronics."""

    error_occurred = qtc.pyqtSignal(tuple)

    def __init__(self, settings=None, port_name='', controls_camera=False):
        """Initialise the controller instance."""
        self.__controls_camera = controls_camera
        mandatory_settings = (('controller', 'serial_port_class'),)

        new_settings, invalid = config_has_sections_and_options(
            self,
            new_settings,
            mandatory_settings
            )
        if invalid:
            if hasattr(invalid, '__len__'):
                for setting in invalid:
                    (error_code,
                     error_msg) = ControllerErrors.INVALID_CONTROLLER_SETTINGS
                     error_msg = error_msg.format(setting)
                     self.error_occurred.emit((error_code, error_msg))
            else:
                self.error_occurred.emit(ControllerErrors.MISSING_SETTINGS)

        self.__settings = new_settings
        serial_name = settings.get('controller', 'serial_port_class')
        self.__serial = get_serial(serial_name)(settings=self.settings,
                                                port_name=port_name)

    @property
    def controls_camera(self):
        """Return whether this controller also manage the camera."""
        return self.__controls_camera

    @controls_camera.setter
    def controls_camera(self, active):
        """Set whether this controller handles the camera."""
        self.__controls_camera = bool(active)

    @property
    def serial(self):
        """Return the serial port instance used."""
        return self.__serial

    @property
    def settings(self):
        """Return the current settings used as a ConfigParser."""
        return self.__settings

    @settings.setter
    def settings(self, new_settings):
        """Set new settings for this controller.

        Parameters
        ----------
        new_settings : dict or ConfigParser
            The new settings. Will be checked for the following
            mandatory sections/options:
        """
        mandatory_settings = (('controller', 'serial_port_class'),)

        new_settings, invalid = config_has_sections_and_options(
            self,
            new_settings,
            mandatory_settings
            )

        if invalid:
            if hasattr(invalid, '__len__'):
                for setting in invalid:
                    (error_code,
                     error_msg) = ControllerErrors.INVALID_CONTROLLER_SETTINGS
                     error_msg = error_msg.format(setting)
                     self.error_occurred.emit((error_code, error_msg))
            else:
                self.error_occurred.emit(ControllerErrors.MISSING_SETTINGS)

        # The next line will also check that new_settings contains
        # appropriate settings for the serial port class used
        self.serial.port_settings = new_settings
        self.__settings = new_settings

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
        energy : float
            Energy to set in eV in order
            to get requested energy
        """
        calibration_coef = ast.literal_eval(
            self.settings['energy_calibration']['coefficients']
            )
        calibration_domain = ast.literal_eval(
            self.settings['energy_calibration']['domain']
            )
        calibration = Polynomial(calibration_coef, domain=calibration_domain,
                                 window=calibration_domain)
        return calibration(energy)
