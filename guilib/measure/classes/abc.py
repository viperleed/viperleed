"""Module abc of viperleed.guilib.measure.classes.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-04-04
Author: Michele Riva
Author: Florian Doerr

This module contains QObjectABCErrors, QObjectWithError,
QObjectWithSettingsABC, HardwareABC and DeviceABC.
QObjectWithError is the base class of QObjectWithSettingsABC,
CalibrationTask, and DataPoints. QObjectWithSettingsABC is the base
class of HardwareABC and MeasurementABC. HardwareABC is the base class
of DeviceABC and SerialABC. DeviceABC is the base class of ControllerABC
and CameraABC.

Note that all of the classes defined in this module are subclasses of
QObject. Any subclass of a QObject cannot inherit from a second QObject,
which means one cannot do: NewClass(QObjectSubclass1, QObjectSubclass2)
"""

from abc import abstractmethod

from PyQt5 import QtCore as qtc

from viperleed.guilib.measure.hardwarebase import emit_error
from viperleed.guilib.measure.hardwarebase import QMetaABC
from viperleed.guilib.measure.hardwarebase import ViPErLEEDErrorEnum
from viperleed.guilib.measure.classes.settings import NoDefaultSettingsError
from viperleed.guilib.measure.classes.settings import NoSettingsError
from viperleed.guilib.measure.classes.settings import ViPErLEEDSettings


class QObjectABCErrors(ViPErLEEDErrorEnum):                                # TODO: make sure all classes attach their class name if they use these errors.
    """Class for errors shared among many ViPErLEED objects."""

    MISSING_SETTINGS = (900,
                        '{} cannot operate without settings. Load an '
                        'appropriate settings file before proceeding.'
                        )
    INVALID_SETTINGS = (901,
                        'Invalid settings for instance of {}. Required '
                        'settings {} missing or value inappropriate. '
                        'Check configuration file.\n{}'
                        )
    INVALID_SETTING_WITH_FALLBACK = (
        902,
        'Invalid settings for instance of {}. Invalid/unreadable '
        'settings value {} for setting {!r}. Using {} instead. '
        'Consider fixing your configuration file.'
        )
    DEFAULT_SETTINGS_CORRUPTED = (903,
                                  'No or multiple default settings '
                                  'found for instance of {!r}.')

class QObjectWithError(qtc.QObject):                                            # TODO: The Measure class was meant to inherit from this class. Due to double inheritance from QObject this is not possible through standard inheritance.
    """Base class of measurement objects with error detection."""

    # Emitted whenver an error has been detected. Contains
    # information about the error which occurred.
    error_occurred = qtc.pyqtSignal(tuple)


class QObjectWithSettingsABC(QObjectWithError, metaclass=QMetaABC):
    """Abstract base class of measurement objects with settings."""

    _mandatory_settings = []

    def __init__(self, *args, **kwargs):
        """Initialise instance."""
        self._settings = ViPErLEEDSettings()
        super().__init__(*args, **kwargs)

    @property
    def settings(self):
        """Return the current settings."""
        return self._settings

    @settings.setter
    def settings(self, new_settings):
        """Set new settings for this instance."""
        self.set_settings(new_settings)

    # @abstractmethod
    # def find_settings(self):                                                    # TODO: abstract or base implementation?
        # """Find appropriate settings for this instance."""
        # # How about making hardwarebase._find_matching_configs a class
        # # specific method that is capable to deal with outdated
        # # settings files instead of just looking for the exact name in
        # # the settings file?

    def are_settings_invalid(self, new_settings):
        """Check if there are any invalid settings.

        Reimplementations can add additional mandatory settings in
        run time.

        Parameters
        ----------
        new_settings : ViPErLEEDSettings
            The new settings.

        Returns
        -------
        invalid_settings : list
            Invalid required_settings of self as a list of strings.
            Each entry can be either '<section>', '<section>/<option>',
            or '<section>/<option> not one of <value1>, <value2>, ...'
        """
        return new_settings.has_settings(*self._mandatory_settings)

    def set_settings(self, new_settings, find_from=None):                                       # TODO: make sure reimplementations and extensions also return bool
        """Set new settings for this instance.

        Parameters
        ----------
        new_settings : dict or ConfigParser or str or Path
            The new settings.
        find_from : str or None, optional
            The string to look for in a configuration file to
            be loaded and returned. If None, no search will be
            performed.

        Returns
        -------
        settings_valid : bool
            True if the new settings given were accepted.

        Raises
        ------
        TypeError
            If new_settings is neither a dict, ConfigParser, str,
            or Path
        TypeError
            If an element of the mandatory settings is None or has
            a length greater than 3.

        Emits
        -----
        QObjectABCErrors.MISSING_SETTINGS
            If the settings have not been found
            or if new_settings was None.
        QObjectABCErrors.INVALID_SETTINGS
            If any element of the new_settings does
            not fit the mandatory settings.
        """
        try:                                                                    # TODO: make method that searches through invalid for old values and replaces deprecated ones
            new_settings = ViPErLEEDSettings.from_settings(new_settings,
                                                           find_from=find_from)
        except(ValueError, NoSettingsError):
            emit_error(self, QObjectABCErrors.MISSING_SETTINGS,
                       type(self).__name__)
            return False
        except NoDefaultSettingsError:
            base.emit_error(self, QObjectABCErrors.DEFAULT_SETTINGS_CORRUPTED,
                            find_from)
            return False
        invalid = self.are_settings_invalid(new_settings)
        if invalid:
            emit_error(self, QObjectABCErrors.INVALID_SETTINGS,
                       type(self).__name__, ', '.join(invalid), '')
            return False

        self._settings = new_settings
        return True


# DISABLE: seems a bug in pylint. This is also an ABC.
# pylint: disable-next=abstract-method
class HardwareABC(QObjectWithSettingsABC):
    """Abstract base class of hardware related objects."""

    # Emitted whenever the busy state of the device changes.
    # Contains the busy state the device changed to.
    busy_changed = qtc.pyqtSignal(bool)

    def __init__(self, *args, **kwargs):
        """Initialise instance."""
        # Busy state of the device.
        self._busy = False
        super().__init__(*args, **kwargs)

    @property
    def busy(self):
        """Return busy state of instance."""
        return self.get_busy()

    @busy.setter
    def busy(self, is_busy):
        """Set busy state of instance."""
        # Note that "busy =" will therefore be
        # performed in the calling thread.
        self.set_busy(is_busy)

    def get_busy(self):
        """Return whether the instance is busy.

        Returns
        -------
        busy : bool
        True if the instence is busy.
        """
        return self._busy

    @qtc.pyqtSlot(bool)
    def set_busy(self, is_busy):
        """Set busy state of the instance.

        Paramaters
        ----------
        is_busy : bool
            True if the instance is busy.

        Emits
        -----
        busy_changed
            If the busy state of the instance changed.
        """
        was_busy = self.busy
        is_busy = bool(is_busy)
        if was_busy is not is_busy:
            self._busy = is_busy
            self.busy_changed.emit(self.busy)


class DeviceABC(HardwareABC):
    """Abstract base class of hardware device objects."""

    # def __init__(self, *args, **kwargs):
        # """Initialise instance."""
        # super().__init__(*args, **kwargs)

    @abstractmethod
    def list_devices(self):                                                     # TODO: Base implementation?
        """List all devices of this class.

        Returns
        -------
        devices : list
            Each element is a string representing the name of a device.
        """
