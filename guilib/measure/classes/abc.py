"""Module abc of viperleed.guilib.measure.classes.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-04-04
Author: Michele Riva
Author: Florian Doerr

This module contains QObjectSettingsErrors, QObjectWithError,
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

from abc import ABCMeta
from abc import abstractmethod
from configparser import ConfigParser
from dataclasses import dataclass, field
from operator import itemgetter
from pathlib import Path

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure.hardwarebase import DEFAULTS_PATH
from viperleed.guilib.measure.hardwarebase import ViPErLEEDErrorEnum
from viperleed.guilib.measure.hardwarebase import emit_error
from viperleed.guilib.measure.classes.settings import NoDefaultSettingsError
from viperleed.guilib.measure.classes.settings import NoSettingsError
from viperleed.guilib.measure.classes.settings import SettingsError
from viperleed.guilib.measure.classes.settings import (
    TooManyDefaultSettingsError
    )
from viperleed.guilib.measure.classes.settings import ViPErLEEDSettings
from viperleed.guilib.measure.dialogs.settingsdialog import SettingsHandler


class QMetaABC(ABCMeta, type(qtc.QObject)):
    """Metaclass common to QObject and ABCMeta allowing @abstractmethod."""


class QObjectSettingsErrors(ViPErLEEDErrorEnum):
    """Class for settings errors of ViPErLEED objects."""

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
                                  'Default settings corrupted. {!r}')

class QObjectWithError(qtc.QObject):                                            # TODO: The Measure class was meant to inherit from this class. Due to double inheritance from QObject this is not possible through standard inheritance.
    """Base class of measurement objects with error detection."""

    # Emitted whenever an error has been detected. Contains
    # information about the error in the form (code, message).
    error_occurred = qtc.pyqtSignal(tuple)


class QObjectWithSettingsABC(QObjectWithError, metaclass=QMetaABC):
    """Abstract base class of measurement objects with settings."""

    _mandatory_settings = []

    def __init__(self, *args, settings=None, **kwargs):
        """Initialise instance."""
        self._settings = ViPErLEEDSettings()
        self._settings_to_load = settings
        self._uses_default_settings = False
        if not settings:
            # Look for a default only if no settings are given.
            _name = type(self).__name__
            self._settings_to_load = self.find_default_settings(_name)
            self._uses_default_settings = True
        super().__init__(*args, **kwargs)

    @property
    def settings(self):
        """Return the current settings."""
        return self._settings

    @settings.setter
    def settings(self, new_settings):
        """Set new settings for this instance."""
        self.set_settings(new_settings)

    def find_default_settings(self, find_from, match_exactly=False):
        """Find default settings for this object.

        Parameters
        ----------
        find_from : SettingsInfo or str
            find_from contains information to look for in the
            configuration files. The way to determine the correct
            settings is up to the reimplementation of
            find_matching_settings_files in self. If it is
            a str it is converted into a SettingsInfo.
        match_exactly : bool, optional
            Whether find_from should be matched exactly.
            Default is False.

        Returns
        -------
        path_to_config : Path
            The path to the only settings file successfully found
            or the best match if no exact match was required.

        Raises
        -----
        NoDefaultSettingsError
            If no default settings were found.
        TooManyDefaultSettingsError
            If multiple matching default settings were found
            and an exact match was asked for.
        """
        # Make a dummy SettingsInfo that will
        # be used to find settings from default name.
        if isinstance(find_from, str):
            find_from = SettingsInfo(find_from)
            find_from.more['name'] = find_from.unique_name
        settings = self.find_matching_settings_files(
            find_from, DEFAULTS_PATH, match_exactly, default=True
            )
        if not settings:
            # No default settings was found.
            raise NoDefaultSettingsError(
                'No default settings found for '
                f'instance of {type(self).__name__}.'
                )
        if match_exactly and len(settings) != 1:
            # Too many default settings were found
            # while trying to find an exact match.
            raise TooManyDefaultSettingsError(
                'Too many default settings that match '
                f'instance of {type(self).__name__}.'
                )
        return settings[0]

    @classmethod
    def find_matching_settings_files(cls, obj_info, directory, match_exactly,
                                     default):
        """Find .ini files for obj_info in the tree starting at directory.

        Parameters
        ----------
        obj_info : SettingsInfo
            The additional information that should be used to find
            appropriate settings.
        directory : str or Path
            The location in which to look for configuration files.
        match_exactly : bool
            Whether obj_info should be matched exactly.
        default : bool
            Whether a default settings is searched or not. If True,
            the matching check for a default settings is performed.

        Returns
        -------
        obj_settings_files : list
            A list of the found settings paths that contain appropriate
            settings sorted by how well the settings match from best to
            worst.
        """
        settings_files = Path(directory).resolve().glob('**/*.ini')
        if not default:
            # Filter out default settings.
            settings_files = [file for file in settings_files
                              if not '_defaults' in str(file)]

        files_and_scores = []
        is_matching = (cls.is_matching_default_settings if default
               else cls.is_matching_settings)
        for settings_file in settings_files:
            config = ConfigParser()
            config.read(settings_file)
            conformity = is_matching(obj_info, config, match_exactly)
            if conformity:
                files_and_scores.append((settings_file, conformity))

        if not files_and_scores:
            return []
        obj_settings_files, _ = zip(
            *sorted(files_and_scores, key=itemgetter(1), reverse=True)
            )
        return obj_settings_files

    @classmethod
    @abstractmethod
    def is_matching_default_settings(cls, obj_info, config, match_exactly):
        """Determine if the default settings file is for this instance.

        Parameters
        ----------
        obj_info : SettingsInfo
            The information that should be used to check 'config'.
        config : ConfigParser
            The settings to check.
        match_exactly : bool
            Whether obj_info should be matched exactly.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used the sort the detected settings.
            Larger values in the tuple indicate a higher degree of
            conformity. The order of the items in the tuple is the
            order of their significance. This return value is used
            to determine the best-matching settings files when
            multiple files are found.
        """

    @classmethod
    @abstractmethod
    def is_matching_settings(cls, obj_info, config, match_exactly):
        """Determine if the settings file is for this instance.

        Parameters
        ----------
        obj_info : SettingsInfo
            The information that should be used to check 'config'.
        config : ConfigParser
            The settings to check.
        match_exactly : bool
            Whether obj_info should be matched exactly.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used the sort the detected settings.
            Larger values in the tuple indicate a higher degree of
            conformity. The order of the items in the tuple is the
            order of their significance. This return value is used
            to determine the best-matching settings files when
            multiple files are found.
        """

    def are_settings_invalid(self, new_settings):
        """Check if there are any invalid settings.

        Reimplementations can add additional mandatory settings at
        runtime. The base implementation will check if all of the
        _mandatory_settings are present in the provided settings.

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

    @abstractmethod
    def get_settings_handler(self):
        """Return a SettingsHandler object for displaying settings.

        This method should be extended in subclasses, i.e., do
        handler = super().get_settings_handler(), and then add
        appropriate sections and/or options to it using the
        handler.add_section, and handler.add_option methods.

        The base-class implementation returns a handler that
        contains the location of the settings file.

        Returns
        -------
        handler : SettingsHandler
            The handler used in a SettingsDialog to display the
            settings of this instance to users.
        """
        if self._uses_default_settings:
            raise SettingsError(
                f'Instance of {type(self).__name__} tried to return '
                'a settings handler using default settings.'
                )
        if not self.settings:
            # Remember to catch this exception before catching
            # SettingsError. Otherwise NoSettingsError will be
            # swallowed up as it is a subclass of SettingsError.
            raise NoSettingsError(
                f'Instance of {type(self).__name__} tried to return '
                'a settings handler without settings.'
                )
        handler = SettingsHandler(self.settings)
        handler.add_static_section('File')
        widget = qtw.QLabel()
        file = self.settings.last_file
        widget.setText(file.stem if file else 'None')
        handler.add_static_option(
            'File', 'config', widget,
            display_name='Settings file',
            tooltip=str(file) if file else None,
            )
        return handler

    def set_settings(self, new_settings):
        """Set new settings for this instance.

        Check and set settings. This method is used in
        the setter of the settings property.

        Parameters
        ----------
        new_settings : dict or ConfigParser or str or Path
            The new settings.

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
        QObjectSettingsErrors.MISSING_SETTINGS
            If the settings have not been found
            or if new_settings was None.
        QObjectSettingsErrors.INVALID_SETTINGS
            If any element of the new_settings does
            not fit the mandatory settings.
        """
        try:                                                                    # TODO: make method that searches through invalid for old values and replaces deprecated ones, make it a method of the ViPErLEEDSettings class
            new_settings = ViPErLEEDSettings.from_settings(new_settings)
        except(ValueError, NoSettingsError):
            emit_error(self, QObjectSettingsErrors.MISSING_SETTINGS,
                       type(self).__name__)
            return False
        invalid = self.are_settings_invalid(new_settings)
        if invalid:
            emit_error(self, QObjectSettingsErrors.INVALID_SETTINGS,
                       type(self).__name__, ', '.join(invalid), '')
            return False

        self._settings = new_settings
        return True


# DISABLE: seems a bug in pylint. This is also an ABC.
# pylint: disable-next=abstract-method
class HardwareABC(QObjectWithSettingsABC):
    """Abstract base class of hardware related objects."""

    # Emitted whenever the busy state of the device changes.
    # Contains the new busy state of the device.
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
        """Set busy state of instance.

        Note that 'busy =' will be performed in the calling
        thread. To avoid this, one can invoke the set_busy slot.
        """
        self.set_busy(is_busy)

    def get_busy(self):
        """Return whether the instance is busy.

        Note that this method is called in the getter of the busy
        property. If this method is overridden, this will also
        change the value that is returned by the busy property.

        Returns
        -------
        busy : bool
            True if the instance is busy.
        """
        return self._busy

    @qtc.pyqtSlot(bool)
    def set_busy(self, is_busy):
        """Set busy state of the instance.

        Parameters
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
        if was_busy != is_busy:
            self._busy = is_busy
            self.busy_changed.emit(self.busy)


class DeviceABC(HardwareABC):
    """Abstract base class of hardware device objects."""

    @abstractmethod
    def list_devices(self):
        """List all devices of this class.

        This method must return a list of SettingsInfo instances. Each
        controller is represented by a single SettingsInfo instance. The
        SettingsInfo object must contain a .unique_name and can contain
        .more information as a dict. The information contained within
        a SettingsInfo must be enough to determine a suitable settings
        file for the device from it.

        Returns
        -------
        devices : list
            Each element is a SettingsInfo instance containing the name
            of a device and additional information as a dict.
        """


@dataclass
class SettingsInfo:
    """A container for information about settings of objects.

    Instances of this class should be constructed in a manner that
    allows the constructing class to find settings from the instance.

    Attributes
    ----------
    unique_name : str
        Unique name identifying the discovered device
    more : dict
        Extra, optional, information about the discovered device.
    """
    unique_name: str
    more: dict = field(default_factory=dict)

    def __post_init__(self):
        """Check that we have a string unique_name."""
        if not isinstance(self.unique_name, str):
            raise TypeError(f'{type(self).__name__}: '
                            'unique_name must be a string')
