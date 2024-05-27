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
from operator import itemgetter
from pathlib import Path

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure.hardwarebase import DEFAULTS_PATH
from viperleed.guilib.measure.hardwarebase import emit_error
from viperleed.guilib.measure.hardwarebase import SettingsInfo
from viperleed.guilib.measure.hardwarebase import ViPErLEEDErrorEnum
from viperleed.guilib.measure.classes.settings import NoSettingsError
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
    NO_DEFAULT_SETTINGS_FOUND = (903,
                                 'No default settings found '
                                 'for instance of {!r}.')
    TOO_MANY_DEFAULT_SETTINGS = (904,
                                 'Too many default settings that '
                                 'match instance of {!r} exactly.')

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

    def find_default_settings(self, find_from, exact_match=False):
        """Find default settings for this object.

        Parameters
        ----------
        find_from : SettingsInfo or str
            find_from contains information to look for in the
            configuration files. The way to determine the correct
            settings is up to the reimplementation of
            find_matching_settings_files in self. If it is
            a str it is converted into a SettingsInfo.
        exact_match : bool, optional
            Whether find_from should be matched exactly.
            Default is False.

        Returns
        -------
        path_to_config : Path or ''
            The path to the only settings file successfully found.
            '' if too many or no settings file were found.

        Emits
        -----
        QObjectSettingsErrors.NO_DEFAULT_SETTINGS_FOUND
            If no default settings were detected.
        QObjectSettingsErrors.TOO_MANY_DEFAULT_SETTINGS
            If too many default settings were detected.
        """
        # Make a dummy SettingsInfo that will
        # be used to find settings from default name.
        if isinstance(find_from, str):
            find_from = SettingsInfo(find_from)
            find_from.more['name'] = find_from.unique_name
        settings = self.find_matching_settings_files(
            find_from, DEFAULTS_PATH, exact_match, default=True
            )
        if not settings:
            # No default settings was found.
            emit_error(self, QObjectSettingsErrors.NO_DEFAULT_SETTINGS_FOUND,
                       find_from.unique_name)
            return ''
        if exact_match and len(settings) != 1:
            # Too many default settings were found
            # while trying to find an exact match.
            emit_error(self, QObjectSettingsErrors.TOO_MANY_DEFAULT_SETTINGS,
                       find_from.unique_name)
            return ''
        return settings[0]

    @classmethod
    def find_matching_settings_files(cls, obj_info, directory, exact_match,
                                     default):
        """Find .ini files for obj_info in the tree starting at directory.

        Parameters
        ----------
        obj_info : SettingsInfo
            The additional information that should be used to find
            appropriate settings.
        directory : str or Path
            The location in which to look for configuration files.
        exact_match : bool
            Whether obj_info should be matched exactly.
            If True, the information is matched exactly.
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
        directory = Path(directory).resolve()
        settings_files = directory.glob('**/*.ini')
        files_and_scores = []
        is_matching = (cls.is_matching_default_settings if default
               else cls.is_matching_settings)
        for settings_file in settings_files:
            config = ConfigParser()
            config.read(settings_file)
            conformity = is_matching(obj_info, config, exact_match)
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
    def is_matching_default_settings(cls, obj_info, config, exact_match):
        """Determine if the default settings file is for this instance.

        Parameters
        ----------
        obj_info : SettingsInfo
            The additional information that should
            be used to check settings.
        config : ConfigParser
            The settings to check.
        exact_match : bool
            Whether obj_info should be matched exactly.
            If True, the information is matched exactly.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used the sort the detected settings.
            Larger values in the tuple indicate a higher degree of
            conformity. The order of the values in the tuple is the
            order of their significance.
        """

    @classmethod
    @abstractmethod
    def is_matching_settings(cls, obj_info, config, exact_match):
        """Determine if the settings file is for this instance.

        Parameters
        ----------
        obj_info : SettingsInfo
            The additional information that should
            be used to check settings.
        config : ConfigParser
            The settings to check.
        exact_match : bool
            Whether obj_info should be matched exactly.
            If True, the information is matched exactly.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used the sort the detected settings.
            Larger values in the tuple indicate a higher degree of
            conformity. The order of the values in the tuple is the
            order of their significance.
        """

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
        handler = SettingsHandler(self.settings)
        handler.add_static_section('File')
        widget = qtw.QLabel()
        widget.setText(
            str(self.settings.last_file).rsplit('\\', maxsplit=1)[-1]
            )
        handler.add_static_option(
            'File', 'config', widget,
            display_name='Settings file',
            tooltip=str(self.settings.last_file),
            )
        return handler

    def set_settings(self, new_settings):
        """Set new settings for this instance.

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
        try:                                                                    # TODO: make method that searches through invalid for old values and replaces deprecated ones
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
        if was_busy is not is_busy:
            self._busy = is_busy
            self.busy_changed.emit(self.busy)


class DeviceABC(HardwareABC):
    """Abstract base class of hardware device objects."""

    @abstractmethod
    def list_devices(self):
        """List all devices of this class.

        This method must return a list of SettingsInfo instances. The
        SettingsInfo class is located in the hardwarebase module. Each
        controller is represented by a single SettingsInfo instance. The
        SettingsInfo object must contain a .unique_name and can contain
        .more information as a dict.

        Returns
        -------
        devices : list
            Each element is a SettingsInfo instance containing the name
            of a device and additional information as a dict.
        """
