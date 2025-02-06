"""Module abc of viperleed.guilib.measure.classes.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2024-04-04
Author: Michele Riva
Author: Florian Doerr

This module contains QObjectSettingsErrors, DeviceABCErrors,
QObjectWithError, QObjectWithSettingsABC, HardwareABC and DeviceABC.
QObjectWithError is the base class of QObjectWithSettingsABC.
QObjectWithSettingsABC is the base class of HardwareABC and
MeasurementABC. HardwareABC is the base class of DeviceABC and
SerialABC. DeviceABC is the base class of ControllerABC and CameraABC.

Note that all of the classes defined in this module are subclasses of
QObject. Any subclass of a QObject cannot inherit from a second QObject,
which means one cannot do: NewClass(QObjectSubclass1, QObjectSubclass2)
"""

from abc import ABCMeta
from abc import abstractmethod
from configparser import ConfigParser
from contextlib import contextmanager
from dataclasses import dataclass, field
from operator import itemgetter
from pathlib import Path

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.measure import hardwarebase as base
from viperleed.guilib.measure.classes.settings import (
    NoDefaultSettingsError,
    NoSettingsError,
    SettingsError,
    TooManyDefaultSettingsError,
    ViPErLEEDSettings,
    )
from viperleed.guilib.measure.dialogs.settingsdialog import SettingsHandler


_UNIQUE = qtc.Qt.UniqueConnection


class QMetaABC(ABCMeta, type(qtc.QObject)):
    """Metaclass common to QObject and ABCMeta allowing @abstractmethod."""


class DeviceABCErrors(base.ViPErLEEDErrorEnum):
    """Errors of ViPErLEED devices."""

    DEVICE_NOT_FOUND = (
        120,
        'Unable to connect to device {}. Make sure the device '
        'has power, is connected, and not in use by another program. You can '
        'try replugging it. Make sure to give it enough time to boot up.'
        )


class QObjectSettingsErrors(base.ViPErLEEDErrorEnum):
    """Settings errors of ViPErLEED objects."""

    MISSING_SETTINGS = (100,
                        'Cannot operate without settings. Load an '
                        'appropriate settings file before proceeding.'
                        )
    INVALID_SETTINGS = (101,
                        'Invalid settings. Required settings {} missing or '
                        'value inappropriate. Check configuration file.\n{}'
                        )
    INVALID_SETTING_WITH_FALLBACK = (
        102,
        'Invalid settings. Invalid/unreadable settings value {} for setting '
        '{!r}. Using {} instead. Consider fixing your configuration file.'
        )
    DEFAULT_SETTINGS_CORRUPTED = (103,
                                  'Default settings corrupted. {} '
                                  'Contact the ViPErLEED team to fix '
                                  'your default settings.')

class QObjectWithError(qtc.QObject):                                            # TODO: The Measure class was meant to inherit from this class. Due to double inheritance from QObject this is not possible through standard inheritance.
    """Base class of measurement objects with error detection."""

    # Emitted whenever an error has been detected. Contains
    # information about the error in the form (code, message).
    error_occurred = qtc.pyqtSignal(tuple)

    emit_error = base.emit_error


class QObjectWithSettingsABC(QObjectWithError, metaclass=QMetaABC):
    """Abstract base class of measurement objects with settings.

    Attributes
    ----------
    uses_default_settings : bool
        uses_default_settings remembers if a default settings was
        loaded and can be used during runtime to adjust behaviour
        according to this.

    Private attributes
    ------------------
    _mandatory_settings : tuple
        The _mandatory_settings are the settings that must be present in
        a configuration file for it to be valid. To see how a configration
        file should look like, look at the .ini files in the _defaults
        folder. Do not edit the files in the _defaults folder!
        Each element of the _mandatory_settings should be a tuple with
        the structure (<section>, <option>, (<value1>, <value2>, ...)).
        Section and option detail the required entries in the .ini file,
        while the values are the states the specified option can assume.
        Every mandatory setting must contain a section, and for values
        to be present, an option must be specified too. The allowed
        tuples are therefore: (<section>,), (<section>, <option>), and
        (<section>, <option>, (<value1>, <value2>, ...)).
        To extend the mandatory settings of a class, unpack the
        _mandatory_settings of the parent class into a new
        _mandatory_settings tuple and add the new settings to
        the new tuple like:
        _mandatory_settings = (
            *parent._mandatory_settings,
            (<section>, <option>, (<value1>, <value2>, ...)),
            ...
            )
        To add _mandatory_settings for an instance at runtime, use:
        self._mandatory_settings = (
            *type(self)._mandatory_settings,
            (<section>, <option>, (<value1>, <value2>, ...)),
            ...
            )
    _settings_to_load : ViPErLEEDSettings
        _settings_to_load are the settings that should be loaded
        into _settings via set_settings. If no settings is given,
        _settings_to_load will automatically be the best-matching
        suitable default settings file.
    """

    _mandatory_settings = ()

    def __init__(self, *args, settings=None, **kwargs):
        """Initialise instance.

        This __init__ should be executed at the beginning of
        subclasses. If this is not possible, it has to be executed
        before the settings of an instance are set, i.e., before
        self.settings = self._settings_to_load, or, equivalently,
        self.set_settings(self._settings_to_load).

        Parameters
        ----------
        *args : object
            Positional arguments passed on to the parent class.
        settings : dict or ConfigParser or str or Path or ViPErLEEDSettings or
                   None, optional
            The object settings. If not given or None,
            a suitable default is searched. Default it None.
        **kwargs : object
            Keyword arguments passed on to the parent class.
        """
        # To get the settings of any QObjectWithSettingsABC use the
        # settings property. In a similar manner, to set settings either
        # use the settings property or the set_settings method.
        self._settings = ViPErLEEDSettings()
        self._settings_to_load = settings
        self.uses_default_settings = False
        if not settings:
            # Look for a default only if no settings are given.
            self._settings_to_load = self.find_default_settings()
            self.uses_default_settings = True
        super().__init__(*args, **kwargs)
        self._delayed_errors = []
        self._delay_errors_timer = qtc.QTimer(parent=self)
        self._delay_errors_timer.setSingleShot(True)
        self._delay_errors_timer.setInterval(10)
        self._delay_errors_timer.timeout.connect(self._report_delayed_errors)

    @contextmanager
    def errors_delayed(self):
        """Temporarily make the next code block to report errors later.

        This context manager is useful when the error_occurred signal
        of this QObjectWithError is not yet connected to its final
        slot, for example, when this object is created.
        """
        base.safe_connect(self.error_occurred,
                          self._on_error_delayed,
                          type=_UNIQUE)
        try:
            yield
        finally:
            self.error_occurred.disconnect(self._on_error_delayed)

    @qtc.pyqtSlot(tuple)
    def _on_error_delayed(self, error):
        """Collect errors to be delayed."""
        self._delayed_errors.append(error)
        self._delay_errors_timer.start()

    @qtc.pyqtSlot()
    def _report_delayed_errors(self):
        """Emit errors that we have accumulated."""
        for error in self._delayed_errors:
            self.emit_error(error)
        self._delayed_errors.clear()

    @property
    def settings(self):
        """Return the current settings."""
        return self._settings

    @settings.setter
    def settings(self, new_settings):
        """Set new settings for this instance."""
        self.set_settings(new_settings)

    def check_creating_settings_handler_is_possible(self):                      # TODO: make private and rather use super().get_settings_handler() for implicit check in subclasses
        """Raise if it is not possible to produce a SettingsHandler."""
        if not self.settings:
            # Remember to catch this exception before catching
            # SettingsError. Otherwise NoSettingsError will be
            # swallowed up as it is a subclass of SettingsError.
            raise NoSettingsError(
                f'Instance of {type(self).__name__} tried to return '
                'a settings handler without settings.'
                )
        if self.uses_default_settings:
            raise SettingsError(
                f'Instance of {type(self).__name__} tried to return '
                'a settings handler using default settings.'
                )

    def find_default_settings(self, find_from=None, match_exactly=False):
        """Find default settings for this object.

        Parameters
        ----------
        find_from : SettingsInfo or None, optional
            find_from contains information to look for in the
            configuration files. The way to determine the correct
            settings is up to the is_settings_for_this_class() and
            is_matching_default_settings() methods. Default is None.
            If it is None, subclasses must attempt to determine suitable
            settings without additional information when searching for
            default settings via is_matching_default_settings().
        match_exactly : bool, optional
            Whether find_from should be matched exactly. False means
            the matching of settings files will be less strict. E.g.,
            firmware versions with lower minors may be allowed.
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
        settings = self.find_matching_settings_files(
            find_from, base.DEFAULTS_PATH, match_exactly,
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
    def find_matching_settings_files(cls, obj_info, directory, match_exactly):
        """Find .ini files for obj_info in the tree starting at directory.

        Parameters
        ----------
        obj_info : SettingsInfo or None
            The additional information that should be used to find
            appropriate settings. If it is None, subclasses must attempt
            to determine suitable settings without additional
            information when searching for default settings via
            is_matching_default_settings(). When looking for user
            settings with is_matching_user_settings(), a TypeError will
            be raised if obj_info is None.
        directory : str or Path
            The location in which to look for configuration files.
            Settings files are searched in directory and all its
            subfolders. If directory is the directory containing the
            default settings, is_matching_default_settings is used to
            determine whether a configuration file is appropriate. If it
            is a different folder, is_matching_user_settings is used
            instead.
        match_exactly : bool
            Whether obj_info should be matched exactly.

        Returns
        -------
        obj_settings_files : list
            A list of the paths to settings files that contain appropriate
            settings sorted by how well the settings match from best to
            worst.
        """
        directory = Path(directory).resolve()
        default = True if directory == base.DEFAULTS_PATH else False
        settings_files = directory.glob('**/*.ini')
        if not default:
            # Filter out default settings.
            settings_files = [file for file in settings_files
                              if '_defaults' not in str(file)]

        files_and_scores = []
        is_matching = (cls.is_matching_default_settings if default
                       else cls.is_matching_user_settings)
        for settings_file in settings_files:
            config = ViPErLEEDSettings.from_settings(settings_file)
            if not cls.is_settings_for_this_class(config):
                continue
            score = is_matching(obj_info, config, match_exactly)
            if score:
                files_and_scores.append((settings_file, score))

        if not files_and_scores:
            return []
        # Sort by score, best first, and discard the scores
        obj_settings_files, _ = zip(
            *sorted(files_and_scores, key=itemgetter(1), reverse=True)
            )
        return obj_settings_files

    @classmethod
    @abstractmethod
    def is_matching_default_settings(cls, obj_info, config, match_exactly):
        """Determine if a default `config` file is for this instance.

        This method must use obj_info to determine if a config is
        suitable for the instance. To perform this check one can
        .get values from config and compare the information from
        obj_info.more to it. In contrast to is_matching_user_settings,
        this method can be called before the connected hardware is
        known. Therefore, all information to compare to must come
        from the object class itself. If match_exactly is True, the
        default has to be a perfect match. Only one default at a time
        can fulfill this criterion.

        Parameters
        ----------
        obj_info : SettingsInfo or None
            The information that should be used to check `config`.
        config : ConfigParser
            The settings to check.
        match_exactly : bool
            Whether obj_info should be matched exactly.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used to sort the detected settings.
            Larger values in the tuple indicate a higher degree of
            conformity. The order of the items in the tuple is the
            order of their significance. This return value is used
            to determine the best-matching settings files when
            multiple files are found. An empty tuple signifies that
            `config` does not match the requirements.
        """

    @classmethod
    @abstractmethod
    def is_matching_user_settings(cls, obj_info, config, match_exactly):
        """Determine if a `config` file is for this instance.

        Subclasses must extend this method by calling
        super().is_matching_user_settings at the beginning.

        This method must use obj_info to determine if a config is
        suitable for the instance. To perform this check one can
        .get values from config and compare the information from
        obj_info.more and obj_info.unique_name to it. If
        match_exactly is True, the config has to be a perfect match
        and must allow full functionality of the connected hardware.

        Parameters
        ----------
        obj_info : SettingsInfo
            The information that should be used to check `config`.
        config : ConfigParser
            The settings to check.
        match_exactly : bool
            Whether obj_info should be matched exactly.

        Returns
        -------
        sorting_info : tuple
            A tuple that can be used to sort the detected settings.
            Larger values in the tuple indicate a higher degree of
            conformity. The order of the items in the tuple is the
            order of their significance. This return value is used
            to determine the best-matching settings files when
            multiple files are found. An empty tuple signifies that
            `config` does not match the requirements.
        """
        if not isinstance(obj_info, SettingsInfo):
            raise TypeError(
                f'obj_info should be of type SettingsInfo '
                f'but is of type {type(obj_info).__name__}.'
                )
        return tuple()

    @classmethod
    @abstractmethod
    def is_settings_for_this_class(cls, config):
        """Determine if a `config` file is for this class.

        This method must use cls attributes to determine if
        a config is suitable for the class.

        Parameters
        ----------
        config : ConfigParser
            The settings to check.

        Returns
        -------
        is_suitable : bool
            True if the settings file is for the class.
        """

    def are_settings_invalid(self, settings):
        """Check if there are any invalid settings.

        Subclasses may add additional mandatory settings at
        runtime. See the documentation of the _mandatory_settings
        attribute for how to do this. The base implementation will
        check if all of the _mandatory_settings are present in the
        provided settings.

        Parameters
        ----------
        settings : ViPErLEEDSettings
            The new settings.

        Returns
        -------
        invalid_settings : list of tuples
            Invalid _mandatory_settings of self as a list of tuples.
            The first entry in each tuple can be either '<section>',
            '<section>/<option>', or
            '<section>/<option> not one of <value1>, <value2>, ...'.
            Further optional entries may be added by subclasses. They
            specify additional information on what is wrong with each
            invalid setting.
        """
        return [(invalid,) for invalid in
                settings.has_settings(*self._mandatory_settings)]

    @abstractmethod
    def get_settings_handler(self):
        """Return a SettingsHandler object for displaying settings.

        This method must be extended in subclasses, i.e., do
        handler = super().get_settings_handler(), and then add
        appropriate sections and/or options to it using the
        handler.add_section and handler.add_option methods.

        The base-class implementation returns a handler that
        contains the location of the settings file.

        Use the QNoDefaultPushButton from the basewidgets module
        in order to prevent any button from being set as the
        default button of the dialog.

        Returns
        -------
        handler : SettingsHandler
            The handler used in a SettingsDialog to display the
            settings of this instance to users.
        """
        self.check_creating_settings_handler_is_possible()
        handler = SettingsHandler(self.settings)
        return handler

    @qtc.pyqtSlot(object)
    def set_settings(self, new_settings):
        """Set new settings for this instance.

        Check and set settings. This method is used in
        the setter of the settings property.

        Parameters
        ----------
        new_settings : dict or ConfigParser or str or Path or ViPErLEEDSettings
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
            If an element of the mandatory settings is None or
            is not a Sequence with length <= 3.

        Emits
        -----
        QObjectSettingsErrors.MISSING_SETTINGS
            If the settings have not been found
            or if new_settings was None.
        QObjectSettingsErrors.INVALID_SETTINGS
            If any element of the new_settings does
            not fit the mandatory settings.
        """
        try:                                                                    # TODO: #242 make method that searches through invalid for old values and replaces deprecated ones, make it a method of the ViPErLEEDSettings class
            new_settings = ViPErLEEDSettings.from_settings(new_settings)
        except (ValueError, NoSettingsError):
            self.emit_error(QObjectSettingsErrors.MISSING_SETTINGS)
            return False
        invalid = self.are_settings_invalid(new_settings)
        for missing, *info in invalid:
            self.emit_error(QObjectSettingsErrors.INVALID_SETTINGS,
                       missing, ' '.join(info))
            return False

        self._settings = new_settings
        return True


# DISABLE: seems a bug in pylint. This is also an ABC.
# pylint: disable-next=abstract-method
class HardwareABC(QObjectWithSettingsABC):
    """Abstract base class of hardware-related objects."""

    # Emitted whenever the busy state of the device changes.
    # Contains the new busy state of the device.
    busy_changed = qtc.pyqtSignal(bool)

    # Emitted right after the hardware connection status has changed.
    # Cointains the new connection status of the hardware.
    connection_changed = qtc.pyqtSignal(bool)

    def __init__(self, *args, **kwargs):
        """Initialise instance."""
        # Busy state of the device.
        self._busy = False
        super().__init__(*args, **kwargs)

    @property
    def busy(self):
        """Return busy state of instance."""
        return self._get_busy()

    @busy.setter
    def busy(self, is_busy):
        """Set busy state of instance.

        Note that 'self.busy =' will be performed in the calling
        thread. To avoid this, one can invoke the set_busy slot.

        Parameters
        ----------
        is_busy : bool
            True if the instance is busy.

        Emits
        -----
        busy_changed
            If the busy state of the instance changed.
        """
        self.set_busy(is_busy)

    def _get_busy(self):
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
        self._busy = is_busy
        if was_busy != is_busy:
            self.busy_changed.emit(self.busy)


class DeviceABC(HardwareABC):
    """Abstract base class of hardware device objects."""

    @abstractmethod
    def list_devices(self):
        """List all devices of this class.

        This method must return a list of SettingsInfo instances. Each
        device is represented by a single SettingsInfo instance. The
        SettingsInfo object must contain a .unique_name, a
        .hardware_interface boolean, and can contain .more information
        as a dict. The information contained within a SettingsInfo must
        be enough to determine settings files that contain the correct
        settings for this device. Subclasses should raise a
        DefaultSettingsError if they fail to create instances from the
        settings in the DEFAULTS_PATH.

        Returns
        -------
        devices : list
            Each item is a SettingsInfo instance containing the name
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
        Unique name identifying the device.
    hardware_interface : bool
        Whether the device has a hardware interface or not.
    more : dict
        Extra, optional, information about the device.
    """
    unique_name: str
    hardware_interface: bool
    more: dict = field(default_factory=dict)

    def __post_init__(self):
        """Check that we have the correct attribute types."""
        if not isinstance(self.unique_name, str):
            raise TypeError(f'{type(self).__name__}: '
                            'unique_name must be a string')
        if not isinstance(self.hardware_interface, bool):
            raise TypeError(f'{type(self).__name__}: '
                            'hardware_interface must be a bool')
