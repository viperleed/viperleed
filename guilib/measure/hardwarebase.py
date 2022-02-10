"""Module hardwarebase of viperleed.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-07-08
Author: Michele Riva
Author: Florian Doerr

This module contains utility functions and classes shared
by multiple ViPErLEED hardware objects.
"""

from abc import ABCMeta
from configparser import ConfigParser
import inspect
from pathlib import Path
import enum
import sys

from PyQt5 import QtCore as qtc

################################## FUNCTIONS ##################################

def class_from_name(package, class_name):
    """Return the serial class given its name.

    Parameters
    ----------
    package : str
        Name of the viperleed.guilib package where
        the class should be looked for. package
        should be the name of the package relative
        to the main measure plugin folder.
    class_name : str
        Name of the class to be returned

    Raises
    ------
    AttributeError
        If package is not a valid package
    ValueError
        If class_name could not be found in package
    """
    try:
        getattr(sys.modules[__package__], package)
    except AttributeError as err:
        raise AttributeError(f"{__package__} does not contain "
                             f"a package named {package}") from err

    package_name = f"{__package__}.{package}"
    try:
        cls = getattr(sys.modules[package_name], class_name)
    except AttributeError as err:
        raise ValueError(
            f"No {package} class named {class_name} found."
            ) from err
    return cls


# too-complex, too-many-branches  -> may be solved by moving
# the str/Path/dict cases into the new HardwareSettings class.
# This whole function can then be a method of the class.
def config_has_sections_and_options(caller, config, mandatory_settings):
    """Make sure settings are fine, and return it as a ConfigParser.

    Dictionaries are converted into ConfigParsers and checked.
    ConfigParsers are only checked. Strings are converted into
    path objects and paths are used to access files from which
    data is read into ConfigParsers and checked.

    Parameters
    ----------
    caller : object
        The object calling this method
    config : dict, ConfigParser, str, path or None
        The configuration to be checked. It is advisable to
        NOT pass a dict, as the return value will be a copy,
        and will not be up to date among objects sharing the
        same configuration.
    mandatory_settings : Sequence
        Each element is a tuple of the forms
        (<section>, )
            Checks that config contains the section <section>.
            No options are checked.
        (<section>, <option>)
            Checks that the config contains <option> in the
            existing section <section>. No values are checked
            for <option>.
        (<section>, <option>, <admissible_values>)
            Checks that the config contains <option> in the
            existing section <section>, and that <option> is
            one of the <admissible_values>. In this case,
            <admissible_values> is a Sequence of strings.

    Returns
    -------
    config : ConfigParser or None
        The same config, but as a ConfigParser, provided
        the settings are OK. Returns None if some of the
        mandatory_settings is invalid or if the original
        config was None.
    invalid_settings : list
        Invalid mandatory_settings of config. None in case
        settings are OK, list if settings are invalid, True
        if config was None.

    Raises
    ------
    TypeError
        If settings is not a dict, str, path, or a
        ConfigParser or mandatory_settings contains invalid
        data (i.e., one of the entries is not a Sequence
        with length <= 3).
    FileNotFoundError
        If the requested configuration file could not be
        found/opened.
    """
    if config is None:
        return None, mandatory_settings

    if not isinstance(config, (dict, ConfigParser, str, Path)):
        raise TypeError(
            f"{caller.__class__.__name__}: invalid type "
            f"{type(config).__name__} for settings. "
            "Should be 'dict', 'str' or 'ConfigParser'."
            )

    if isinstance(config, str):
        config = Path(config)
    if isinstance(config, Path):
        if not config.is_file():
            raise FileNotFoundError(f"Could not open/read file: {config}. "
                                    "Check if this file exists and is in the "
                                    "correct folder.")
        tmp_config = ConfigParser(comment_prefixes='/', allow_no_value=True)
        tmp_config.read(config)
        config = tmp_config  # pylint: disable=redefined-variable-type

    if isinstance(config, dict):
        tmp_config = ConfigParser()
        tmp_config.read_dict(config)
        config = tmp_config

    invalid_settings = []
    for setting in mandatory_settings:
        if not setting or len(setting) > 3:
            raise TypeError(f"Invalid mandatory setting {setting}. "
                            f"with length {len(setting)}. Expected "
                            "length <= 3.")
        # (<section>,)
        if len(setting) == 1:
            if not config.has_section(setting[0]):
                invalid_settings.append(setting)
            continue

        # (<section>, <option>) or (<section>, <option>, <admissible>)
        section, option = setting[:2]
        if not config.has_option(section, option):
            invalid_settings.append(setting)
            continue

        # (<section>, <option>, <admissible>)
        if len(setting) == 3:
            admissible_values = setting[2]
            value = config.get(section, option)
            if value not in admissible_values:
                invalid_settings.append(setting)

    if invalid_settings:
        config = None
    return config, invalid_settings


def emit_error(sender, error, *msg_args, **msg_kwargs):
    """Emit a ViPErLEEDErrorEnum-like error.

    Parameters
    ----------
    sender : object
        The instance in which the error occurred. Must have
        an error_occurred pyqtSignal.
    error : tuple or ViPErLEEDErrorEnum
        The error to be emitted. Should be a 2-tuple of the
        form (error_code, error_message).
    *msg_args : object
        Extra info to be inserted in the error message as
        positional arguments to str.format
    **msg_kwargs : object
        Extra info to be inserted in the error message as
        keyword arguments.

    Raises
    ------
    TypeError
        If sender does not have an error_occurred signal,
        or if error is not a Sequence.
    ValueError
        If error is not a 2-element tuple
    """
    if not hasattr(sender, 'error_occurred'):
        raise TypeError(f"Object {sender} has no error_occurred signal "
                        "to emit. Probably an inappropriate type")
    try:
        _ = len(error)
    except TypeError as err:
        raise TypeError(f"Invalid error {error} cannot be emitted") from err

    if len(error) != 2:
        raise ValueError(f"Invalid error {error} cannot be emitted")

    if not msg_args and not msg_kwargs:
        sender.error_occurred.emit(error)
        return

    error_code, error_msg = error
    error_msg = error_msg.format(*msg_args, **msg_kwargs)
    sender.error_occurred.emit((error_code, error_msg))


def get_devices(package):
    """Return all supported and available devices from a package.

    Supported devices are those for which a driver class exists.
    Each class should implement a .list_devices() method.

    Parameters
    ----------
    package : str
        Name of the viperleed.measure package for which supported,
        available devices should be listed. Typically 'controller'
        or 'camera'.

    Returns
    -------
    devices : dict
        Dictionary of available, supported devices. Keys are
        names of available devices, values are the driver classes
        that can be used to handle the devices (one class per
        device).

    Raises
    ------
    AttributeError
        If package is not one of the guilib.measure subpackages.
    """
    try:
        getattr(sys.modules[__package__], package)
    except AttributeError as err:
        raise AttributeError(f"{__package__} does not contain "
                             f"a package named {package}") from err

    package_name = f"{__package__}.{package}"
    devices = {}
    for _, cls in inspect.getmembers(sys.modules[package_name],
                                     inspect.isclass):
        if hasattr(cls, 'list_devices'):
            dummy_instance = cls()
            dev_list = dummy_instance.list_devices()
            for dev_name in dev_list:
                devices[dev_name] = cls
    return devices


################################### CLASSES ###################################


class QMetaABC(ABCMeta, type(qtc.QObject)):
    """Metaclass common to QObject and ABCMeta allowing @abstractmethod."""


class ViPErLEEDErrorEnum(tuple, enum.Enum):
    """Base class for ViPErLEED hardware errors.

    Each Enum value is a tuple of the form
        (error_code, error_description)
    where error_code is an integer, and error_description a string
    that will be shown to the user when the error occurs. The string
    may also contain formatting directives, that should be .format()ted
    by the controller or by the GUI.
    """
    @classmethod
    def as_dict(cls):
        """Return a dictionary of error codes and error names."""
        return {el[0]: el.name for el in cls}

    @classmethod
    def from_code(cls, code):
        """Return the ViPErLEEDErrorEnum with a given code.

        Parameters
        ----------
        code : int
            The code of the error

        Returns
        -------
        ViPErLEEDErrorEnum.value
            The hardware error itself, which can be unpacked
            into a tuple of the form (code, error_description)

        Raises
        ------
        TypeError
            If code is not an int
        AttributeError
            If code is an unknown error code
        """
        if not isinstance(code, int):
            raise TypeError(f"{cls.__name__}.from_code(code): type of code "
                            f"should be 'int', not {type(code).__name__!r}.")
        try:
            error_name = cls.as_dict()[code]
        except KeyError as err:
            raise AttributeError(
                f"Unknown {cls.__name__} error with code {code}"
                ) from err
        return getattr(cls, error_name)

    @property
    def code(self):
        """Return the code of this error."""
        return self.value[0]
