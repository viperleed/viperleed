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

from PyQt5 import (QtWidgets as qtw, QtCore as qtc)

from viperleed.guilib import measure as vpr_measure


DEFAULT_CONFIG_PATH = (Path(inspect.getfile(vpr_measure)).parent
                       / '_defaults')


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
        config = tmp_config

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
        elif len(setting) == 1:
            # (<section>,)
            if not config.has_section(setting[0]):
                invalid_settings.append(setting[0])
                continue
        elif len(setting) in (2, 3):
            # (<section>, <option>) or (<section>, <option>, <admissible>)
            section, option = setting[:2]
            if not config.has_option(section, option):
                invalid_settings.append("/".join(setting))
                continue
            if len(setting) == 3:
                # (<section>, <option>, <admissible>)
                admissible_values = setting[2]
                value = config.get(section, option)
                if value not in admissible_values:
                    invalid_settings.append("/".join(setting))

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
    except TypeError:
        raise TypeError(f"Invalid error {error} cannot be emitted")

    if len(error) != 2:
        raise ValueError(f"Invalid error {error} cannot be emitted")

    if not msg_args and not msg_kwargs:
        sender.error_occurred.emit(error)
        return

    error_code, error_msg = error
    error_msg = error_msg.format(*msg_args, **msg_kwargs)
    sender.error_occurred.emit((error_code, error_msg))


def get_device_config(device_name, directory=DEFAULT_CONFIG_PATH,
                      parent_widget=None, prompt_if_invalid=True):
    """Return the configuration file for a specific device."""
    directory = Path(directory)
    config_files = [f for f in directory.glob('**/*')
                    if f.is_file() and f.suffix == '.ini']
    device_config_files = []
    for config_name in config_files:
        with open(config_name, 'r') as config_file:
            if device_name in config_file.read():
                device_config_files.append(config_name)

    if device_config_files and len(device_config_files) == 1:
        # Found exactly one config file
        return device_config_files[0]

    if not prompt_if_invalid:
        return None

    if not device_config_files:
        msg_box = qtw.QMessageBox(parent=parent_widget)
        msg_box.setWindowTitle("No settings file found")
        msg_box.setText(
            f"Directory {directory} and its subfolders do not contain "
            f"any settings file for device {device_name}. Select "
            "a different directory."
            )
        msg_box.setIcon(msg_box.Critical)
        msg_box.addButton(msg_box.Cancel)
        btn = msg_box.addButton("Select path", msg_box.ActionRole)
        msg_box.exec_()
        if msg_box.clickedButton() is btn:
            new_path = qtw.QFileDialog.getExistingDirectory(
                parent=parent_widget,
                caption="Choose directory of device settings",
                directory=str(directory)
                )
            if new_path:
                return get_device_config(device_name, new_path)
        return None

    # Found multiple config files that match.
    # Let the use pick which one to use
    names = [f.name for f in device_config_files]
    dropdown = vpr_measure.dialogs.DropdownDialog(
        "Found multiple settings files",
        "Found multiple settings files for device "
        f"{device_name} in {directory} and subfolders.\n"
        "Select which one should be used:",
        names, parent=parent_widget
        )
    config = None
    if dropdown.exec_() == dropdown.Apply:
        config = device_config_files[names.index(dropdown.selection)]
    return config

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


class QMetaABC(type(qtc.QObject), ABCMeta):
    """Metaclass common to QObject and ABCMeta allowing @abstractmethod."""

    pass


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
