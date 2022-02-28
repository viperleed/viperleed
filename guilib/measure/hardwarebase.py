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
import inspect
from pathlib import Path
import enum
import sys

from PyQt5 import (QtWidgets as qtw, QtCore as qtc)

from viperleed.guilib import measure as vpr_measure


DEFAULTS_PATH = (Path(inspect.getfile(vpr_measure)).parent                      # TODO: get it from system settings
                 / '_defaults')


################################## FUNCTIONS ##################################

def class_from_name(package, class_name):
    """Return a class given its name.

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


def get_device_config(device_name, directory=DEFAULTS_PATH,
                      prompt_if_invalid=True, parent_widget=None):
    """Return the configuration file for a specific device.

    Only configuration files with a .ini suffix are considered.

    Parameters
    ----------
    device_name : str
        String to be looked up in the configuration files
        to identify that a file is meant for the device.
    directory : str or Path, optional
        The base of the directory tree in which the
        configuration file should be looked for. The
        search is recursive. Default is the path to
        the _defaults drectory.
    prompt_if_invalid : bool, optional
        In case the search for a config file failed, pop up
        a dialog asking the user for input. The search is
        considered failed in case no config file was found,
        or if multiple config files matched the criterion.
        Defaul is True.
    parent_widget : QWidget or None, optional
        The parent widget of the pop up. Unless parent_widget
        is None, the pop up will be blocking user events for
        parent_widget till the user closes it. Used only
        if prompt_if_invalid is True. Default is None.

    Returns
    -------
    path_to_config : Path or None
        The path to the only configuration file succesfully
        found. None if no configuraiton file was found (either
        because prompt_if_invalid is False, or because the
        user dismissed the pop-up).
    """
    directory = Path(directory)
    config_files = [f for f in directory.glob('**/*')
                    if f.is_file() and f.suffix == '.ini']
    device_config_files = []
    for config_name in config_files:
        with open(config_name, 'r', encoding='utf-8') as config_file:
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
                return get_device_config(device_name, directory=new_path)
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


def safe_connect(signal, slot, **kwargs):
    """Connect signal to slot swallowing TypeError."""
    if not callable(slot) and not isinstance(slot, qtc.pyqtSignal):
        raise TypeError(f"{slot} is not a valid slot")
    for arg in kwargs:
        if arg not in ("type", "no_receiver_check"):
            raise TypeError("safe_connect got an unexpeced "
                            f"keyord argument {arg!r}")
    try:
        signal.connect(slot, **kwargs)
    except TypeError:
        # typically happens when type=QtCore.Qt.UniquConnection
        # and there is already a valid connection between signal
        # and slot.
        pass


def safe_disconnect(signal, *slot):
    """Disconnect signal from optional slot swallowing TypeError."""
    if (slot and not callable(slot[0])
            and not isinstance(slot[0], qtc.pyqtSignal)):
        raise TypeError(f"{slot} is not a valid slot")
    if len(slot) > 1:
        raise TypeError("safe_disconnect: slot should be a "
                        "single slot when given.")
    try:
        signal.disconnect(*slot)
    except TypeError:
        # Was not connected to begin with
        pass


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
