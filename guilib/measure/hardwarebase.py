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
from dataclasses import dataclass, field
import enum
import inspect
from pathlib import Path
import re
import sys

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.guilib.dialogs.dropdowndialog import DropdownDialog

# TODO: not nice. Also, there's two places where the _defaults
# path is used. Here and in classes.settings. However, due to circular
# imports I cannot use either for the other one. Probably better
# to move the path in __init__?
DEFAULTS_PATH = Path(__file__).parent / '_defaults'


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
    try:
        sender.parent()
    except RuntimeError:
        # C++ QObject deleted
        return

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


_W_DASH_DOT_RE = re.compile(r'[^-\w.]')
_MULTI_DASH_RE = re.compile(r'[-]+')
_MULTI_UNDERSCORE_RE = re.compile(r'[_]+')

def as_valid_filename(name):
    """Return a valid filename from name."""
    # - Spaces to underscores;
    # - Keep alphanumeric, dashes, underscores, dots;
    # - Keep ASCII only
    # - Replace multiple dashes/underscores with singles;
    # - Remove dashes/underscores from the ends;
    fname = name.replace(' ', '_')
    fname = _W_DASH_DOT_RE.sub('', fname)
    fname = fname.encode('ascii', 'ignore').decode('ascii')
    fname = _MULTI_DASH_RE.sub('-', fname)
    return _MULTI_UNDERSCORE_RE.sub('_', fname).strip('-_')


_MATCH_SQUARES = re.compile(r'(\[)(\[*(?:[^\[\]]*|\[[^\]]*\])*\]*)(\])')
# All capture groups appear in re.split
# (\[)                 # open bracket capture group
# (                    # capture group
#     \[*              # maybe an open bracket
#     (?:              # maybe a non-capture group containing
#         [^\[\]]*     # any except '[' ']'
#         |            # or
#         \[           # an open bracket,
#         [^\]]*       # anything that is not a closed bracket
#         \]           # a closed bracket
#     )*
#     \]*              # maybe a closed bracket
# )
# (\])                 # closed bracket capture group

def _device_name_re(name):
    """Return a re.Pattern object that matches name tolerantly.

    Returns
    -------
    pattern : re.Pattern
        A compiled regular expression to use for matching a name
        tolerantly. The regular expression will match <text> if:
        (1) <text> contains exactly name, (2) <text> contains name
        with its square brackets removed, (3) <text> contains name,
        but possibly misses parts of name included in square brackets;
        the text in brackets may be missing, or may be different;
        (4) <text> may contain parts in between square brackets
        where there are spaces in name, if name does not contain any
        quare bracket.

        Matching only occurs toward the end of <text>, apart from an
        arbitrary numbers of spaces
    """
    # (1) Match exactly
    _patterns = [re.escape(name)]

    # Now things depend on whether we have '[' in name:
    if '[' not in name:
        # Allow for squares to appear anywhere where there is a space
        _patterns.append(re.escape(name).replace('\\ ', '\\ (\\[.*\\]\\ )?'))
    else:
        # (2) Same as name, apart from '[' or ']' characters
        _patterns.append(re.escape(name.replace('[', '').replace(']', '')))
        # (3) All parts in name, but possibly excluding those in squares
        #     This is a bit more involved as we want to replace only
        #     outer square brackets with an optional capture group, which
        #     can capture any character. The group can appear zero or once.
        #     Moreover, we need to escape all characters out of the group,
        #     and replace ' ' with ' ?' in case parts are missing.
        _pattern3 = []
        _in_brackets = False
        _space = re.escape(' ')
        for part in _MATCH_SQUARES.split(name):
            if not part:
                continue
            if not _in_brackets and part not in '[]':
                _pattern3.append(re.escape(part).replace(_space, _space + '?'))
            elif _in_brackets and part != ']':
                _pattern3.append(f"({_MATCH_SQUARES.pattern})?")
            if part == '[':
                _in_brackets = True
            elif part == ']':
                _in_brackets = False
        _patterns.append(''.join(_pattern3))
    # All patterns should match towards the end of the line,
    # and can have any characters before
    return re.compile("|".join(r'^.*' + p + r'\s*$' for p in _patterns))


def _get_object_config_not_found(obj_cls, obj_info, **kwargs):
    """Return a device config when it was not found in the original path.

    This function is only called as part of get_object_config
    in case the no config file was found with the original arguments
    given to get_object_config. Prompts the user for an action.

    Parameters
    ----------
    obj_cls : object
        The class of the object to get settings for.
    obj_info : str, SettingsInfo
        If obj_info is a string, it is the string to be looked up in the
        configuration files to identify that a file is meant for the
        device. If it is a SettingsInfo, the way to determine the correct
        settings is up to the reimplementation of find_matching_configs
        in obj_cls.
    **kwargs : dict
        The same arguments given to get_object_config

    Returns
    -------
    path_to_config : Path, "", or None
        The path to the only configuration file successfully
        found. None or "" if no configuration file was found.
        None is always returned if prompt_if_invalid is False,
        or because the user dismissed the pop-up. The empty
        string is returned if an alternative option was given
        as a third_btn_text parameter, but the user anyway
        dismissed the dialog.
    """
    parent_widget = kwargs.get("parent_widget", None)
    directory = kwargs.get("directory", DEFAULTS_PATH)
    third_btn_text = kwargs.get("third_btn_text", "")

    if isinstance(obj_info, str):
        obj_name = obj_info
    else:
        obj_name = obj_info.unique_name

    msg_box = qtw.QMessageBox(parent=parent_widget)
    msg_box.setWindowTitle("No settings file found")
    msg_box.setText(
        f"Directory {directory} and its subfolders do not contain "
        f"any settings file for device {obj_name}. Select "
        "a different directory."
        )
    msg_box.setIcon(msg_box.Warning)
    btn = msg_box.addButton("Select path", msg_box.ActionRole)
    third_btn = None
    if third_btn_text:
        third_btn = msg_box.addButton(third_btn_text, msg_box.AcceptRole)
    msg_box.addButton(msg_box.Cancel)
    msg_box.exec_()
    msg_box.setParent(None)  # py garbage collector will take care
    _clicked = msg_box.clickedButton()
    if msg_box.clickedButton() is btn:
        new_path = qtw.QFileDialog.getExistingDirectory(
            parent=parent_widget,
            caption="Choose directory of device settings",
            directory=str(directory)
            )
        if new_path:
            kwargs["directory"] = new_path
            return get_object_config(obj_cls, obj_info, **kwargs)
    if third_btn and _clicked is not third_btn:
        # User had another option but dismissed the dialog
        return ""
    return None


def get_object_config(obj_cls, obj_info, **kwargs):                             # TODO: This currently may execute in secondary threads. Split settings search from reporting errors and multiple-settings-found.
    """Return the configuration file for a specific device.

    Only configuration files with a .ini suffix are considered.
    This method must be executed in the main GUI thread.

    Parameters
    ----------
    obj_cls : object
        The class of the object to get settings for.
    obj_info : str, SettingsInfo
        If obj_info is a string, it is the string to be looked up in the
        configuration files to identify that a file is meant for the
        device. If it is a SettingsInfo, the way to determine the correct
        settings is up to the reimplementation of find_matching_configs
        in obj_cls.
    **kwargs : dict, optional
        directory : str or Path, optional
            The base of the directory tree in which the
            configuration file should be looked for. The
            search is recursive. Default is the path to
            the _defaults directory.
        tolerant_match : bool, optional
            Whether obj_info as str should be looked up tolerantly
            or not. If False, obj_info is matched exactly,
            otherwise parts of obj_info within square brackets
            are ignored. Default is True.
        prompt_if_invalid : bool, optional
            In case the search for a config file failed, pop up
            a dialog asking the user for input. The search is
            considered failed in case no config file was found,
            or if multiple config files matched the criterion.
            Default is True.
        parent_widget : QWidget or None, optional
            The parent widget of the pop up. Unless parent_widget
            is None, the pop up will be blocking user events for
            parent_widget till the user closes it. Used only
            if prompt_if_invalid is True. Default is None.
        third_btn_text : str, optional
            If given and not empty, the string is used as the text
            for an extra button with "accept role" in the did-not-
            find-a-config case. This function may return an empty
            string only if this is given and if the user dismisses
            the dialog. Default is an empty string.

    Returns
    -------
    path_to_config : Path, "", or None
        The path to the only configuration file successfully
        found. None or "" if no configuration file was found.
        None is always returned if prompt_if_invalid is False,
        or because the user dismissed the pop-up. The empty
        string is returned if an alternative option was given
        as a third_btn_text parameter, but the user anyway
        dismissed the dialog.
    """
    directory = kwargs.get("directory", DEFAULTS_PATH)
    tolerant_match = kwargs.get("tolerant_match", True)
    prompt_if_invalid = kwargs.get("prompt_if_invalid", True)
    parent_widget = kwargs.get("parent_widget", None)

    device_config_files = obj_cls.find_matching_configs(
                            obj_info, directory, tolerant_match
                            )

    if device_config_files and len(device_config_files) == 1:
        # Found exactly one config file
        return device_config_files[0]

    if not prompt_if_invalid:
        return None

    if not device_config_files:
        return _get_object_config_not_found(obj_cls, obj_info, **kwargs)

    # Found multiple config files that match.
    # Let the user pick which one to use
    if isinstance(obj_info, str):
        obj_name = obj_info
    else:
        obj_name = obj_info.unique_name
    names = [f.name for f in device_config_files]
    dropdown = DropdownDialog(
        "Found multiple settings files",
        "Found multiple settings files for "
        f"{obj_name} in {directory} and subfolders.\n"
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
        names of available devices, values are tuples containing
        the driver classes that can be used to handle the devices
        (one class per device) and the SettingsInfo from list_devices.

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
            for device in dev_list:
                devices[device.unique_name] = (cls, device)
    return devices


def safe_connect(signal, slot, **kwargs):
    """Connect signal to slot swallowing TypeError."""
    if not callable(slot) and not isinstance(slot, qtc.pyqtSignal):
        raise TypeError(f"{slot} is not a valid slot")
    for arg in kwargs:
        if arg not in ("type", "no_receiver_check"):
            raise TypeError("safe_connect got an unexpected "
                            f"keyword argument {arg!r}")
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


@dataclass
class SettingsInfo:
    """A container for information about settings of objects.

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


_DOTS_OR_DIGITS = re.compile(r"[.0-9]+")


class Version:
    """Simple class to easily handle versions."""

    def __init__(self, major, minor=-1, patch=-1):
        """Initialize instance.

        Parameters
        ----------
        major : int, float, str, or Version
            When an integer, interpret as the "major" part of
            version; when a string, it should be of the form
            "<major>.<minor>.<patch>", with minor and patch
            optional, and each field be int-able; when a float,
            interpret the whole as the major and the decimal
            as the minor.
        minor : int, optional
            Minor portion of the version. This is discarded if major
            is given as a string, a float, or a Version. If negative
            or not given, it is interpreted as "zero".
        minor : int, optional
            Patch portion of the version. This is discarded if major
            is given as a string, a float, or a Version. If negative
            or not given, it is interpreted as "zero".

        Raises
        ------
        TypeError
            If any of the inputs is not an integer. If major is given
            as a string, the exception is raised if one of the parts
            is not int-able.
        ValueError
            If major is an empty string, if it does not conform to
            the <number>[.<number>[.<number>]] syntax, if <major>
            is negative, if <minor> is not given but <patch> is.
        """
        _name = self.__class__.__name__
        if isinstance(major, float):
            major = str(major)
        if isinstance(major, str):
            major, minor, patch = self.__from_string(major)
        elif isinstance(major, Version):
            major, minor, patch = major

        if not all(isinstance(p, int) for p in (major, minor, patch)):
            raise TypeError(f"{_name}: Invalid type for one of the parts")
        if major < 0:
            raise ValueError(f"{_name}: version cannot be negative")

        # pylint: disable=chained-comparison
        # pylint would probably like the next one to be
        # minor < 0 < patch, but it feels less readable
        if patch > 0 and minor < 0:
            raise ValueError(f"{_name}: Cannot have 'patch' without 'minor'")
        self.__has_parts = (True, minor > 0, patch > 0)
        self._parts = (major, max(minor, 0), max(patch, 0))

    @property
    def major(self):
        """Return the major part of self."""
        return self._parts[0]

    @property
    def minor(self):
        """Return the minor part of self (if present) or RuntimeError."""
        if not self.__has_parts[1]:
            raise RuntimeError(f"{self.__class__.__name__} has no minor.")
        return self._parts[1]

    @property
    def patch(self):
        """Return the minor part of self (if present) or RuntimeError."""
        if not self.__has_parts[2]:
            raise RuntimeError(f"{self.__class__.__name__} has no patch.")
        return self._parts[2]

    def __eq__(self, other):
        """Return True if self == other, NotImplemented otherwise."""
        if not isinstance(other, self.__class__):
            try:
                other = self.__class__(other)
            except TypeError:
                return NotImplemented
        if self._parts == other._parts:
            return True
        return NotImplemented

    def __iter__(self):
        """Return an iterable version of self."""
        return iter(self._parts)

    def __lt__(self, other):
        """Return whether self < other."""
        if not isinstance(other, self.__class__):
            try:
                other = self.__class__(other)
            except TypeError:
                return NotImplemented
        return self._parts < other._parts

    def __le__(self, other):
        """Return whether self <= other."""
        if not isinstance(other, self.__class__):
            try:
                other = self.__class__(other)
            except TypeError:
                return NotImplemented
        return self._parts <= other._parts

    def __gt__(self, other):
        """Return whether self > other."""
        if not isinstance(other, self.__class__):
            try:
                other = self.__class__(other)
            except TypeError:
                return NotImplemented
        return self._parts > other._parts

    def __ge__(self, other):
        """Return whether self > other."""
        if not isinstance(other, self.__class__):
            try:
                other = self.__class__(other)
            except TypeError:
                return NotImplemented
        return self._parts >= other._parts

    def __str__(self):
        """Return a string representation of self."""
        return ".".join(str(p)
                        for p, had_part in zip(self._parts, self.__has_parts)
                        if had_part)

    def __repr__(self):
        """Return a representation of self."""
        return f"{self.__class__.__name__}({self})"

    def __from_string(self, _as_string):
        """Return major, minor, patch from a string."""
        _name = self.__class__.__name__
        if not _as_string:
            # Empty version
            raise ValueError(f"{_name}: version cannot be an empty string")
        _as_string = _as_string.replace('v', '')
        if not _DOTS_OR_DIGITS.match(_as_string):
            raise ValueError(
                f"{_name}: version can contain only digits and dots"
                )
        major, *rest = (int(v) for v in _as_string.split('.'))
        try:
            minor = rest.pop(0)
        except IndexError:
            minor = -1
        try:
            patch = rest.pop(0)
        except IndexError:
            patch = -1
        return major, minor, patch


_RESULT_UNKNOWN = object()

class QMainThreadDispatcher(qtc.QObject):
    """Class to always execute a function in the main GUI thread.

    Instances are callable like so:
    >>> dispatcher = QMainThreadDispatcher(func, *args, **kwargs)
    >>> dispatcher(*other_args, **other_kwargs)

    Attributes
    ----------
    func : callable
        The function to be called in the main GUI thread
    args : tuple
        The positional arguments (excluding self/cls) for the
        call to func
    has_result : bool
        Whether the call to func has returned already, and
        its result is available in the .result attribute.
    kwargs : dict
        The keyword arguments for the call to func
    result : object
        The latest result of the call to func(*args, **kwargs)

    Signals
    -------
    result_arrived()
        Emitted whenever the call to func(*args, **kwargs) has been
        completed in the GUI thread.
    """

    __dispatch = qtc.pyqtSignal()
    result_arrived = qtc.pyqtSignal()

    def __init__(self, func, *args, **kwargs):
        """Create instance, move it to the main thread and call func there.

        Parameters
        ----------
        func : callable
            The function to be executed. It can be accessed (and
            freely edited) via the .func attribute.
        *args : object, optional
            Positional-only arguments for the call to func. They can
            later be accessed as a tuple via the .args attribute. If
            func is an instance or class method, do not pass self/cls
            as part of args, as this is taken care of automatically.
        call_immediately : bool, optional
            Whether func(*args, **kwargs) should be called right at
            the end of the initialization. Default is True.
        **kwargs : object, optional
            Keyword-only arguments for the call to func. They can
            later be accessed as a dictionary via the .kwargs attribute

        Returns
        -------
        None.
        """
        self.func, self.args, self.kwargs = func, args, kwargs
        self.__result = _RESULT_UNKNOWN
        super().__init__()

        if self.thread().currentThread() != qtw.qApp.thread():
            self.moveToThread(qtw.qApp.thread())

        self.__dispatch.connect(self.__call_func)
        if kwargs.pop('call_immediately', True):
            self.__dispatch.emit()

    @property
    def has_result(self):
        """Return whether the function call has already returned."""
        return self.__result is _RESULT_UNKNOWN

    @property
    def result(self):
        """Return the result of func(*args, **kwargs)."""
        return self.__result

    def __call__(self, *args, **kwargs):
        """Call func(*args, **kwargs)."""
        self.__result = _RESULT_UNKNOWN
        if args:
            self.args = args
        if kwargs:
            self.kwargs = kwargs
        self.__dispatch.emit()

    @qtc.pyqtSlot()
    def __call_func(self):
        """Call the method in the main thread, store the result."""
        assert self.thread().currentThread() == qtw.qApp.thread()
        self.__result = self.func(*self.args, **self.kwargs)
        self.result_arrived.emit()

    def moveToThread(self, new_thread):  # pylint: disable=invalid-name
        """Change thread affinity of self to new_thread."""
        self.__result = _RESULT_UNKNOWN
        super().moveToThread(new_thread)
