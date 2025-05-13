"""Module settings of viperleed.guilib.measure.classes.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2022-02-16
Author: Michele Riva

This module defines the ViPErLEEDSettings class, a ConfigParser
subclass used to read/write and handle configurations for hardware
equipment and measurements. It also defines the SystemSettings
subclass of ViPErLEEDSettings, a singleton providing thread-safe
access to system settings.
"""

import ast
from configparser import (ConfigParser, MissingSectionHeaderError,
                          DuplicateSectionError, DuplicateOptionError,
                          SectionProxy, _UNSET)
from collections import defaultdict
from collections.abc import Sequence
import copy
import os
from pathlib import Path
import sys

from wrapt import synchronized  # thread-safety decorator

from viperleed.guilib.measure.dialogs.settingsdialog import SettingsHandler
from viperleed.guilib.measure.widgets.pathselector import PathSelector


# TODO: not nice. Also, there's two places where the _defaults
# path is used. Here and in hardwarebase. However, due to circular
# imports I cannot use either for the other one. Probably better
# to move the path in __init__?
_MEASURE_PATH = Path(__file__).parent.parent
SYSTEM_CONFIG_PATH = _MEASURE_PATH / "_defaults" / "_system_settings.ini"


def interpolate_config_path(filenames):
    """Interpolate filenames with system settings.

    Replaces "__CONFIG__" at the beginning of filenames
    with the path contained in the system-wide settings,
    unless (i) there is no valid system-wide setting, or
    (ii) a filename contains "__CONFIG__" multiple times.

    Parameters
    ----------
    filenames : Sequence
        File names to be interpolated. The interpolation is
        done in-place.
    """
    # NB: the next lines should NOT use SystemSettings, as
    # interpolate_config_path is used in ViPErLEEDSettings,
    # the ancestor of SystemSettings. This leads to infinite
    # recursion.
    cfg = ConfigParser()
    cfg.read(SYSTEM_CONFIG_PATH)
    _sys_path = cfg.get('PATHS', 'configuration', fallback=None)
    if not _sys_path:
        return

    for i, fname in enumerate(filenames):
        if isinstance(fname, os.PathLike):
            fname = os.fspath(fname)
        if not isinstance(fname, (str, bytes)):
            continue
        if (fname.startswith("__CONFIG__/")
                and fname.count("__CONFIG__") == 1):
            filenames[i] = fname.replace("__CONFIG__", _sys_path)


class SettingsError(Exception):
    """Base exception for all settings-related errors."""


class DefaultSettingsError(SettingsError):
    """Exception raised when the default settings are corrupted."""


class MissingSettingsFileError(SettingsError):
    """Exception raised when failed to read settings file(s)."""


class NoSettingsError(SettingsError):
    """Exception raised when failed to read settings file(s)."""


class NoDefaultSettingsError(DefaultSettingsError):
    """Exception raised when no default settings file was found."""


class NotASequenceError(SettingsError):
    """Exception raised when getsequence fails to return a sequence."""


class TooManyDefaultSettingsError(DefaultSettingsError):
    """Exception raised when too many default settings files were found."""


class ViPErLEEDSettings(ConfigParser):
    """Class for read/write and handling ViPErLEED settings.

    The class remembers the path to the last file read, which
    can be updated by calling .update_file(), and remembers all
    the comments, which are rewritten to file right below each
    section header.

    The initialization arguments are changed as follows relative
    to ConfigParser:
    * allow_no_value : fixed to True
    * comment_prefixes : fixed to "#" and ";"
    * inline_comment_prefixes : fixed to None
    * empty_lines_in_values : fixed to False
    * strict : fixed to True
    * interpolate_paths : bool, optional
        New keyword-only argument used while reading. If True,
        the parser will try to replace an initial "__CONFIG__"
        in the path(s) to be read with the default path of the
        configuration files found in the system-wide settings.
        Defaults to True
    """

    def __init__(self, *args, **kwargs):
        """Initialize instance with custom defaults."""
        # ConfigParser has 3 positional-or-keyword arguments. We
        # want to set the third (allow_no_value) always to True.
        args = args[:2]  # skip positional allow_no_value
        kwargs["allow_no_value"] = True

        kwargs["comment_prefixes"] = "#;"
        kwargs["inline_comment_prefixes"] = ""  # no in-line comments

        # Disallow empty lines in values, as this has the potential
        # to create mess; also make the parser strict (no duplicates)
        kwargs["empty_lines_in_values"] = False
        kwargs["strict"] = True

        self.__interpolate_paths = kwargs.pop("interpolate_paths", True)
        super().__init__(*args, **kwargs)

        self.__comments = defaultdict(list)

        # Keep track of the last .ini file loaded, if read from
        # a file on disk, and the path to the "directory" containing
        # it. This is necessary to resolve special "./" paths in
        # the config file. __base_dir should be set from the outside
        # using the .base_dir property to the "<path/to/archive.zip>"
        # when the content of the config was read from an archive.
        self.__last_file = ""
        self.__base_dir = ""

    def __str__(self):
        """Return a simple string representation of self."""
        return str(self.as_dict())

    @classmethod
    def from_settings(cls, settings):
        """Return a ViPErLEEDSettings from the settings passed.

        Parameters
        ----------
        settings : dict or ConfigParser or str or Path or ViPErLEEDSettings
            The settings to load from.
            If a ViPErLEEDSettings, no copy is made.

        Returns
        -------
        loaded_settings : ViPErLEEDSettings
            The ViPErLEEDSettings after loading settings.

        Raises
        ------
        ValueError
            If there are no settings to create a ViPErLEEDSettings from.
        NoSettingsError
            If settings is invalid.
        """
        if not settings:
            raise ValueError(f'{cls.__name__}: cannot create from nothing.')

        if isinstance(settings, cls):
            return settings

        config = cls()
        if isinstance(settings, (dict, ConfigParser)):
            config.read_dict(settings)
            return config

        try:
            config.read(settings)
        except (TypeError, MissingSettingsFileError) as exc:
            raise NoSettingsError(f'{cls.__name__}: could not '
                                  'load settings.') from exc
        return config

    def __bool__(self):
        """Return True if there is any section."""
        return bool(self.sections())

    @property
    def base_dir(self):
        """Return the directory from which the config was read."""
        return self.__base_dir

    @base_dir.setter
    def base_dir(self, new_dir):
        """Set the directory from which the config was read."""
        self.__base_dir = str(new_dir)

    @property
    def comments(self):
        """Return all comments."""
        return self.__comments

    @property
    def last_file(self):
        """Return the path to the last file read."""
        return self.__last_file

    def as_dict(self):
        """Return a dictionary version of self."""
        return {sec_name: dict(sec) for sec_name, sec in self.items()}

    def getsequence(self, section, option, *, fallback=_UNSET, **kwargs):
        """Return a section/option as a sequence."""
        conv = ast.literal_eval
        try:
            converted = self._get_conv(section, option, conv,
                                       fallback=fallback, **kwargs)
        except (ValueError, TypeError, SyntaxError,
                MemoryError, RecursionError) as err:
            raise NotASequenceError(f"Could not convert {section}/{option} "
                                    "into a sequence") from err
        if converted is not fallback and not isinstance(converted, Sequence):
            raise NotASequenceError(f"Could not convert {section}/{option} "
                                    "into a sequence")
        return converted

    def has_settings(self, *required_settings):
        """Check whether self has given settings.

        Parameters
        ----------
        *required_settings : tuple
            Each required setting should be a tuple of the forms
            (<section>, )
                Check that config contains the section <section>.
                No options are checked.
            (<section>, <option>)
                Check that the config contains <option> in the
                existing section <section>. No values are checked
                for <option>.
            (<section>, <option>, <admissible_values>)
                Check that the config contains <option> in the
                existing section <section>, and that <option> is
                one of the <admissible_values>. In this case,
                <admissible_values> is a Sequence of strings.

        Returns
        -------
        invalid_settings : list
            Invalid required_settings of self as a list of strings.
            Each entry can be either '<section>', '<section>/<option>',
            or '<section>/<option> not one of <value1>, <value2>, ...'

        Raises
        ------
        TypeError
            If required_settings contains invalid data (i.e., one
            of the entries is not a Sequence with length <= 3).
        """
        invalid_settings = []
        for setting in required_settings:
            if not setting or len(setting) > 3:
                raise TypeError(f"Invalid mandatory setting {setting}. "
                                f"with length {len(setting)}. Expected "
                                "length <= 3.")
            # (<section>,)
            if len(setting) == 1:
                if not self.has_section(setting[0]):
                    invalid_settings.append(setting[0])
                continue

            # (<section>, <option>) or (<section>, <option>, <admissible>)
            section, option = setting[:2]
            if not self.has_option(section, option):
                invalid_settings.append(f"{section}/{option}")
                continue

            # (<section>, <option>, <admissible>)
            if len(setting) == 3:
                admissible_values = setting[2]
                value = self[section][option]
                if value not in admissible_values:
                    invalid_settings.append(
                        "/".join(setting[:2])
                        + " not one of " + ", ".join(admissible_values)
                        )
        return invalid_settings

    def read(self, filenames, encoding=None):
        """Read and parse a filename or an iterable of filenames.

        This extension of the ConfigParser implementation stores
        the file name of the last successfully read-in file.
        Differently from the base implementation, it also complains
        if any of the files to be read was not.

        Parameters
        ----------
        filenames : str or PathLike or Sequence
            Path(s) to the files to be read. Should any of the
            filenames start with the special "__CONFIG__/" string,
            this will be replaced with the system-wide standard
            location for configuration files.
        encoding : str or None, optional
            Encoding to be used to open the files. Default is None.

        Returns
        -------
        read_ok : list
            List of successfully read files.

        Raises
        ------
        MissingSettingsFileError
            If any of the files in filename was not read successfully.
        """
        if isinstance(filenames, (str, bytes, os.PathLike)):
            filenames = [filenames]
        if self.__interpolate_paths:
            interpolate_config_path(filenames)

        read_ok = super().read(filenames, encoding=encoding)
        if read_ok:
            self.__last_file = Path(read_ok[-1]).resolve()
            self.base_dir = self.__last_file.parent

        if len(read_ok) != len(filenames):
            # Need to run through again, as super().read
            # converts path-like entries to string before
            # inserting them in read_ok
            missing = []
            for fname in filenames:
                if isinstance(fname, os.PathLike):
                    fname = os.fspath(fname)
                if fname not in read_ok:
                    missing.append(fname)
            raise MissingSettingsFileError(
                "Could not read the following settings file(s):"
                ', '.join(missing)
                )
        return read_ok

    def read_again(self):
        """Read once more the last file, returning True if successful."""
        if not self.last_file or not self.last_file.exists():
            return False
        try:
            self.read(self.last_file)
        except MissingSettingsFileError:
            return False
        return True

    def read_file(self, f, source=None):
        """Read from a file-like object."""
        fname = ""
        try:
            fname = f.name
        except AttributeError:
            pass

        if fname:
            self.__last_file = Path(fname).resolve()
            self.base_dir = self.__last_file.parent
        super().read_file(f, source)

    def update_file(self):
        """Update the last file read with the data in self."""
        if not self.__last_file:
            raise RuntimeError(f"{self.__class__.__name__}: no known "
                               "last file read to be updated.")
        with open(self.__last_file, 'w', encoding='utf-8') as fproxy:
            self.write(fproxy)

    def write(self, fp, space_around_delimiters=True):
        """Write an .ini-format representation of the configuration state.

        Parameters
        ----------
        fp : io.TextIOBase
            The open file-like object to write to. Should be
            open in "write-text" mode.
        space_around_delimiters : bool
            If `space_around_delimiters' is True (the default),
            delimiters between keys and values are surrounded
            by spaces.
        """
        # Start with comments at the beginning
        fp.write("\n".join(self.__comments[None]))
        if self.__comments[None]:
            # Add newline before the first section
            # header if there were any comments
            fp.write("\n")

        # And store this file's name
        try:
            self.__last_file = Path(fp.name).resolve()
        except AttributeError:
            pass
        else:
            self.base_dir = self.__last_file.parent

        super().write(fp, space_around_delimiters)

    # pylint: disable=too-complex,too-many-locals
    # pylint: disable=too-many-branches,too-many-statements
    # The code below is essentially a (simplified) copy of the
    # one in the configparser standard library for RawConfigParser.
    def _read(self, fp, fpname):
        """Override original to preserve comments."""
        elements_added = set()
        cursect = None                        # None, or a dictionary
        sectname = None
        optname = None
        lineno = 0
        indent_level = 0
        err = None                            # None, or an exception
        for lineno, line in enumerate(fp, start=1):
            # Keep track of full-line comments, but prevent creating
            # duplicates when reading two files with the same sections
            # and the same comments.
            if self.__store_if_comment(line, sectname):
                continue

            value = line.strip()
            if not value:
                # empty line marks end of value
                indent_level = sys.maxsize
                continue

            # continuation line?
            first_nonspace = self.NONSPACECRE.search(line)
            cur_indent_level = first_nonspace.start() if first_nonspace else 0
            if (cursect is not None and optname and
                    cur_indent_level > indent_level):
                # pylint: disable=unsubscriptable-object
                # bug? cursect is subscriptable if we enter this
                cursect[optname].append(value)
                continue

            # a section header or option header?
            indent_level = cur_indent_level
            _match = self.SECTCRE.match(value)

            # no section header in the file?
            if not _match and cursect is None:
                raise MissingSectionHeaderError(fpname, lineno, line)

            if _match:  # is a section header
                sectname = _match.group('header')
                if sectname in self._sections:
                    if sectname in elements_added:  # strict by default
                        raise DuplicateSectionError(sectname, fpname, lineno)
                    cursect = self._sections[sectname]
                    elements_added.add(sectname)
                elif sectname == self.default_section:
                    cursect = self._defaults
                else:
                    cursect = self._dict()
                    self._sections[sectname] = cursect
                    self._proxies[sectname] = SectionProxy(self, sectname)
                    elements_added.add(sectname)
                # So sections can't start with a continuation line
                optname = None
                continue

            # an option line?
            _match = self._optcre.match(value)
            if _match:
                optname, optval = _match.group('option', 'value')
                if not optname:
                    err = self._handle_error(err, fpname, lineno, line)
                optname = self.optionxform(optname.rstrip())
                if (sectname, optname) in elements_added:
                    raise DuplicateOptionError(sectname, optname,
                                               fpname, lineno)
                elements_added.add((sectname, optname))
                # This check is fine because the OPTCRE cannot
                # match if it would set optval to None
                if optval is not None:
                    optval = optval.strip()
                    cursect[optname] = [optval]
                else:
                    # valueless option handling
                    cursect[optname] = None
                continue

            # a non-fatal parsing error occurred. set up the
            # exception but keep going. the exception will be
            # raised at the end of the file and will contain a
            # list of all bogus lines
            err = self._handle_error(err, fpname, lineno, line)

        self._join_multiline_values()
        # if any parsing errors occurred, raise an exception
        if err:
            raise err
    # pylint: enable=too-complex,too-many-locals
    # pylint: enable=too-many-branches,too-many-statements

    def _write_section(self, fp, section_name, section_items, delimiter):
        """Write a single section to the specified 'fp'."""
        # Header
        fp.write(f"[{section_name}]\n")

        # All comments (newlines were not stripped)
        if self.__comments[section_name]:
            fp.write("".join(self.__comments[section_name]))

        # And all the options/values
        for key, value in section_items:
            value = self._interpolation.before_write(self, section_name,
                                                     key, value)
            if value is not None or not self._allow_no_value:
                value = delimiter + str(value).replace('\n', '\n\t')
            else:
                value = ""
            fp.write(f"{key}{value}\n")
        fp.write("\n")

    def __store_if_comment(self, line, sectname):
        """Store a line as comment if it is one.

        Used while reading from file. Comment lines are stored if
        they are not duplicates.

        Parameters
        ----------
        line : str
            The line to be checked and stored
        sectname : str or None
            The name of the section where the comment appeared.
            None if the line appears before any section header.

        Returns
        -------
        is_comment : bool
            True if the line was a comment (whether it was stored
            or not).
        """
        for prefix in self._comment_prefixes:
            if line.strip().startswith(prefix):
                if line not in self.__comments[sectname]:
                    self.__comments[sectname].append(line)
                return True
        return False


@synchronized
class SystemSettings(ViPErLEEDSettings):
    """A thread-safe singleton for accessing/editing system settings.

    Deep-copy is possible, but will return a ViPErLEEDSettings instance
    with a copy of all the system settings at the time of copying.
    """

    __instance = None
    __non_null = (
        ('PATHS', 'configuration'),
        ('PATHS', 'measurements'),
        )

    # __mandatory can later be extended as __non_null + (...)
    # Can also decide to add depending on the Version of viperleed
    __mandatory = __non_null

    # __non_mandatory entries are settings that can be added, but
    # that are not enforced. They should be used if the setting is
    # only required by a part of the package that is not essential
    # for it's core functionality.
    __non_mandatory = (
        ('PATHS', 'arduino_cli'),
        ('PATHS', 'firmware'),
        )

    def __new__(cls, *args, **kwargs):
        """Return an uninitialized instance of cls."""
        if cls.__instance is None:
            new_instance = super().__new__(cls, *args, **kwargs)
            new_instance._initialized = False
            cls.__instance = new_instance
        return cls.__instance

    def __deepcopy__(self, memo):
        """Return a ViPErLEEDSettings with the same contents as self."""
        fake_copy = ViPErLEEDSettings()
        fake_copy.read_dict(self)
        return copy.deepcopy(fake_copy, memo)

    def __init__(self, *args, **kwargs):
        """Initialize instance once."""
        if self._initialized:  # pylint: disable=no-member
            return
        super().__init__(*args, **kwargs)
        self.__handler = None

        self.read(SYSTEM_CONFIG_PATH)
        self.__check_mandatory_settings()

    @property
    def paths(self):
        """Return the 'PATHS' section of self."""
        return self['PATHS']

    @property
    def valid(self):
        """Return whether all settings are valid."""
        return all(self[sec][opt] for sec, opt in self.__non_null)

    @property
    def settings(self):
        """Return self."""
        return self

    def get_settings_handler(self):
        """Return a SettingsHandler instance for displaying self."""
        if self.__handler is not None:
            return self.__handler

        handler = self.__handler = SettingsHandler(self)
        handler.make_from_config()

        # Add some informative text to the entries
        _infos = (
        ('PATHS', 'configuration',
         '<nobr>This is the directory that contains all the</nobr> '
         'configuration files for your devices and measurements. '
         'It must be set before you can run any measurement.'),
        ('PATHS', 'measurements',
         '<nobr>This is the default folder where all your</nobr> '
         'measurements will be automatically saved. IN THE FUTURE '
         'you will be able to decide if you want to be asked each '
         'time a measurement starts.'),
        ('PATHS', 'arduino_cli',
         '<nobr>This is the folder in which the Arduino</nobr> '
         'command-line interface is installed.'),
        ('PATHS', 'firmware',
         '<nobr>This is the folder containing archives with</nobr> '
         'firmware for ViPErLEED controllers.'),

        )
        for section, option, info in _infos:
            handler[section][option].set_info_text(info)

        # Make sure the PATHS are "PathSelector"s even with empty fields
        for option in handler['PATHS'].values():
            if not isinstance(option.handler_widget, PathSelector):
                option.handler_widget = PathSelector(select_file=False)

        # Make it clear which settings are mandatory by
        # adding a '*' at the beginning of the label text
        for section, option in self.__mandatory:
            label = handler[section][option].label
            unstarred_name = label.text()
            if unstarred_name.startswith('* '):
                continue
            label.setText('* ' + unstarred_name)

        return self.__handler

    def __check_mandatory_settings(self):
        """Check, and possibly add missing settings.

        Missing, non-null settings are added as empty strings.

        Raises
        ------
        RuntimeError
            If any mandatory setting is missing (after potentially
            filling the non-null ones)
        """
        invalid = self.has_settings(*self.__non_null, *self.__non_mandatory)

        # We're missing settings. Let's add them back...
        for missing in invalid:
            section, option = missing.split('/')
            if section not in self:
                self.add_section(section)
            self.set(section, option, '')

        if invalid:
            # ...and save changes
            if self.last_file:
                self.update_file()
            else:
                # We even do not have a system settings file
                self.write(SYSTEM_CONFIG_PATH)

        # Now check all those that must be there
        missing = self.has_settings(*self.__mandatory)
        if missing:
            raise RuntimeError(
                f"System settings file at {SYSTEM_CONFIG_PATH} is "
                "missing the following mandatory sections/options: "
                "; ".join(missing)
                )
