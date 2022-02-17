"""Module settings of viperleed.guilib.measure.classes.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2022-02-16
Author: Michele Riva

This module defines the ViPErLEEDSettings class, a ConfigParser
subclass used to read/write and handle configurations for hardware
equipment and measurements.
"""

from configparser import (ConfigParser, MissingSectionHeaderError,
                          DuplicateSectionError, DuplicateOptionError,
                          SectionProxy)
from collections import defaultdict
import inspect
import os
from pathlib import Path
import sys

from viperleed import guilib as gl

_MEASURE_PATH = Path(__file__).parent.parent
SYSTEM_CONFIG_PATH = _MEASURE_PATH / "_defaults" / "_system_settings.ini"


def get_system_config():
    """Return a ConfigParser loaded with system settings."""
    config = ConfigParser()
    config.read(SYSTEM_CONFIG_PATH)
    return config

def _interpolate_config_path(filenames):
    """Interpolate filenames with system settings.

    Replaces "__CONFIG__" at the baginning of filenames
    with the path contained in the system-wide settings,
    unless (i) there is no valid system-wide setting, or
    (ii) a filename contains "__CONFIG__" multiple times.

    Parameters
    ----------
    filenames : Sequence
        File names to be interpolated. The interpolation is
        done in-place.
    """
    cfg = get_system_config()
    _sys_path = cfg.get('settings', 'configuration_path', fallback=None)
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


class MissingSettingsFileError(Exception):
    """Exception raised when failed to read settings file(s)."""


class NoSettingsError(Exception):
    """Exception raised when failed to read settings file(s)."""


class NoDefaultSettingsError(Exception):
    """Exception raised when no default settings file was found."""


class ViPErLEEDSettings(ConfigParser):
    """Class for read/write and handling ViPErLEED settings.

    The class rememebers the path to the last file read, which
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
        the parser will try to replace an intial "__CONFIG__"
        in the path(s) to be read with the default path of the
        configuration files found in the system-wide settings.
        Defaults to True
    """

    def __init__(self, *args, **kwargs):
        """Init instance with custom defaults."""
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
        self.__last_file = ""

    @classmethod
    def from_settings(cls, settings, find_from=None):
        """Return a ViPErLEEDSettings from the settings passed.

        Parameters
        ----------
        settings : str or os.PathLike or ViPErLEEDSettings or None
            The settings to load from. If settings is None,
            find_from will be used to search for a default settings.
            If a ViPErLEEDSettings, no copy is made.
        find_from : str or None, optional
            The string to look for in a configuration file to
            be loaded and returned. If None, no search will be
            performed.

        Returns
        -------
        loaded_settings : ViPErLEEDSettings
            The ViPErLEEDSettings after loading settings.

        Raises
        ------
        ValueError
            If settings is not a ViPErLEEDSettings, it is
            False-y and find_from is False-y.
        NoSettingsError
            If settings is invalid and find_from was not given.
        NoDefaultSettingsError
            If settings is invalid, and looking for default settings
            failed.
        """
        if settings is None and not find_from:
            raise ValueError(f"{cls.__name__}: cannot create from nothing.")
        if isinstance(settings, cls):
            return settings

        if not settings and not find_from:
            raise ValueError(f"{cls.__name__}: cannot create from nothing.")

        config = cls()
        if isinstance(settings, (dict, ConfigParser)):
            config.read_dict(settings)
            return config

        try:
            config.read(settings)
        except (TypeError, MissingSettingsFileError):
            pass
        else:
            return config

        if not find_from:
            raise NoSettingsError(
                f"{cls.__name__}: could not load settings, and "
                "no find_from was provided to look for a default."
                )

        # Failed to read from settings. Try with find_from.
        get_cfg = gl.measure.hardwarebase.get_device_config
        settings = get_cfg(find_from, prompt_if_invalid=False)
        if settings:
            try:
                config.read(settings)
            except MissingSettingsFileError:
                pass
            else:
                return config
        raise NoDefaultSettingsError(
            f"{cls.__name__}: could not find (or load) a suitable default."
            )

    def __bool__(self):
        """Return True if there is any section."""
        return bool(self.sections())

    @property
    def comments(self):
        """Return all comments."""
        return self.__comments

    @property
    def last_file(self):
        """Return the path to the last file read."""
        return self.__last_file

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

        Raises
        ------
        TypeError
            If required_settings contains invalid data (i.e., one
            of the entries is not a Sequence with length <= 3).
        """
        if not self:
            return required_settings

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
                invalid_settings.append("/".join(setting))
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
        the file name of the last succesfully read-in file.
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
            _interpolate_config_path(filenames)

        read_ok = super().read(filenames, encoding=encoding)
        if read_ok:
            self.__last_file = Path(read_ok[-1])

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

    def read_file(self, f, source=None):
        """Read from a file-like object."""
        fname = ""
        try:
            fname = f.name
        except AttributeError:
            pass

        if fname:
            self.__last_file = Path(fname)
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
            self.__last_file = Path(fp.name)
        except AttributeError:
            pass

        super().write(fp, space_around_delimiters)

    # pylint: disable=too-complex,too-many-locals
    # pylint: disable=too-many-branches,too-many-statements
    # The code below is essentially a (simplified) copy of the
    # one in the configparser standard library for RawConfigParser.
    def _read(self, fp, fpname):
        """Reimplement to preserve comments."""
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
        """Write a single section to the specified `fp'."""
        # Header
        fp.write(f"[{section_name}]\n")

        # All comments
        if self.__comments[section_name]:
            fp.write("\n".join(self.__comments[section_name]) + "\n")

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
            The line to be checked an stored
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
