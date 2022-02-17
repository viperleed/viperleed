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
                          DuplicateSectionError, DuplicateOptionError)
from collections import defaultdict
import io
import os
from pathlib import Path

from viperleed import guilib as gl


SYSTEM_CONFIG_PATH = (Path(inspect.getfile(gl.measure)).parent
                      / "_defaults" / "_system_settings.ini")


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
        fname = str(fname)
        if (fname.startswith("__CONFIG__/")
                and fname.count("__CONFIG__") == 1):
            filenames[i] = fname.replace("__CONFIG__", _sys_path)


class MissingSettingsFileError(Exception):


class ViPErLEEDSettings(ConfigParser):
    """Class for read/write and handling ViPErLEED settings.

    The class rememebers the path to the last file read, which
    can be updated by calling .update_file(), and remembers all
    the comments, which are rewritten to file right below each
    section header.
    """

    def __init__(self, *args, **kwargs):
        """Initalize instance with custom defaults."""
        # ConfigParser has 3 positional-or-keyword arguments. We
        # want to set the third (allow_no_value) always to False.
        args = args[:2]  # skip positional allow_no_value
        kwargs["allow_no_value"] = False

        kwargs["comment_prefixes"] = "#;"
        kwargs["inline_comment_prefixes"] = ""  # no in-line comments

        # Disallow empty lines in values, as this has
        # the potential to create mess; and make it strict
        kwargs["empty_lines_in_values"] = False
        kwargs["strict"] = True

        super().__init__(*args, **kwargs)

        self.__comments = defaultdict(list)
        self.__last_file = ""

    def __bool__(self):
        """Return True if there is any section."""
        return bool(self.sections)

    def has_settings(self, required_settings):
        """Check whether self has given settings.

        Parameters
        ----------
        required_settings : Sequence
            Each element is a tuple of the forms
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

    def read(self, filenames, encoding=None, interpolate=True):
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
        if interpolate:
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
        """Read from a file-like object.

        Parameters
        ----------
        f : file-like
            The `f' argument must be iterable, returning one line
            at a time.
        source : str or None, optional
            The name of the file being read. If None, it is taken
            from f.name. If `f' has no `name' attribute, `<???>' is
            used. Default is None.

        Returns
        -------
        None.
        """
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
        with open(self.__last_file, 'w') as fproxy:
            self.write(fproxy)

    def write(self, fp, space_around_delimiters=True):
        """Write an .ini-format representation of the configuration state.

        Parameters
        ----------
        fp : file-like
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

    def _read(self, fileproxy, fpname):
        """Reimplement to preserve comments."""
        elements_added = set()
        cursect = None                        # None, or a dictionary
        sectname = None
        optname = None
        lineno = 0
        indent_level = 0
        err = None                            # None, or an exception
        for lineno, line in enumerate(fileproxy, start=1):
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
                optname, vi, optval = _match.group('option', 'vi', 'value')
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
