"""Version class to keep track of version numbers internally.

Taken originally from guilib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-04-15'
__license__ = 'GPLv3+'

import re

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
        self.__has_parts = (True, minor >= 0, patch >= 0)
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

    def __from_string(self, _as_string):                                        # TODO: complain if _as_string has more then three dotted fields
        """Return major, minor, patch from a string."""
        _name = self.__class__.__name__
        if not _as_string:
            # Empty version
            raise ValueError(f"{_name}: version cannot be an empty string")
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
