# -*- coding: utf-8 -*-
"""Module _reader of viperleed.tleedmlib.files.parameters.

Created on 2023-10-16

@author: Michele Riva (@michele-riva)

This module is based on the original version of readPARAMETERS from
@fkraushofer (2020). This module exists only to avoid cyclic imports
between ._read and ._write, as both need access to one of the classes
defined here (ParametersReader, RawLineParametersReader), as well as
some functionality from the other module.

Defines context-manager, iterator classes for reading parameters from
a PARAMETERS file.
"""

from contextlib import AbstractContextManager
import logging
from pathlib import Path
import re

from viperleed.tleedmlib.base import strip_comments

from .errors import ParameterNotRecognizedError
from ._known_parameters import from_alias
from ._utils import Assignment


_LOGGER = logging.getLogger('tleedm.files.parameters')


class ParametersReader(AbstractContextManager):
    """A context manager that iterates the contents of a PARAMETERS file."""

    def __init__(self, filename, noisy=True):
        """Initialize instance.

        Parameters
        ----------
        filename : str or Path
            Path to the file to be read.
        noisy : bool, optional
            Whether the reader will emit logging messages and raise
            errors if unknown or malformed parameters are encountered.
        """
        self.filename = Path(filename)
        self.noisy = noisy
        self._file_obj = None

    def __enter__(self):
        """Enter context."""
        self._file_obj = self.filename.open('r', encoding='utf-8')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close file, then return control to caller to handle exceptions."""
        try:
            self._file_obj.close()
        except AttributeError:
            pass
        return super().__exit__(exc_type, exc_value, traceback)

    def __iter__(self):
        """Yield parameters and assignments from a PARAMETERS file."""
        for line_nr, line in enumerate(self._file_obj, start=1):
            param, *rest = self._read_one_line(line, line_nr)
            if not param:
                continue
            yield (param, *rest)

    def _read_one_line(self, line, line_nr):
        """Return a parameter and other custom information from one line."""
        line = strip_comments(line)

        # We treat STOP differently, as its presence, even without
        # an equal sign, is interpreted as 'please, stop right now'
        if self._stop_requested(line):
            return 'STOP', True

        if '=' not in line:
            self._warn_about_missing_equals(line, line_nr)
            return '', None

        try:
            param, assignment = self._parse_line(line)
        except ParameterNotRecognizedError:
            if not self.noisy:
                return '', None
            raise
        if not param:
            return '', None
        return param.upper(), assignment

    def _parse_line(self, line):
        """Return a parameter string and an Assignment from line."""
        # Syntax is 'PARAMETER flag flag ... = value value value ...'
        left_side, values_str = line.split('=', maxsplit=1)
        if not left_side:  # Nothing left of '='
            return '', None
        param, *flags = left_side.split()
        param = from_alias(param)
        values_str = values_str.strip()
        if not values_str:
            raise ParameterNotRecognizedError(parameter=param)
        assignment = Assignment(values_str=values_str,
                                parameter=param,
                                flags_str=flags)
        return param, assignment

    def _stop_requested(self, line):
        """Return whether line contains a STOP condition."""
        # SEARCH_KILL is legacy name
        line = line.upper()
        for param in ['SEARCH_KILL', 'STOP']:
            if (line.startswith(param)
                    and not re.match(fr'\s*{param}\s*=\s*[F](ALSE)?', line)):
                return True
        return False

    def _warn_about_missing_equals(self, line, line_nr):
        """Warn the user if line contains a known parameter but no '='."""
        if not line or not self.noisy:
            return
        param, *_ = line.split()
        try:
            param = from_alias(param)
        except ParameterNotRecognizedError:
            return
        _LOGGER.warning(f'PARAMETERS file, line {line_nr}: found {param} '      # TODO: should we actually raise?
                        'in a line without an "=" sign. Assignment will be '
                        f'SKIPPED.\n    Faulty line: {line}')

