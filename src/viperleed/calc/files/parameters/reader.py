"""Module reader of viperleed.calc.files.parameters.

This module is based on the original version of readPARAMETERS from
@fkraushofer (2020). This module exists only to avoid cyclic imports
between .read and .write, as both need access to one of the classes
defined here (ParametersReader, RawLineParametersReader), as well as
some functionality from the other module.

Defines context-manager, iterator classes for reading parameters from
a PARAMETERS file.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-16'
__license__ = 'GPLv3+'

import re

from viperleed.calc.lib.string_utils import strip_comments

from .errors import MissingEqualsError
from .errors import ParameterHasNoValueError
from .errors import ParameterNotRecognizedError
from ..file_reader import InputFileReader
from .known_parameters import did_you_mean
from .known_parameters import from_alias
from .utils import Assignment

class ParametersReader(InputFileReader):
    """A context manager that iterates the contents of a PARAMETERS file."""

    def __next__(self):
        """Return the next understandable information in the file."""
        for line in self._file_obj:
            self._current_line += 1
            param, *rest = self._read_one_line(line)
            if not param:
                continue
            return (param, *rest)
        raise StopIteration

    def _complain_about_line_parse_errors(self, line, exc):
        """Re-raise an exception occurred while line was being parsed."""
        if not isinstance(exc, ParameterNotRecognizedError):
            raise exc
        # Suggest a close match for an unknown parameter
        param, *_ = self._tokenize_line(line)
        try:
            close_match = did_you_mean(param)
        except ParameterNotRecognizedError:
            # No close match. Stick to the original exception
            raise exc from None
        new_message = (f'{exc.message} (at line {self._current_line}). '
                       f'Did you mean {close_match!r}?')
        raise ParameterNotRecognizedError(exc.parameter, message=new_message)

    def _complain_about_missing_equals(self, line):
        """Warn the user if line contains a known parameter but no '='."""
        if not line or not self.noisy:
            return
        param, *_ = line.split()
        try:
            param = from_alias(param)
        except ParameterNotRecognizedError:
            return
        _err_msg = (f'-- line {self._current_line} -- Found {param} '
                    'in a line without an "=" sign. Assignment will '
                    f'be SKIPPED.\n    Faulty line: {line}')
        raise MissingEqualsError(param, message=_err_msg)

    def _read_one_line(self, line):
        """Return a parameter and other custom information from one line."""
        line = strip_comments(line)

        # We treat STOP differently, as its presence, even without
        # an equal sign, is interpreted as 'please, stop right now'
        if self._stop_requested(line):
            return 'STOP', Assignment('True', 'STOP')

        if '=' not in line:
            self._complain_about_missing_equals(line)
            return '', None

        try:
            param, assignment = self._parse_line(line)
        except (ParameterNotRecognizedError, ParameterHasNoValueError) as exc:
            if not self.noisy:
                return '', None
            self._complain_about_line_parse_errors(line, exc)
        if not param:
            return '', None
        return param.upper(), assignment

    def _parse_line(self, line):
        """Return a parameter string and an Assignment from line."""
        param, flags, values_str = self._tokenize_line(line)
        if not param:
            return '', None
        param = from_alias(param)
        values_str = values_str.strip()
        if not values_str:
            raise ParameterHasNoValueError(parameter=param)
        assignment = Assignment(values_str=values_str,
                                parameter=param,
                                flags_str=flags)
        return param, assignment

    def _stop_requested(self, line):
        """Return whether line contains a STOP condition."""
        # SEARCH_KILL is legacy name
        line = line.upper()
        stop_false = {alias: re.compile(fr'\s*{alias}\s*=\s*[F](ALSE)?').match
                      for alias in ('SEARCH_KILL', 'STOP')}
        return any(line.startswith(stop) and not is_false(line)
                   for stop, is_false in stop_false.items())

    @staticmethod
    def _tokenize_line(line):
        """Split up line into tokens."""
        # Syntax is 'PARAMETER flag flag ... = value value value ...'
        left_side, values_str = line.split('=', maxsplit=1)
        if not left_side:  # Nothing left of '='
            return left_side, (), values_str
        param, *flags = left_side.split()
        return param, flags, values_str


# To me this pylint complaint does not make much sense
# here. The public methods come from the parent.
# pylint: disable-next=too-few-public-methods
class RawLineParametersReader(ParametersReader):
    """A ParametersReader that also returns lines exactly as they were read.

    Differently from a ParametersReader, a RawLineParametersReader
    returns all the lines, whether they contain a parameter or not.
    """

    def __next__(self):
        """Return the next acceptable information in the file."""
        self._current_line += 1
        return self._read_one_line(next(self._file_obj))

    def _read_one_line(self, line):
        """Return a parameter, and the whole raw line it was in."""
        param, *_ = super()._read_one_line(line)
        return param, line
