"""Module file_reader of viperleed.calc.classes.

This module is a more abstract version of a settings file reader. It is based on
@michele-riva's work on the PARAMETERS file reader, which in turn is based on
the old readPARAMETERS from @fkraushofer.

Defines context-manager, iterator classes for reading parameters from a file.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-03'
__license__ = 'GPLv3+'

from abc import abstractmethod
from collections.abc import Iterator
from contextlib import AbstractContextManager
from pathlib import Path
import re

from viperleed.calc.lib.string_utils import strip_comments

from .errors import MissingEqualsError
from .errors import ParameterHasNoValueError
from .errors import ParameterNotRecognizedError
from .known_parameters import did_you_mean
from .known_parameters import from_alias
from .utils import Assignment

class SettingsFileReader(AbstractContextManager, Iterator):
    """A context manager that iterates the contents of a settings file
    (e.g PARAMETERS)."""

    def __init__(self, filename, noisy=True):                                   # TODO: it would be nice to support passing file contents via a StringIO or similar
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
        self._current_line = 0

    def __enter__(self):
        """Enter context."""
        self._file_obj = self.filename.open('r', encoding='utf-8')
        self._current_line = 0
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close file, then return control to caller to handle exceptions."""
        try:
            self._file_obj.close()
        except AttributeError:
            pass
        return super().__exit__(exc_type, exc_value, traceback)

    def __next__(self):
        """Return the next understandable information in the file."""
        for line in self._file_obj:
            self._current_line += 1
            param, *rest = self._read_one_line(line)
            if not param:
                continue
            return (param, *rest)
        raise StopIteration

    @abstractmethod
    def _read_one_line(self, line):
        """Return a parameter and other custom information from one line."""
        pass

    @abstractmethod
    def _parse_line(self, line):
        """Return a parameter string and an Assignment from line."""
        pass

    @staticmethod
    @abstractmethod
    def _tokenize_line(line):
        """Split up line into tokens."""
        pass
