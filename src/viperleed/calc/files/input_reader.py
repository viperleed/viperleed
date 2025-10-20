"""Module input_reader of viperleed.calc.files.

This module is a more abstract version of an input-file reader. It is based on
@michele-riva's work on the PARAMETERS file reader, which in turn is based on
the old readPARAMETERS from @fkraushofer.

Defines context-manager, iterator classes for reading information from a file
as well as from TextIOBase streams.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-03'
__license__ = 'GPLv3+'

from abc import ABC
from abc import abstractmethod
from io import TextIOBase
from collections.abc import Iterator
from contextlib import AbstractContextManager
import logging
from pathlib import Path

from viperleed.calc.constants import LOG_VERY_VERBOSE

_LOGGER = logging.getLogger(__name__)

# TODO: Add subclasses for the InputReaders as suggested by @michele-riva in
#       PR #266. E.g. DelimitedFileReader, BlockFileReader, etc.

class ShouldSkipLineError(Exception):
    """Exception raised when a line's contents should be skipped."""


class InputReader(Iterator):
    """Common base class for all input readers."""

    def __init__(self, noisy=True):
        """Initialize base-class instance.

        Parameters
        ----------
        noisy : bool, optional
            Whether the reader will emit logging messages and raise
            errors if unknown or malformed lines are encountered.

        Returns
        -------
        None.
        """
        self.noisy = noisy
        self._current_line = 0
        super().__init__()

    def __next__(self):
        """Return the next understandable information in the file."""
        for line in self.stream:
            self._current_line += 1
            try:
                return self._read_one_line(line)
            except ShouldSkipLineError:
                if self.noisy:
                    line = line.rstrip('\n')
                    _LOGGER.log(LOG_VERY_VERBOSE,
                                f'Skipping line {self._current_line} '
                                 f'in input file: {line!r}.')
                continue
        raise StopIteration

    @property
    @abstractmethod
    def stream(self):
        """Return the input stream."""

    @abstractmethod
    def _read_one_line(self, line):
        """Return understandable information from `line`.

        This method is guaranteed to be called once on each
        line read from self.stream.

        Parameters
        -------------
        line : str
            A single line read from `self.stream.`

        Returns
        --------
        info : object
            Understandable information read from `line`. This is the
            same object returned whenever this reader is iterated over.

        Raises
        -------
        ShouldSkipLineError
            If `line` does not contain any valuable information that is
            worth returning while iterating over `self.stream`.
        """


# pylint: disable-next=too-few-public-methods   # Inherited from parent
class InputStreamReader(ABC, InputReader):
    """Class for reading an input file from a stream."""

    def __init__(self, source, noisy=True):
        """Initialize instance.

        Parameters
        ----------
        source : TextIOBase
            Input stream to be read.
        noisy : bool, optional
            Whether the reader will emit logging messages and raise
            errors if unknown or malformed lines are encountered.

        Raises
        ------
        TypeError
            If `source` is not a TextIOBase instance.
        """
        if not isinstance(source, TextIOBase):
            raise TypeError('Input source must be a TextIOBase type object.')
        self._source = source
        super().__init__(noisy=noisy)

    @property
    def stream(self):
        """Return the input stream."""
        return self._source


class InputFileReader(AbstractContextManager, InputReader):
    """A context manager that iterates the contents of an input file."""

    def __init__(self, filename, noisy=True):
        """Initialize instance.

        Parameters
        ----------
        filename : str or Path
            Path to the file to be read.
        noisy : bool, optional
            Whether the reader will emit logging messages and raise
            errors if unknown or malformed lines are encountered.

        Returns
        -------
        None.
        """
        self._filename = Path(filename)
        self._file_obj = None
        super().__init__(noisy=noisy)

    @property
    def stream(self):
        """Return the input file stream."""
        return self._file_obj

    def __enter__(self):
        """Enter context."""
        self._file_obj = self._filename.open('r', encoding='utf-8')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close file, then return control to caller to handle exceptions."""
        try:
            self._file_obj.close()
        except AttributeError:
            pass
        return super().__exit__(exc_type, exc_value, traceback)
