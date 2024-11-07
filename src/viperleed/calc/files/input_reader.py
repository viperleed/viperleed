"""Module input_reader of viperleed.calc.files.

This module is a more abstract version of an input-file reader. It is based on
@michele-riva's work on the PARAMETERS file reader, which in turn is based on
the old readPARAMETERS from @fkraushofer.

Defines context-manager, iterator classes for reading parameters from a file
as well as from TextIOBase streams.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-03'
__license__ = 'GPLv3+'

from abc import abstractmethod
from io import TextIOBase
from collections.abc import Iterator
from contextlib import AbstractContextManager
import logging
from pathlib import Path

_LOGGER = logging.getLogger(__name__)

# TODO: Add subclasses for the InputReaders as suggested by @michele-riva in
#       PR #266. E.g. DelimitedFileReader, BlockFileReader, etc.

class ShouldSkipLineError(Exception):
    """Exception raised when a line's contents should be skipped."""


class InputReader(Iterator):
    """Common base class for all input readers."""

    def __init__(self, noisy=True) -> None:
        self.noisy = noisy
        super().__init__()

    def __next__(self):
        """Return the next understandable information in the file."""
        for line in self.stream:
            try:
                return self._read_one_line(line)
            except ShouldSkipLineError:
                # Intentional sub-debug logging level as this is expected
                # behavior and may spam the logs.
                _LOGGER.log(level=5,
                            msg=f'Skipping line in input file: "{line}".')
                continue
        raise StopIteration

    @property
    @abstractmethod
    def stream(self):
        """Return the input stream."""
        pass

    @abstractmethod
    def _read_one_line(self, line):
        """Return understandable information from `line`.

        This method is guaranteed to be called once on each
        line read from the file.

        Parameters
        -------------
        line : str
            A single line read from `self._file_obj.`

        Returns
        --------
        info : object
            Understandable information read from `line`. This is the
            same object returned whenever this reader is iterated over.

        Raises
        -------
        ShouldSkipLineError
            If `line` does not contain any valuable information that is
            worth returning while iterating over the file.
        """
        pass


class InputStreamReader(InputReader):
    """Class for reading an input file from a stream."""
    def __init__(self, source, noisy=True):
        """Initialize instance.

        Parameters
        ----------
        source : Stream
            Input stream to be read.
        noisy : bool, optional
            Whether the reader will emit logging messages and raise
            errors if unknown or malformed parameters are encountered.
        """
        if not isinstance(source, TextIOBase):
            raise TypeError('Input source must be a stream-like object.')
        self._source = source
        self.noisy = noisy
        self._current_line = 0
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
            errors if unknown or malformed parameters are encountered.
        """
        self._filename = Path(filename)
        if not self._filename.is_file():
            raise FileNotFoundError(f'File {self._filename} does not exist.')
        self._file_obj = None
        self._current_line = 0
        super().__init__(noisy=noisy)

    @property
    def stream(self):
        """Return the input file stream."""
        return self._file_obj

    def __enter__(self):
        """Enter context."""
        self._file_obj = self._filename.open('r', encoding='utf-8')
        self._current_line = 0
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close file, then return control to caller to handle exceptions."""
        if self._file_obj is None:
            return super().__exit__(exc_type, exc_value, traceback)
        try:
            self._file_obj.close()
        except AttributeError:
            pass
        return super().__exit__(exc_type, exc_value, traceback)
