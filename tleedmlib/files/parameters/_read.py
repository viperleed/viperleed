# -*- coding: utf-8 -*-
"""Module _read of viperleed.tleedmlib.files.parameters.

Created on Tue Aug 18 16:56:39 2020

@author: Florian Kraushofer (@fkraushofer)
@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)

Initial version by @fkraushofer in 2020, major rewrite by @amimre
and @michele-riva in June 2023. This module used to be part of
parameters.py. Refactored in October 2023.

Functions for reading from a PARAMETERS file and for updating
an Rparams object at runtime from a user-modified PARAMETERS file.
"""

from contextlib import AbstractContextManager
import logging
from pathlib import Path
import re

from viperleed.tleedmlib.base import strip_comments
from viperleed.tleedmlib.classes import rparams

from .errors import ParameterNotRecognizedError
from ._known_parameters import from_alias
from ._interpret import ParameterInterpreter
from ._utils import Assignment
from ._write import modifyPARAMETERS


_LOGGER = logging.getLogger('tleedm.files.parameters')


def readPARAMETERS(filename='PARAMETERS'):
    """Return an Rparams with the raw contents read from a PARAMETERS file.

    Parameters
    ----------
    filename : str or Path, optional
        The file to be read. The default is 'PARAMETERS'.

    Returns
    -------
    rpars : Rparams
        Object storing parameters for current run. Contains the
        raw parameters read in this function in its `.readParams`
        attribute. The parameters read are not interpreted. For
        that, call `interpretPARAMETEERS` passing the the same
        `rpars` object.

    Raises
    ------
    FileNotFoundError
        If filename does not exist.
    ParameterNotRecognizedError
        If one of the parameters read from filename is not
        a known one, or if a parameter read has no value.
    """
    filename = Path(filename).resolve()
    if not filename.is_file():
        _LOGGER.error('PARAMETERS file not found.')
        raise FileNotFoundError(filename)

    with filename.open('r', encoding='utf-8') as param_file:
        lines = param_file.readlines()

    # read PARAMETERS:
    rpars = rparams.Rparams()
    for line in lines:
        line = strip_comments(line)
        for param in ['STOP', 'SEARCH_KILL']:
            if (line.upper().startswith(param)
                    and not re.match(fr'\s*{param}\s*=\s*[Ff](alse)?', line)):
                _LOGGER.warning(
                    f'PARAMETERS file: {param} was set at start of '
                    f'program. Modifying PARAMETERS to disable {param}; '
                    're-insert it if you actually want to stop.'
                    )
                modifyPARAMETERS(rpars, param,
                                 comment='Disabled at program start',
                                 path=filename.parent,
                                 suppress_ori=True)
        if '=' not in line:
            continue  # ignore all lines that don't have an '=' sign            # TODO: we should probably still check whether the line starts with something that looks like a parameter and warn. Can easily happen to forget an "=".
        param, value = line.split('=', maxsplit=1)  # parameter at left of '='
        if not param:
            continue
        # get rid of spaces and check the leftmost entry.
        param, *flags = param.split()
        param = from_alias(param)
        value = value.strip()
        if not value:
            rpars.setHaltingLevel(1)
            raise ParameterNotRecognizedError(parameter=param)
        rpars.readParams[param].append((flags, value))
    return rpars


def updatePARAMETERS(rpars, filename='PARAMETERS', update_from=''):
    """Update `rpars` from file, ignoring non-SEARCH-related parameters.

    The following parameters are considered:
    SEARCH_CONVERGENCE, SEARCH_KILL, STOP

    Parameters
    ----------
    rpars : Rparams
        Parameters for current run. Its SEARCH-related members
        are updated if the contents of `filename` has changed.
    filename : str or Path, optional
        The file to be read. The default is 'PARAMETERS'.
    update_from : str or Path, optional
        Path to the directory in which to look for `filename`.
        The default is an empty string, corresponding to the
        current directory.

    Raises
    ------
    FileNotFoundError
        If no `filename` is found in `update_from`.
    ParameterError
        If one of the parameters to be updated has an invalid value.
    """
    filename = Path(update_from, filename).resolve()
    if not filename.is_file():
        _LOGGER.error('updatePARAMETERS routine: PARAMETERS file not found.')
        raise FileNotFoundError(filename)

    # Note that no slab is given to the interpreter. It is not needed
    # for STOP, SEARCH_KILL and SEARCH_CONVERGENCE. Also, note that
    # we don't complain (again) about faulty parameters while reading
    interpreter = ParameterInterpreter(rpars)
    with ParametersReader(filename, noisy=False) as param_file:
        for param, assignment in param_file:
            if param == 'STOP':
                rpars.STOP = assignment
            if rpars.STOP:
                return  # No need to continue reading
            if param == 'SEARCH_CONVERGENCE':
                interpreter.interpret_search_convergence(assignment=assignment,
                                                         is_updating=True)


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
            line = strip_comments(line)

            # We treat STOP differently, as its presence, even without
            # an equal sign, is interpreted as 'please, stop right now'
            if self._stop_requested(line):
                yield 'STOP', True

            if '=' not in line:
                self._warn_about_missing_equals(line, line_nr)
                continue

            try:
                param, assignment = self._parse_line(line)
            except ParameterNotRecognizedError:
                if not self.noisy:
                    continue
                raise
            if not param:
                continue
            yield param, assignment

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
