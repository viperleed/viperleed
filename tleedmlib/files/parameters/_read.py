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
    try:  # pylint: disable=too-many-try-statements
        with filename.open('r', encoding='utf-8') as param_file:
            lines = param_file.readlines()
    except FileNotFoundError:
        _LOGGER.error('updatePARAMETERS routine: PARAMETERS file not found.')
        raise

    # note no slab is given to the interpreter.
    # Slab is not needed for STOP, SEARCH_KILL and SEARCH_CONVERGENCE
    interpreter = ParameterInterpreter(rpars)
    for line in lines:
        line = strip_comments(line)
        for param in ['SEARCH_KILL', 'STOP']:  # SEARCH_KILL is legacy name
            if line.upper().startswith(param):
                if not re.match(fr'\s*{param}\s*=\s*[Ff](alse)?', line):
                    rpars.STOP = True
                    return  # We can stop right away

        # Ignore all lines that don't have an '=' sign at all
        if '=' not in line:
            continue

        param, value_str = line.split('=', maxsplit=1)
        if not param:  # Nothing left of '='
            continue

        param, *flags = param.split()
        try:
            param = from_alias(param)
        except ParameterNotRecognizedError:
            continue
        values = value_str.rstrip().split()
        if not values:
            # don't complain *again* in updatePARAMETERS
            continue
        if param == 'SEARCH_CONVERGENCE':
            new_assignment = Assignment(values_str=value_str,
                                        parameter=param,
                                        flags_str=flags)
            interpreter.interpret_search_convergence(assignment=new_assignment,
                                                     is_updating=True)
