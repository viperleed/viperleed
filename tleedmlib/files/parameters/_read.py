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

import ast
from collections.abc import Sequence
from dataclasses import dataclass, field
from functools import partialmethod
import logging
from pathlib import Path
import re
import shutil

import numpy as np

from viperleed.tleedmlib.base import (strip_comments, splitSublists,
                                      readVector, readIntRange,
                                      recombineListElements)
from viperleed.tleedmlib.classes import rparams
from viperleed.tleedmlib.files.woods_notation import readWoodsNotation
from viperleed.tleedmlib import periodic_table
from viperleed.tleedmlib.sections._sections import TLEEDMSection as Section

from .errors import (
    ParameterError, ParameterValueError, ParameterParseError,
    ParameterIntConversionError, ParameterFloatConversionError,
    ParameterBooleanConversionError, ParameterNotRecognizedError,
    ParameterNumberOfInputsError, ParameterRangeError,
    ParameterUnknownFlagError, ParameterNeedsFlagError
    )
from ._known_parameters import _KNOWN_PARAMS, _PARAM_ALIAS
from ._interpret import ParameterInterpreter
from ._utils import Assignment
from ._write import modifyPARAMETERS


logger = logging.getLogger('tleedm.files.parameters')


def readPARAMETERS(filename='PARAMETERS'):
    """
    Reads a PARAMETERS file and returns an Rparams object with the
    raw input, without interpretation.

    Parameters
    ----------
    filename : str, optional
        The file to be read. The default is 'PARAMETERS'.

    Returns
    -------
    rpars : Rparams
        Object storing parameters for current run. Will contain
        parameters read in this function, and some 'global'
        parameters defined at runtime.
    """
    filename = Path(filename).resolve()
    try:
        rf = filename.open('r', encoding='utf-8')
    except FileNotFoundError as err:
        logger.error('PARAMETERS file not found.')
        raise err
    # read PARAMETERS:
    rpars = rparams.Rparams()
    for line in rf:
        line = strip_comments(line)
        for param in ['STOP', 'SEARCH_KILL']:
            if (line.upper().startswith(param)
                    and not re.match(fr'\s*{param}\s*=\s*[Ff](alse)?', line)):
                logger.warning(
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
        if (param not in _KNOWN_PARAMS and
                param.lower().replace('_', '') in _PARAM_ALIAS):
            param = _PARAM_ALIAS[param.lower().replace('_', '')]
        if param not in _KNOWN_PARAMS:
            raise ParameterNotRecognizedError(parameter=param)
        value = value.strip()
        if not value:
            raise ParameterNotRecognizedError(parameter=param)
        if param not in rpars.readParams:
            rpars.readParams[param] = []
        rpars.readParams[param].append((flags, value))
    rf.close()
    return rpars


def updatePARAMETERS(rp, filename='PARAMETERS', update_from=''):
    """
    Reads PARAMETERS file again, but ignores everything not concerning the
    search or STOP. Updates the given Rparams object accordingly.

    Parameters
    ----------
    rp : Rparams
        Parameters for current run, as defined previously. Will be updated if
        parameters have changed.
    filename : str, optional
        The file to be read. The default is 'PARAMETERS'.
    update_from : str or path-like, optional                                    # TODO: missing doc

    Returns
    -------
    None.
    """
    update_from = Path(update_from)
    try:
        with open(update_from / filename, 'r', encoding='utf-8') as rf:
            lines = rf.readlines()
    except FileNotFoundError:
        logger.error('updatePARAMETERS routine: PARAMETERS file not found.')
        raise

    # note no slab is given to the interpreter.
    # Slab is not needed for STOP, SEARCH_KILL and SEARCH_CONVERGENCE
    interpreter = ParameterInterpreter(rp)
    for line in lines:
        line = strip_comments(line)
        for param in ['SEARCH_KILL', 'STOP']:  # SEARCH_KILL is legacy name
            if line.upper().startswith(param):
                if not re.match(fr'\s*{param}\s*=\s*[Ff](alse)?', line):
                    rp.STOP = True
                    continue  # if need to STOP, we don't need continue interpreting the line
        if '=' not in line:
            continue  # ignore all lines that don't have an '=' sign at all
        param, value_str = line.split('=', maxsplit=1)
        if param:
            # get rid of spaces and check the leftmost entry.
            param, *flags = param.split()
        else:
            flags = []
        param_alias = param.lower().replace('_', '')
        if param not in _KNOWN_PARAMS and param_alias in _PARAM_ALIAS:
            param = _PARAM_ALIAS[param_alias]
        if param not in _KNOWN_PARAMS:
            continue
        values = value_str.rstrip().split()
        if not values:
            # don't complain *again* in updatePARAMETERS
            continue
        if param == 'SEARCH_CONVERGENCE':
            new_assignment = Assignment(values_str=value_str,
                                        parameter=param,
                                        flags_str=' '.join(flags))
            interpreter.interpret_search_convergence(assignment=new_assignment,
                                                     is_updating=True)
