"""Module read of viperleed.calc.files.parameters.

Initial version by @fkraushofer in 2020, major rewrite by @amimre
and @michele-riva in June 2023. This module used to be part of
parameters.py. Refactored in October 2023.

Functions for reading from a PARAMETERS file and for updating
an Rparams object at runtime from a user-modified PARAMETERS file.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-18'
__license__ = 'GPLv3+'

import logging
from pathlib import Path

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.lib.context import execute_in_dir
from viperleed.calc.lib.string_utils import parent_name

from .errors import MissingEqualsError
from .interpret import ParameterInterpreter
from .reader import ParametersReader
from .write import comment_out


_LOGGER = logging.getLogger(parent_name(__name__))


def read(filename='PARAMETERS'):
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
        that, call `parameters.interpret` passing the the same
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

    rpars = Rparams()
    comment_out_stop = False
    with ParametersReader(filename, noisy=True) as param_file:
        while True:
            try:
                param, assignment = next(param_file)
            except StopIteration:
                break
            except MissingEqualsError as exc:
                _LOGGER.warning(exc)
                rpars.setHaltingLevel(2)
                continue
            if param == 'STOP':
                comment_out_stop = True
            rpars.readParams[param].append(assignment)

    if comment_out_stop:
        _LOGGER.warning(
            'PARAMETERS file: STOP was set at start of '
            'program. Modifying PARAMETERS to disable STOP; '
            're-insert it if you actually want to stop.'
            )
        with execute_in_dir(filename.parent):
            comment_out(rpars, 'STOP', comment='Disabled at program start')
    return rpars


def update(rpars, filename='PARAMETERS', update_from=''):
    """Update `rpars` from file.

    The following parameters are considered:
    SEARCH_CONVERGENCE, STOP (and its legacy names).

    Parameters
    ----------
    rpars : Rparams
        Parameters for current run. Its attributes are updated
        if the contents of `filename` has changed.
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
        _LOGGER.error('parameters.update: PARAMETERS file not found.')
        raise FileNotFoundError(filename)

    # Note that no slab is given to the interpreter. It is not needed
    # for STOP (+legacy names) and SEARCH_CONVERGENCE. Also, note that
    # we don't complain (again) about faulty parameters while reading
    interpreter = ParameterInterpreter(rpars)
    with ParametersReader(filename, noisy=False) as param_file:
        for param, assignment in param_file:
            if param == 'STOP':
                # pylint: disable=no-member
                # Method is added dynamically
                interpreter.interpret_stop(assignment)
            if rpars.STOP:
                return  # No need to continue reading
            if param == 'SEARCH_CONVERGENCE':
                interpreter.interpret_search_convergence(assignment=assignment,
                                                         is_updating=True)
