# -*- coding: utf-8 -*-
"""Module _modify of viperleed.tleedmlib.files.parameters.

Created on Tue Aug 18 16:56:39 2020

@author: Florian Kraushofer (@fkraushofer)
@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)

Initial version by @fkraushofer in 2020, major rewrite by @amimre
and @michele-riva in June 2023. This module used to be part of
parameters.py. Refactored in October 2023.

Functions for writing and editing a PARAMETERS file.
"""

import logging
from pathlib import Path
import shutil

from viperleed.tleedmlib.base import strip_comments


_LOGGER = logging.getLogger('tleedm.files.parameters')


def modifyPARAMETERS(rp, modpar, new='', comment='', path='',
                     suppress_ori=False, include_left=False):
    """
    Looks for 'modpar' in the PARAMETERS file, comments that line out, and
    replaces it by the string specified by 'new'

    Parameters
    ----------
    rp : Rparams
        The run parameters object.
    modpar : str
        The parameter that should be modified.
    new : str, optional
        The new value for the parameter. If not passed (default),
        the parameter will be commented out without replacement.
    comment : str, optional
        A comment to be added in the new line in PARAMETERS.
    path : str or Path, optional
        Where to find the PARAMETERS file that should be modified.
        Default is an empty string, i.e., the current directory.
    suppress_ori : bool, optional
        If True, no 'PARAMETERS_original' file will be created.
    include_left : bool, optional
        If False (default), 'new' will be interpreted as only the
        string that should go on the right-hand side of the equal
        sign. If True, the entire line will be replace by 'new'.

    Returns
    -------
    None.
    """
    _path = Path(path)
    file = _path / 'PARAMETERS'
    oriname = f'PARAMETERS_ori_{rp.timestamp}'
    ori = _path / oriname
    if oriname not in rp.manifest and not suppress_ori and file.is_file():
        try:
            shutil.copy2(file, ori)
        except Exception:
            _LOGGER.error(
                'modifyPARAMETERS: Could not copy PARAMETERS file to '
                'PARAMETERS_ori. Proceeding, original file will be lost.'
                )
        rp.manifest.append(oriname)
    if 'PARAMETERS' not in rp.manifest and _path == Path():
        rp.manifest.append('PARAMETERS')
    output = ''
    headerPrinted = False

    try:
        with file.open('r', encoding='utf-8') as parameters_file:
            lines = parameters_file.readlines()
    except FileNotFoundError:
        lines = []
    except Exception:
        _LOGGER.error('Error reading PARAMETERS file.')
        raise

    found = False
    for line in lines:
        if '! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #' in line:
            headerPrinted = True
        valid = False
        param = ''
        stripped_line = strip_comments(line)
        if '=' in stripped_line:
            # Parameter is defined left of '='
            param, *_ = stripped_line.split('=')
            if param:
                valid = True
                param, *_ = param.split()
        else:
            for p in ['STOP', 'SEARCH_KILL']:
                if stripped_line.upper().startswith(p):
                    valid = True
                    param = p
        if valid and param == modpar:
            found = True
            if new and f'{modpar} = {new.strip()}' == stripped_line:
                output += line
            elif new:
                output += f'!{line[:-1]} ! line automatically changed to:\n'
                if not include_left:
                    output += f'{modpar} = '
                output += f'{new} ! {comment}\n'
            else:
                comment = comment or 'line commented out automatically'
                output += f'!{line.rstrip():<34} ! {comment}\n'
        else:
            output += line
    if new and not found:
        if not headerPrinted:
            output += """

! ######################################################
! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #
! ######################################################
"""
        output += f'\n{modpar} = {new}'
        if comment:
            output += f' ! {comment}'
    try:
        with file.open('w', encoding='utf-8') as parameters_file:
            parameters_file.write(output)
    except Exception:
        _LOGGER.error('modifyPARAMETERS: Failed to write PARAMETERS file.')
        raise
