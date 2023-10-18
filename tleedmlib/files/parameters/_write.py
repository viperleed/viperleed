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

from collections.abc import Sequence
from contextlib import AbstractContextManager
import logging
from pathlib import Path
import shutil

import numpy as np

from viperleed.tleedmlib.base import strip_comments
from viperleed.tleedmlib.files.woods_notation import writeWoodsNotation


_LOGGER = logging.getLogger('tleedm.files.parameters')

# Some weird string that will unlikely be the new value for
# a parameter. Used as default for modify to signal that
# we want to comment out a line
_COMMENT_OUT = '__comment_out_the_parameter_requested_alkjv__'


def comment_out(rpars, modpar, comment='', path='', suppress_ori=False):
    """Comment out modpar in the PARAMETERS file."""
    modify(rpars, modpar, new=_COMMENT_OUT, comment=comment,
           path=path, suppress_ori=suppress_ori, include_left=False)


def modify(rp, modpar, new=_COMMENT_OUT, comment='', path='',
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
        The new value for the parameter. If not passed, the
        parameter will be commented out without replacement.
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
                'parameters.modify: Could not copy PARAMETERS file to '
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
            if (new != _COMMENT_OUT
                    and f'{modpar} = {new.strip()}' == stripped_line):
                output += line
            elif new != _COMMENT_OUT:
                output += f'!{line[:-1]} ! line automatically changed to:\n'
                if not include_left:
                    output += f'{modpar} = '
                output += f'{new} ! {comment}\n'
            else:
                comment = comment or 'line commented out automatically'
                output += f'!{line.rstrip():<34} ! {comment}\n'
        else:
            output += line
    if new != _COMMENT_OUT and not found:
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
        _LOGGER.error('parameters.modify: Failed to write PARAMETERS file.')
        raise


# This is almost a dataclass (but not quite), all the
# information is static and generated only once via
# private methods. We don't really need public ones.
# pylint: disable-next=too-few-public-methods
class ModifiedParameterValue:
    """A container for a single modified PARAMETER.

    Attributes
    ----------
    already_written : bool
        Whether this parameter was already written to file.
    comment : str
        An optional comment that should accompany this
        modified parameter when it will be written to file
    fmt_value : str
        A formatted version of the new value. This normally
        includes only the part on the right side of the equals
        sign. `fmt_value` is the empty string if this parameter
        should be `only_comment_out`.
    line : str
        A formatted version of the full line to be written
        to file when modifying the parameter. This contains
        also the comment, if given, and is `\n`-terminated.
        `line` is the empty string if the parameter should
        be `only_comment_out`.
    only_comment_out : bool
        Whether this parameter should only be commented out.
    param : str
        The PARAMETERS parameter to be modified.
    """

    _simple_params = {   # Those with a specific format string
        'V0_IMAG': '.4f',
        'BULK_LIKE_BELOW': '.4f',
        }
    _string_params = {   # Those for which str is enough
        'N_BULK_LAYERS',
        'SYMMETRY_FIX',  # For now
        }
    _vector_params = {   # (fmt_for_element, brackets)
        'THEO_ENERGIES': ('{:.4g}', ''),
        }
    _wood_or_matrix = {'SYMMETRY_CELL_TRANSFORM', 'SUPERLATTICE'}

    def __init__(self, param, rpars_or_value,
                 comment='', only_comment_out=False):
        """Initialize instance.

        Parameters
        ----------
        param : str
            The PARAMETERS parameter whose value is modified.
        rpars_or_value : object
            If an Rparams object, the new value is taken from
            the current value of the corresponding attribute,
            otherwise, it is treated as the raw version of the
            new value.
        comment : str, optional
            Extra comments to add to the line that should be
            written.
        only_comment_out : bool, optional
            Whether the parameter should only be removed by
            commenting, and not actually modified.

        Returns
        -------
        None.
        """
        self.only_comment_out = only_comment_out
        self.already_written = False
        self.param = param
        self.comment = comment
        self._raw_value = None
        self.fmt_value = ''
        self.line = ''
        if not self.only_comment_out:
            self._update_(rpars_or_value)

    def _update_(self, rpars_or_value):
        """Gather the value, its formatted version and a line."""
        try:
            rpars_or_value.BULK_REPEAT
        except AttributeError:  # Assume it's a value
            self._raw_value = rpars_or_value
        else:                   # Assume it's an Rparams
            self._raw_value = self._get_raw_value(rpars_or_value)
        self.fmt_value = self._format_value()
        self.line = self._format_line()

    def _format_value(self):
        """Return a formatted version of raw_value for param."""
        param = self.param
        if param in self._string_params:
            return str(self._raw_value)

        fmt_string = self._simple_params.get(param, None)
        if fmt_string:
            return f'{self._raw_value:{fmt_string}}'

        if param in self._wood_or_matrix:
            return self._format_wood_or_matrix_value_()

        if param in self._vector_params:
            args = self._vector_params[param]
            return self._format_vector(self._raw_value, *args)

        try:
            formatter = getattr(self, f'_format_{param.lower()}_value')
        except AttributeError:
            raise ValueError(f'No formatting available for {param}') from None
        return formatter()

    def _format_line(self):
        """Return a formatted version of the new line for a PARAMETER."""
        if '=' in self.fmt_value:
            fmt_line = self.fmt_value
        else:
            fmt_line = f'{self.param} = {self.fmt_value}'
        if self.comment:
            fmt_line = f'{fmt_line.rstrip():<35} ! {self.comment}'
        return fmt_line + '\n'

    def _format_beam_incidence_value(self):
        """Return a formatted version of BEAM_INCIDENCE."""
        theta, phi = self._raw_value
        return f'THETA {theta:.4f}, PHI {phi:.4f}'

    def _format_bulk_repeat_value(self):
        """Return a formatted version of BULK_REPEAT."""
        bulk_repeat = self._raw_value
        if isinstance(bulk_repeat, (np.ndarray, Sequence)):
            return self._format_vector(bulk_repeat, '{:.5f}', brackets='[]')
        return f'{bulk_repeat:.5f}'  # Z distance

    def _format_layer_cuts_value(self):
        """Return a formatted version of LAYER_CUTS."""
        layer_cuts = self._raw_value
        return ' '.join(f'{cut:.4f}' for cut in layer_cuts)                     # TODO: BETTER!

    def _format_lmax_value(self):
        """Return a formatted version of LMAX."""
        min_, max_ = self._raw_value
        if min_ == max_:
            return str(min_)
        return f'{min_}-{max_}'

    def _format_wood_or_matrix_value_(self):
        """Return a formatted version of an integer matrix value."""
        formatted = writeWoodsNotation(self._raw_value)
        if formatted:
            formatted = f'= {formatted}'
        else:
            matrix = self._raw_value.round().astype(int)
            # Disable pylint, as the *ravel way looks cleaner
            # than an f-string with explicit matrix elements
            # pylint: disable-next=consider-using-f-string
            formatted = 'M = {} {}, {} {}'.format(*matrix.ravel())
        return f'{self.param} {formatted}'

    @staticmethod
    def _format_vector(vector, elem_fmt, brackets=''):
        """Return a formatted version of a 1D vector."""
        elements = ' '.join(elem_fmt.format(v) for v in vector)
        if not brackets:
            return elements
        left, right = brackets
        return left + elements + right

    def _get_raw_value(self, rpars):
        """Return the current value of a parameter from rpars."""
        try:
            return getattr(rpars, self.param)
        except AttributeError:
            pass
        getter = getattr(self, f'_get_raw_value_{self.param.lower()}', None)
        if getter is None:
            raise ValueError(f'No raw-value getter for {self.param}')
        return getter(rpars)

    @staticmethod
    def _get_raw_value_beam_incidence(rpars):
        """Return theta and phi."""
        return rpars.THETA, rpars.PHI
