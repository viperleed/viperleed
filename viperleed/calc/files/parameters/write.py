"""Module write of viperleed.calc.files.parameters.

Initial version by @fkraushofer in 2020, major rewrite by @amimre
and @michele-riva in June 2023. This module used to be part of
parameters.py. Refactored in October 2023.

Functions for writing and editing a PARAMETERS file.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-18'
__license__ = 'GPLv3+'

from collections.abc import Sequence
from contextlib import AbstractContextManager
import logging
from pathlib import Path
import shutil

import numpy as np

from viperleed.calc.lib.base import strip_comments
from viperleed.calc.lib.woods_notation import writeWoodsNotation

from .reader import RawLineParametersReader


_LOGGER_NAME, _ = __name__.rsplit('.', maxsplit=1)
_LOGGER = logging.getLogger(_LOGGER_NAME)


def comment_out(rpars, modpar, comment='', path='', suppress_ori=False):
    """Comment out modpar in the PARAMETERS file."""
    editor = ParametersFileEditor(
        rpars, path=path,
        save_existing_parameters_file=not suppress_ori
        )
    with editor:
        editor.comment_out_parameter(modpar, comment=comment)


def modify(rpars, modpar, new=None, comment='', path='', suppress_ori=False):
    """Change the value for `modpar` in the PARAMETERS file.

    The lines that contains `modpar` are commented out, and
    replaced by the string specified by `new`, or the value
    currently stored in `rpars`.

    Parameters
    ----------
    rpars : Rparams
        The run parameters object.
    modpar : str
        The parameter that should be modified.
    new : object or None, optional
        The new value for the parameter. If None, the current
        value of `modpar` in `rpars` is used instead. Default
        is None.
    comment : str, optional
        A comment to be added in the new line in PARAMETERS.
    path : str or Path, optional
        Where to find the PARAMETERS file that should be modified.
        Default is an empty string, i.e., the current directory.
    suppress_ori : bool, optional
        If True, no 'PARAMETERS_original' file will be created.

    Returns
    -------
    new_value : str
        String value of the parameter as inserted in the
        PARAMETERS file.
    """
    editor = ParametersFileEditor(
        rpars, path=path,
        save_existing_parameters_file=not suppress_ori
        )
    with editor:
        new_param = editor.modify_param(modpar, new_value=new, comment=comment)
    return new_param.fmt_value


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
        'BULK_LIKE_BELOW': '.4f',
        'LAYER_CUTS': '.4f',
        'V0_IMAG': '.4f',
        }
    _string_params = {   # Those for which str is enough
        'LMAX',
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
            self._update(rpars_or_value)

    def _update(self, rpars_or_value):
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
            return self._format_wood_or_matrix_value()

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
        value = self.fmt_value
        line = value if '=' in value else f'{self.param} = {value}'
        if self.comment:
            line = f'{line.rstrip():<35} ! {self.comment}'
        return line + '\n'

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

    def _format_wood_or_matrix_value(self):
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


class ParametersFileEditor(AbstractContextManager):
    """A context manager for editing a PARAMETERS file."""

    _filename = 'PARAMETERS'
    _header = """

! ######################################################
! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #
! ######################################################

"""

    def __init__(self, rpars, path='', save_existing_parameters_file=True):
        """Initialize instance."""
        self._path = Path(path).resolve()
        self._rpars = rpars
        self._should_save_existing_file = save_existing_parameters_file

        self._to_modify = {}  # Parameters to modify/comment out
        self._has_header = False
        self._write_param_file = None  # The file object to write to

    @property
    def _file(self):
        """Return the path to the PARAMETERS file to modify."""
        return self._path / self._filename

    def __enter__(self):
        """Prepare to modify a PARAMETERS file."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Modify PARAMETERS file before exiting."""
        self.write_modified_parameters()
        return super().__exit__(exc_type, exc_value, traceback)

    def comment_out_parameter(self, param, comment=''):
        """Mark param as to be completely commented out."""
        commented = ModifiedParameterValue(param, self._rpars,
                                           comment=comment,
                                           only_comment_out=True)
        self._to_modify[param] = commented
        return commented

    def modify_param(self, param, new_value=None, comment=''):
        """Mark param as to be modified from the current value in rpars."""
        if new_value is None:
            new_value = self._rpars
        mod_param = ModifiedParameterValue(param, new_value, comment=comment)
        self._to_modify[param] = mod_param
        return mod_param

    def save_existing_parameters_file(self):
        """Save a copy of an existing PARAMETERS file as 'ori_<timestamp>'."""
        if not self._should_save_existing_file:
            return

        oriname = f'{self._filename}_ori_{self._rpars.timestamp}'
        if not self._file.is_file() or oriname in self._rpars.manifest:
            self._should_save_existing_file = False
            return

        try:
            shutil.copy2(self._file, self._path / oriname)
        except OSError:
            _LOGGER.error(
                f'parameters.modify: Could not copy {self._filename} file to '
                f'{self._filename}_ori. Proceeding. Original file will be lost'
                )
        self._rpars.manifest.append(oriname)
        self._should_save_existing_file = False

    def write_modified_parameters(self):
        """Write a new PARAMETERS, modifying all the requested lines."""
        self.save_existing_parameters_file()
        cwd = Path().resolve()
        if self._filename not in self._rpars.manifest and self._path == cwd:
            self._rpars.manifest.append(self._filename)

        if not self._to_modify:
            return

        # Collect the contents of the PARAMETERS
        # file before re-opening it for writing
        lines = self._collect_lines()
        self._has_header = False
        _head_mark = '! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #'
        with self._file.open('w', encoding='utf-8') as self._write_param_file:
            for param, raw_line in lines:
                if _head_mark in raw_line:
                    self._has_header = True
                self._write_one_line(param, raw_line)
            self._write_missing_parameters()
        self._to_modify.clear()

    def _collect_lines(self):
        """Return lines read from the PARAMETERS file."""
        if not self._file.is_file():
            return ()

        reader = RawLineParametersReader(self._file, noisy=False)
        try:  # pylint: disable=too-many-try-statements
            with reader:
                return tuple(reader)
        except Exception:
            _LOGGER.error(f'Error reading {self._filename} file.')
            raise

    def _get_comment_line_for(self, modified, raw_line):
        """Return a commented version of `raw_line`."""
        if modified.only_comment_out:
            comment = modified.comment or 'line commented out automatically'
            return f'! {raw_line.rstrip():<33} ! {comment}\n'
        return f'! {raw_line.rstrip():<33} ! line automatically changed to:\n'

    def _is_unchanged(self, modified, raw_line):
        """Return whether `modified` is already present in `raw_line`."""
        return strip_comments(modified.line) == strip_comments(raw_line)

    def _write_one_line(self, param, raw_line):
        """Write one `raw_line` containing param."""
        if not param:  # Empty, comment, or invalid line
            self._write_param_file.write(raw_line)
            return

        modified = self._to_modify.get(param, None)
        if modified is None or self._is_unchanged(modified, raw_line):
            # No modification needed
            self._write_param_file.write(raw_line)
            return

        commented_line = self._get_comment_line_for(modified, raw_line)
        self._write_param_file.write(commented_line)
        self._write_new_parameter_line(modified)

    def _write_new_parameter_line(self, modified):
        """Write one line for `modified`, if necessary."""
        if modified.only_comment_out or modified.already_written:
            return
        self._write_param_file.write(modified.line)
        modified.already_written = True

    def _write_missing_parameters(self):
        """Write all the parameters that were not found."""
        to_be_written = [p for p in self._to_modify.values()
                         if not p.already_written]
        if not to_be_written:
            return
        if not self._has_header:
            self._write_param_file.write(self._header)
        self._write_param_file.writelines(p.line for p in to_be_written)
