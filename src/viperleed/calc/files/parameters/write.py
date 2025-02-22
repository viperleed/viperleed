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

import numpy as np

from viperleed.calc.lib.string_utils import parent_name
from viperleed.calc.lib.string_utils import strip_comments
from viperleed.calc.lib.woods_notation import writeWoodsNotation

from .known_parameters import _ACCEPTS_MULTIPLE_ASSIGNMENTS
from .reader import RawLineParametersReader
from .utils import Assignment


_LOGGER = logging.getLogger(parent_name(__name__))


def comment_out(rpars, modpar, comment='', original=None):
    """Comment out modpar in the PARAMETERS file."""
    with ParametersFileEditor(rpars) as editor:
        editor.comment_out_parameter(modpar, comment, original)


def modify(rpars, modpar, new=None, comment='', original=None):
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
    original : Assignment, optional
        The user-given assignment for `modpar`, before modification.
        This argument is only used if `modpar` is a parameter for which
        all the user-given definitions are interpreted, rather than
        only the last one. Used to discern which of the assignments
        should be modified. It is a mandatory argument if `modpar`
        was given multiple times by the user.

    Returns
    -------
    new_value : ModifiedParameterValue
        Information about the new value of the parameter as inserted
        in the PARAMETERS file.
    """
    with ParametersFileEditor(rpars) as editor:
        new_param = editor.modify_param(modpar, new, comment, original)
    return new_param


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
        modified parameter when it will be written to file.
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
        'DOMAIN',
        'N_BULK_LAYERS',
        'SYMMETRY_FIX',  # For now
        }
    _vector_params = {   # (fmt_for_element, brackets)
        'THEO_ENERGIES': ('{:.4g}', ''),
        }
    _wood_or_matrix = {'SYMMETRY_CELL_TRANSFORM', 'SUPERLATTICE'}

    def __init__(self, param, rpars_or_value, original,
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
        original : Assignment or None
            The original value of the parameter that is to
            be modified. May be None if the parameter had
            no assignment (i.e., it has to be written anew).
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
        self.original = original
        self.comment = comment
        self._raw_value = None
        self.fmt_value = ''
        self.line = ''
        if not self.only_comment_out:
            self._update(rpars_or_value)

    def to_assignment(self):
        """Return an Assignment for this modified parameter.

        Returns
        -------
        assignment : Assignment or None
            An assignment for the edited value of this parameter.
            None if the parameter was only commented out.
        """
        if self.only_comment_out:
            return None
        flags = self.original.flags_str if self.original else ''
        return Assignment(self.fmt_value,
                          self.param,
                          raw_line=strip_comments(self.line),
                          flags_str=flags)

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
        flags = ('' if not self.original or not self.original.flags_str
                 else f' {self.original.flags_str}')
        line = (value if '=' in value
                else f'{self.param}{flags} = {value}')
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

    def __init__(self, rpars, path=''):
        """Initialize instance."""
        self._path = Path(path).resolve()
        self._rpars = rpars

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

    def comment_out_parameter(self, param, comment='', original=None):
        """Mark param as to be completely commented out."""
        assignments = self._infer_original_assignments(param, original)
        return self._comment_out_multiple_assignments(param,
                                                      assignments,
                                                      comment)

    def modify_param(self, param, new_value=None, comment='', original=None):   # Here we should probably comment out repeated single-assignment parameters except for the last one
        """Mark param as to be modified from the current value in rpars."""
        *to_comment, to_edit = self._infer_original_assignments(param, original)
        # Comment out any parameter that the user may have given
        # more than once (and for which we consider only the last
        # value).
        _duplicates = 'Duplicate values given. This was ignored.'
        self._comment_out_multiple_assignments(param, to_comment, _duplicates)

        # Finally, actually edit the one that needs to be edited
        if new_value is None:
            new_value = self._rpars
        modified = ModifiedParameterValue(param, new_value, to_edit, comment)
        self._to_modify[(param, to_edit)] = modified
        return modified

    def write_modified_parameters(self):
        """Write a new PARAMETERS, modifying all the requested lines."""
        if not self._to_modify:
            return

        # Collect the contents of the PARAMETERS
        # file before re-opening it for writing
        lines = self._collect_lines()
        self._has_header = False
        _head_mark = '! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #'
        with self._file.open('w', encoding='utf-8') as self._write_param_file:
            for param, assignment, raw_line in lines:
                if _head_mark in raw_line:
                    self._has_header = True
                self._write_one_line(param, assignment, raw_line)
            self._write_missing_parameters()
        self._to_modify.clear()

        # Finally, mark the file as being modified
        self._rpars.files_to_out.add(self._filename)

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

    def _comment_out_multiple_assignments(self, param, assignments, comment):
        """Mark multiple assignments for param as to be commented."""
        kwargs = {
            'param': param,
            'rpars_or_value': self._rpars,  # Not used anyway
            'comment': comment,
            'only_comment_out': True,
            }
        commented = []
        for original_assignment in assignments:
            comment_kwargs = kwargs.copy()
            comment_kwargs['original'] = original_assignment
            _edited = ModifiedParameterValue(**comment_kwargs)
            self._to_modify[(param, original_assignment)] = _edited
            commented.append(_edited)
        return commented

    def _get_comment_line_for(self, modified, raw_line):
        """Return a commented version of `raw_line`."""
        if modified.only_comment_out:
            comment = modified.comment or 'line commented out automatically'
            return f'! {raw_line.rstrip():<33} ! {comment}\n'
        return f'! {raw_line.rstrip():<33} ! line automatically changed to:\n'

    def _infer_original_assignments(self, param, original):
        """Find the relevant original assignments for param."""
        # Notice the use of .get instead of __getitem__: since
        # readParams is a defaultdict, using __getitem__ would
        # create a non-existing item. This may later screw up
        # with checks like "param in readParams".
        assignments = tuple(self._rpars.readParams.get(param, ()))
        if param not in _ACCEPTS_MULTIPLE_ASSIGNMENTS:
            return assignments or (None,)
        if len(assignments) == 1:
            return assignments
        if not original:
            raise TypeError(f'Cannot edit {param}: found multiple assignment '
                            'lines. Specify which line needs editing by '
                            'passing the \'original\' argument.')
        try:
            assignments.index(original)
        except ValueError:
            raise ValueError(f'Original {original} not found '
                             f'for {param} among the user-given '
                             f'ones: {assignments}') from None
        return (original,)

    def _is_unchanged(self, modified, raw_line):
        """Return whether `modified` is already present in `raw_line`."""
        return strip_comments(modified.line) == strip_comments(raw_line)

    def _update_read_params(self, modified):
        """Edit readParams to reflect the changes in `modified`."""
        read_ = self._rpars.readParams[modified.param]
        try:
            read_.remove(modified.original)
        except ValueError:
            pass
        if not modified.only_comment_out:
            read_.append(modified.to_assignment())
        if not read_:
            del self._rpars.readParams[modified.param]

    def _write_one_line(self, param, assignment, raw_line):
        """Write one `raw_line` containing param."""
        if not param:  # Empty, comment, or invalid line
            self._write_param_file.write(raw_line)
            return

        modified = self._to_modify.get((param, assignment), None)
        if modified is None or self._is_unchanged(modified, raw_line):
            # No modification needed
            self._write_param_file.write(raw_line)
            return

        commented_line = self._get_comment_line_for(modified, raw_line)
        self._write_param_file.write(commented_line)
        self._update_read_params(modified)
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
        for modified in to_be_written:
            self._update_read_params(modified)
        self._write_param_file.writelines(p.line for p in to_be_written)
