# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 16:56:39 2020

@author: Florian Kraushofer
@author: Alexander M. Imre
@author: Michele Riva
Initial version by Florian Kraushofer in 2020, major rewrite by Alexander Imre
in June 2023.

Functions for reading from and writing to the PARAMETERS file
"""

from collections.abc import Sequence
from dataclasses import dataclass, field
from functools import partialmethod
import logging
from pathlib import Path
import re
import shutil

import numpy as np


from viperleed.tleedmlib.periodic_table import PERIODIC_TABLE
from viperleed.tleedmlib.base import (strip_comments, splitSublists,
                                      readVector, readIntRange,
                                      recombineListElements)
from viperleed.tleedmlib.classes import rparams
from viperleed.tleedmlib.files.parameter_errors import (
    ParameterError, ParameterValueError, ParameterParseError,
    ParameterIntConversionError, ParameterFloatConversionError,
    ParameterBooleanConversionError, ParameterNotRecognizedError,
    ParameterNumberOfInputsError, ParameterRangeError,
    ParameterUnknownFlagError, ParameterNeedsFlagError
    )
from viperleed.tleedmlib.files.woods_notation import readWoodsNotation
from viperleed.tleedmlib.sections._sections import TLEEDMSection as Section


logger = logging.getLogger('tleedm.files.parameters')


# list of allowed parameters
_KNOWN_PARAMS = (                                                               # TODO: IntEnum?
    'ATTENUATION_EPS', 'AVERAGE_BEAMS', 'BEAM_INCIDENCE', 'BULKDOUBLING_EPS',
    'BULKDOUBLING_MAX', 'BULK_LIKE_BELOW', 'BULK_REPEAT', 'DOMAIN',
    'DOMAIN_STEP', 'ELEMENT_MIX', 'ELEMENT_RENAME', 'FILAMENT_WF',
    'FORTRAN_COMP', 'HALTING', 'INTPOL_DEG', 'IV_SHIFT_RANGE', 'LAYER_CUTS',
    'LAYER_STACK_VERTICAL', 'LMAX', 'LOG_LEVEL', 'LOG_SEARCH', 'N_BULK_LAYERS',
    'N_CORES', 'OPTIMIZE', 'PARABOLA_FIT', 'PHASESHIFT_EPS',
    'PHASESHIFTS_CALC_OLD', 'PHASESHIFTS_OUT_OLD', 'PLOT_IV', 'RUN',
    'R_FACTOR_LEGACY', 'R_FACTOR_SMOOTH', 'R_FACTOR_TYPE', 'S_OVL',
    'SCREEN_APERTURE', 'SEARCH_BEAMS', 'SEARCH_CONVERGENCE', 'SEARCH_CULL',
    'SEARCH_MAX_GEN', 'SEARCH_POPULATION', 'SEARCH_START', 'SITE_DEF',
    'SUPERLATTICE', 'SUPPRESS_EXECUTION', 'SYMMETRIZE_INPUT', 'SYMMETRY_BULK',
    'SYMMETRY_CELL_TRANSFORM', 'SYMMETRY_EPS', 'SYMMETRY_FIND_ORI',
    'SYMMETRY_FIX', 'TENSOR_INDEX', 'TENSOR_OUTPUT', 'THEO_ENERGIES',
    'TL_VERSION', 'TL_IGNORE_CHECKSUM',
    'T_DEBYE', 'T_EXPERIMENT', 'V0_IMAG', 'V0_REAL',
    'V0_Z_ONSET', 'VIBR_AMP_SCALE', 'ZIP_COMPRESSION_LEVEL',
    )


# parameters that can be optimized in FD optimization
_OPTIMIZE_OPTIONS = {'theta', 'phi', 'v0i',
                     'a', 'b', 'c', 'ab', 'abc', 's_ovl'}


# _PARAM_ALIAS keys should be all lowercase, with no underscores
_PARAM_ALIAS = {
    'bulklike': 'BULK_LIKE_BELOW',
    'bulksymmetry': 'SYMMETRY_BULK',
    'compiler': 'FORTRAN_COMP',
    'fortrancompile': 'FORTRAN_COMP',
    'fortrancompiler': 'FORTRAN_COMP',
    'fdoptimize': 'OPTIMIZE',
    'fdoptimization': 'OPTIMIZE',
    'logdebug' : 'LOG_LEVEL',
    'plotrfactor': 'PLOT_IV',
    'plotrfactors': 'PLOT_IV',
    'ignorechecksum': 'TL_IGNORE_CHECKSUM',
    'ivplot': 'PLOT_IV',
    'overlap': 'S_OVL',
    'mtoverlap': 'S_OVL',
    'compression_level': 'ZIP_COMPRESSION_LEVEL',
    'compression': 'ZIP_COMPRESSION_LEVEL',
    }


for known_param in _KNOWN_PARAMS:
    _PARAM_ALIAS[known_param.lower().replace('_', '')] = known_param


# Bool parameters for which to create interpret...() methods automatically.
# key is parameter name, value is tuple of keyword arguments that should be
# passed to interpret_bool_parameter, in order.
_SIMPLE_BOOL_PARAMS = {
    'LOG_SEARCH' : (),
    'PHASESHIFTS_CALC_OLD' : (),
    'PHASESHIFTS_OUT_OLD' : (),
    'R_FACTOR_LEGACY' : (),
    'SUPPRESS_EXECUTION' : (),
    'SYMMETRIZE_INPUT' : (),
    'SYMMETRY_FIND_ORI' : (),
    'TL_IGNORE_CHECKSUM' : (),
    'LAYER_STACK_VERTICAL' : ({False: 'c', True: 'z'},),
    }


@dataclass(frozen=True)
class NumericBounds:
    """A container for bounds of numeric parameters.

    Attributes
    ----------
    type_ : {int, float}, optional
        Which numeric type these bounds are for.
    range : tuple, optional
        Minimum and maximum value for numeric parameter.
        Either limit can be None, signifying no limitation.
        Default is (None, None).
    accept_limits : tuple
        Whether the min/max endpoints of the range are
        acceptable values. Default is (True, True).
    out_of_range_event : {'fail', 'coerce', 'modulo'}, optional
        How to behave in case a number goes out of bounds.
        Default is 'fail'.
    """
    type_: type = float
    range_: tuple = (None, None)
    accept_limits: tuple = (True, True)
    out_of_range_event: str = 'fail'
    _msg_: str = field(init=False)

    def __post_init__(self):
        """Process assigned attributes.

        Raises
        ------
        ValueError
            If out_of_range_event is not one of the acceptable values
        ValueError
            If out_of_range_event is 'modulo' and some of the bounds
            were not specified
        ValueError
            If type_ is not one of the acceptable values
        """
        if self.out_of_range_event not in ('fail', 'modulo', 'coerce'):
            raise ValueError(f'Invalid {self.out_of_range_event=}. '
                             'Must be "fail", "coerce", or "modulo"')
        if self.type_ not in (int, float):
            raise ValueError('type_ must be int or float')
        is_modulo = self.out_of_range_event == 'modulo'
        if is_modulo and any(limit is None for limit in self.range_):
            raise ValueError('out_of_range_event "modulo" needs both limits')
        if is_modulo and not all(self.accept_limits):
            raise ValueError('out_of_range_event "modulo" needs both limits '
                             'to be acceptable.')

        # Use object.__setattr__ for frozen data-class
        object.__setattr__(self, '_msg_', self._make_message())

    @property
    def coerce(self):
        """Return whether the out-of-bounds reaction is 'coerce'."""
        return self.out_of_range_event == 'coerce'

    @property
    def fail(self):
        """Return whether the out-of-bounds reaction is 'fail'."""
        return self.out_of_range_event == 'fail'

    @property
    def max(self):
        """Return the lower bound."""
        return self.range_[1]

    @property
    def min(self):
        """Return the lower bound."""
        return self.range_[0]

    @property
    def closed_range(self):
        """A range within which values should fall in, including extremes."""
        increment = 1 if self.type_ is int else 1e-4  # float
        _min, _max = self.range_
        if _min is not None and not self.accept_limits[0]:
            _min += increment
        if _max is not None and not self.accept_limits[1]:
            _max -= increment
        return _min, _max

    @property
    def unlimited(self):
        """Return whether there is any bound."""
        return all(limit is None for limit in self.range_)

    def format_out_of_range(self, value):
        """Return a string for an out-of-range value."""
        return f'Value {value} is {self._msg_}'

    def is_in_range(self, val):
        """Return whether a value is within the acceptable bounds."""
        _min, _max = self.closed_range
        return (_min is None or val >= _min,
                _max is None or val <= _max,)

    def make_in_range(self, value):
        """Return a version of value that is within range."""
        in_range = self.is_in_range(value)
        if self.unlimited or all(in_range):
            return value
        if self.fail and not all(in_range):
            raise RuntimeError

        if self.coerce:
            _, value, _ = sorted((value, *self.closed_range))
            return value

        assert self.out_of_range_event == 'modulo'
        _min, _max = self.range_
        return (value - _min) % (_max - _min) + _min

    def _make_message(self):
        """Return a message for value-out-of-bounds."""
        if self.unlimited:
            return ''

        accept_min, accept_max = self.accept_limits
        range_chars = ({True: '[', False: ']'}, {True: ']', False: '['})
        if all(limit is not None for limit in self.range_):
            return (f'not in range {range_chars[0][accept_min]}{self.min}, '
                    f'{self.max}{range_chars[1][accept_max]}')
        if self.min is not None and accept_min:
            return f'less than {self.min}'
        if self.min is not None:
            return f'less than or equal to {self.min}'

        assert self.max is not None
        if accept_max:
            return f'larger than {self.max}'
        return f'larger than or equal to {self.max}'


_POSITIVE_INT = NumericBounds(type_=int, range_=(1, None))
_POSITIVE_FLOAT = NumericBounds(type_=float, range_=(0, None))


# Numerical parameters for which to create interpret...() methods
# automatically. Key is parameter name, value is a NumericBounds
_SIMPLE_NUMERICAL_PARAMS = {
    # Positive-only integers
    'BULKDOUBLING_MAX' : _POSITIVE_INT,
    'N_CORES' : _POSITIVE_INT,
    'SEARCH_MAX_GEN' : _POSITIVE_INT,
    'TENSOR_INDEX' : _POSITIVE_INT,
    # Positive-only floats
    'T_DEBYE' : _POSITIVE_FLOAT,
    'T_EXPERIMENT' : _POSITIVE_FLOAT,
    'V0_IMAG' : _POSITIVE_FLOAT,
    'TL_VERSION' : _POSITIVE_FLOAT,
    'S_OVL' : _POSITIVE_FLOAT,
    # Other floats
    'V0_Z_ONSET' : NumericBounds(),
    'ATTENUATION_EPS' : NumericBounds(range_=(1e-6, 1),
                                      accept_limits=(True, False)),
    'BULKDOUBLING_EPS' : NumericBounds(range_=(1e-4, None)),
    'BULK_LIKE_BELOW': NumericBounds(range_=(0, 1),
                                     accept_limits=(False, False)),
    'SCREEN_APERTURE' : NumericBounds(range_=(0, 180)),
    # Other integers
    'HALTING' : NumericBounds(type_=int, range_=(1, 3)),
    'N_BULK_LAYERS' : NumericBounds(type_=int, range_=(1, 2)),
    'R_FACTOR_SMOOTH' : NumericBounds(type_=int, range_=(0, 999)),
    'R_FACTOR_TYPE' : NumericBounds(type_=int, range_=(1, 2)),
    'ZIP_COMPRESSION_LEVEL' : NumericBounds(type_=int, range_=(0, 9))
    }


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

    interpreter = ParameterInterpreter(rp)
    interpreter.slab = rp.slab
    for line in lines:
        line = strip_comments(line)
        for param in ['SEARCH_KILL', 'STOP']:  # SEARCH_KILL is legacy name
            if line.upper().startswith(param):
                if not re.match(fr'\s*{param}\s*=\s*[Ff](alse)?', line):
                    rp.STOP = True
                    continue  # if need to STOP, we don't need continue interpreting the line
        if '=' not in line:                                                     # TODO: perhaps we should rather look at the first entry of a .split() and see if, by chance, it is a valid parameter. Then warn and continue reading.
            continue  # ignore all lines that don't have an '=' sign at all
        param, value_str = line.split('=', maxsplit=1)
        if param:
            # get rid of spaces and check the leftmost entry.
            param, *_ = param.split()
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
                                        parameter=param)
            interpreter.interpret_search_convergence(assignment=new_assignment,
                                                     is_updating=True)


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
        Object storing parameters for current run. Will contain parameters
        read in this function, and some 'global' parameters defined at runtime.
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
                modifyPARAMETERS(rpars, param, comment='Disabled at program '
                                 'start', path=filename.parent,
                                 suppress_ori=True)
        if '=' not in line:
            continue     # ignore all lines that don't have an '=' sign at all
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
        The new value for the parameter. If not passed (default), the parameter
        will be commented out without replacement.
    comment : str, optional
        A comment to be added in the new line in PARAMETERS.
    path : str or Path, optional
        Where to find the PARAMETERS file that should be modified.
        Default is an empty string, i.e., the current directory.
    suppress_ori : bool, optional
        If True, no 'PARAMETERS_original' file will be created.
    include_left : bool, optional
        If False (default), 'new' will be interpreted as only the string that
        should go on the right-hand side of the equal sign. If True, the entire
        line will be replace by 'new'.

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
            logger.error(
                'modifyPARAMETERS: Could not copy PARAMETERS file '
                'to PARAMETERS_ori. Proceeding, original file will be lost.')
        rp.manifest.append(oriname)
    if 'PARAMETERS' not in rp.manifest and _path == Path():
        rp.manifest.append('PARAMETERS')
    output = ''
    headerPrinted = False

    try:
        with open(file, 'r', encoding='utf-8') as rf:
            plines = rf.readlines()
    except FileNotFoundError:
        plines = []
    except Exception as err:
        logger.error('Error reading PARAMETERS file.')
        raise err
    found = False
    for line in plines:
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
        with open(file, 'w', encoding='utf-8') as wf:
            wf.write(output)
    except Exception:
        logger.error('modifyPARAMETERS: Failed to write PARAMETERS file.')
        raise


def interpretPARAMETERS(rpars, slab=None, silent=False):
    """Interpret rpars.readParams to actual values.

    Parameters
    ----------
    rpars : Rparams
        Object storing parameters for current run. Created previously by
        readPARAMETERS, and should already contain raw string data.
    slab : Slab, optional
        Slab object with elements and atomic position data. If not passed, some
        parameters will not be interpreted.
    silent : bool, optional
        If True, less output will be printed. The default is False.

    Raises
    ------
    ParameterError
        If any parameter interpretation fails.
    """
    interpreter = ParameterInterpreter(rpars)
    interpreter.interpret(slab=slab, silent=silent)


@dataclass(frozen=True)
class Assignment:
    """Class to store the flags and values of a parameter.

    Attributes
    ----------
    values_str : str
        The right-hand side of the parameter assignment. Can also be
        passed as a tuple of strings upon instantiation. In this case
        elements will be joined.
    flags_str : str
        The left-hand side of the parameter assignment, excluding the
        parameter itself. Can also be passed as a tuple of strings upon
        instantiation. In this case elements will be joined.
    parameter : str
        The PARAMETER this assignment corresponds to.
    flags : tuple
        The flags of the parameter as a tuple of strings.
        This is automatically generated from the string
        arguments given at instantiation.
    values : tuple
        The values of the parameter as a tuple of strings.
        This is automatically generated from the string
        arguments given at instantiation.
    """
    values_str: str
    parameter: str
    flags_str: str = ''
    flags: tuple = field(init=False)
    values: tuple = field(init=False)

    def __post_init__(self):
        """Split out left- and right-hand sides into flags and values."""
        flags = self._unpack_assignment_side(self.flags_str)                # TODO: are there any instances where we can have more than one flag per line?? If not we should raise always in that case!!
        values = self._unpack_assignment_side(self.values_str)

        object.__setattr__(self, 'flags', flags)
        object.__setattr__(self, 'values', values)

        # Make sure values_str and flags_str are actually strings:
        # we also accept Sequence of strings at __init__
        if not isinstance(self.values_str, str):
            object.__setattr__(self, 'values_str', ' '.join(self.values_str))
        if not isinstance(self.flags_str, str):
            object.__setattr__(self, 'flags_str', ' '.join(self.flags_str))

    @property
    def flag(self):
        """Return the leftmost flag as a string."""
        return self.flags[0] if self.flags else ''

    @property
    def value(self):
        """Return the leftmost value as a string."""
        return self.values[0] if self.values else ''

    @property
    def other_flags(self):
        """Return all the flags except for the first one."""
        return self.flags[1:]

    @property
    def other_values(self):
        """Return all the values except for the first one (as strings)."""
        return self.values[1:]

    def _unpack_assignment_side(self, side):
        """Split the side left or right of the equal into bits."""
        if isinstance(side, str):
            return tuple(side.split())
        if isinstance(side, Sequence):
            return tuple(side)
        raise TypeError


_BELOW_DEBUG = 2


class ParameterInterpreter:
    """Class to interpret parameters from the PARAMETERS file.

    To add a new parameter, add a method with the name 'interpret_<param>'.
    """

    domains_ignore_params = {
        'BULK_REPEAT', 'ELEMENT_MIX', 'ELEMENT_RENAME',
        'LAYER_CUTS', 'N_BULK_LAYERS', 'SITE_DEF', 'SUPERLATTICE',
        'SYMMETRY_CELL_TRANSFORM', 'TENSOR_INDEX', 'TENSOR_OUTPUT'
        }

    grouplist = [                                                               # TODO: take from elsewhere
        'p1', 'p2', 'pm', 'pg', 'cm', 'rcm', 'pmm', 'pmg', 'pgg', 'cmm',
        'rcmm', 'p4', 'p4m', 'p4g', 'p3', 'p3m1', 'p31m', 'p6', 'p6m']

    bool_synonyms = {
        True: {'true', 'yes', '1', 't', 'y', 'on', 'enable', '.true.'},
        False: {'false', 'no', '0', 'f', 'n', 'off', 'disable', '.false.'}
        }

    def __init__(self, rpars):
        """Initialize interpreter instance from an Rparams."""
        self.rpars = rpars
        self._slab = None
        self.param_names = []  # In precedence order

        # Some flags
        self._search_conv_read = False
        self.silent = False

    @property
    def slab(self):
        """Return the slab currently used for interpretation (or None)."""
        return self._slab

    @slab.setter
    def slab(self, slab):
        """Set the slab to be used for interpretation."""
        self._slab = slab

    def interpret(self, slab, silent=False):
        """Interpret all known parameters using slab."""
        self.slab = slab
        self._search_conv_read = False
        self.silent = silent

        _backup_log_level = logger.level
        if self.silent:
            logger.setLevel(logging.ERROR)

        self._update_param_order()

        for param, assignment in self._get_param_assignments():
            # Check if we are doing a domain calculation                        # TODO: this is a logical error! We may not know if we are doing a domain calc before yet. If a parameter is specified that should be ignored (e.g. BULK_REPEAT) prior to RUN this may crash!
            _is_domain_calc = 4 in self.rpars.RUN or self.rpars.domainParams
            if _is_domain_calc and param in self.domains_ignore_params:
                # skip in domain calculation
                continue

            self._interpret_param(param, assignment)
            logger.log(_BELOW_DEBUG,
                       f'Successfully interpreted parameter {param}')

        # finally set the log level back to what it was
        logger.setLevel(_backup_log_level)

    # ----------------  Helper methods for interpret() ----------------
    def _get_param_assignments(self):
        """Yield parameters and assignments for each PARAMETER read."""
        flat_params = (
            (p_name, *flags_and_value)
            for p_name in self.param_names
            for flags_and_value in self.rpars.readParams[p_name]
            )
        for param, flags, right_side in flat_params:
            assignment = Assignment(values_str=right_side,                      # TODO: use directly in readPARAMETERS
                                    flags_str=flags,
                                    parameter=param)
            yield param, assignment

    def _interpret_param(self, param, assignment):
        """Interpret the value of a single PARAMETER if known, or complain."""
        interpreter = getattr(self, f'interpret_{param.lower()}', None)
        if interpreter is None:
            # Complain about unknown parameter
            self.rpars.setHaltingLevel(2)
            raise ParameterNotRecognizedError(param)
        interpreter(assignment)

    def _update_param_order(self):
        """Define order in which parameters should be read."""
        ordered_params = 'LOG_LEVEL', 'RUN'                                     # TODO: I think DOMAIN should be in here too
        self.param_names = [p for p in ordered_params
                            if p in self.rpars.readParams]
        self.param_names.extend(
            p for p in _KNOWN_PARAMS
            if p in self.rpars.readParams and p not in self.param_names
            )

    # ---------- Methods for interpreting simple parameters -----------
    def interpret_bool_parameter(self, param, assignment,
                                 allowed_values=None,
                                 var_name=None, return_only=False):
        """Set a parameter to a boolean value.

        Parameters
        ----------
        param : str
            The name of the parameter in PARAMETERS.
        assignment : Assignment
            The assignment of the parameter, read from PARAMETERS.
        allowed_values : dict, optional
            Additional string values which should be interpreted as
            False or True. E.g. by default, 'f' and 'false' are False,
            't' and 'true' are True. Keys should be True and False,
            values are Sequences of acceptable strings.
        var_name : str, optional
            The variable name in the Rparams class, if it differs
            from 'param'.
        return_only: bool, optional
            If True, only return the value of the parameter, but do
            not set it.

        Returns
        -------
        value : bool
            The value of the parameter.
        """
        if var_name is None:
            var_name = param.upper()
        if allowed_values is None:
            allowed_values = {}

        _bool_synonyms = self.bool_synonyms.copy()
        for option, values in allowed_values.items():
            try:
                _bool_synonyms[option].update(v.lower() for v in values)
            except KeyError:
                raise ValueError(f'Unexpected {option=} '
                                 'in allowed_values') from None

        # Make sure there is no intersection between the two sets
        if (allowed_values
                and set(allowed_values[True]) & set(allowed_values[False])):
            raise ValueError('The sets of allowed values for '
                             'True and False must not overlap.')

        # Check if the value is in the allowed ones
        str_value = assignment.value.lower()
        try:
            value = next(bool_ for bool_, synonyms in _bool_synonyms.items()
                         if str_value in synonyms)
        except StopIteration:  # Value is invalid
            self.rpars.setHaltingLevel(1)
            raise ParameterBooleanConversionError(param,
                                                  assignment.value) from None

        if not return_only:
            setattr(self.rpars, var_name, value)
        return value

    def interpret_numerical_parameter(self, param, assignment,
                                      var_name=None, return_only=False,         # I don't understand the need for var_name. We already have param. Also, we could make param into a keyword, as we have .parameter in Assignment.
                                      bounds=NumericBounds()):                  # TODO: ideally one could default to actually using the limits known from self.rpars.get_limits!
        """Set a parameter to a numeric (int or float) value.

        Parameters
        ----------
        param : str
            The name of the parameter in PARAMETERS
        assignment : Assignment
            The assignment of the parameter, read from PARAMETERS.
        var_name : str, optional
            The variable name in the Rparams class, if it differs from
            'param'. Default is None.
        return_only: bool, optional
            If True, only return the value of the parameter, but do not
            set it.
        bounds : NumericBounds, optional
            Acceptable limits for value, and what to do in case an
            out-of-bounds value is given. Defaults is an unbounded
            float.

        Returns
        ----------
        value : int or float
            The interpreted value.

        Raises
        ------
        ParameterFloatConversionError
            When conversion of a float-type value to float fails.
        ParameterIntConversionError
            When conversion of a int-type value to int fails.
        ParameterRangeError
            When a the value falls outside the constraints of bounds.
        """
        type_ = bounds.type_
        try:
            # The inner float is necessary for, e.g., 1e6 as int
            value = type_(float(assignment.value))
        except ValueError:
            self.rpars.setHaltingLevel(1)
            if type_ is float:
                exc = ParameterFloatConversionError
            else:
                exc = ParameterIntConversionError
            raise exc(parameter=param) from None

        in_range = all(bounds.is_in_range(value))
        if not in_range and bounds.fail:
            self.rpars.setHaltingLevel(1)
            out_of_range = bounds.format_out_of_range(value)
            raise ParameterRangeError(param, message=out_of_range)

        if not in_range:
            in_range_value = bounds.make_in_range(value)
            out_of_range = bounds.format_out_of_range(value)
            logger.warning(f'PARAMETERS file: {param}: {out_of_range}. '
                           f'Value will be set to {in_range_value}.')
            value = in_range_value

        if var_name is None:
            var_name = param.upper()
        if not return_only:
            setattr(self.rpars, var_name, value)
        return value

    @classmethod
    def _make_methods(cls, wrapped, kwargs_names, new_methods_info):
        """Dynamically add methods for this class.

        Parameters
        ----------
        wrapped : callable
            The callable that will be wrapped to create methods.
            The a call to new_method(*args, **kwargs) becomes a
            wrapped(param, *args, **wrapped_kwargs, **kwargs)
            call, where param and wrapped_kwargs are generated
            here using kwargs_names and new_methods_info.
        kwargs_names : Sequence
            Names (str) of the keyword arguments that may be
            replaced. At most len(kwargs_names) keyword arguments
            will be given to wrapped.
        new_methods_info : dict
            Keys are names of the parameters for which methods are to
            be generated. They are used as the the first positional
            argument to wrapped. Also, the new methods will be
            named "interpret_<key>". Values are sequences of the
            values of keyword arguments to be passed on to wrapped,
            in the same order.

        Raises
        ------
        TypeError
            If too many keyword argument values are given in any of
            the new_methods_info
        """
        for param, kwargs_values in new_methods_info.items():
            method_name = f'interpret_{param.lower()}'
            if not isinstance(kwargs_values, Sequence):
                kwargs_values = (kwargs_values,)
            if len(kwargs_values) > len(kwargs_names):
                raise TypeError(
                    f'Too many keyword arguments for {method_name}. '
                    f'Expected at most values for {kwargs_names}'
                    )
            kwargs = dict(zip(kwargs_names, kwargs_values))
            method = partialmethod(wrapped, param, **kwargs)
            setattr(cls, method_name, method)

    @classmethod
    def make_boolean_interpreter_methods(cls):
        """Dynamically generate bool-setting methods."""
        cls._make_methods(cls.interpret_bool_parameter,
                          ('allowed_values',),
                          _SIMPLE_BOOL_PARAMS)

    @classmethod
    def make_numerical_interpreter_methods(cls):
        """Dynamically generate int/float-setting methods."""
        cls._make_methods(cls.interpret_numerical_parameter,
                          ('bounds',),
                          _SIMPLE_NUMERICAL_PARAMS)

    # ----- Methods to make sure no extra flags/values are given ------
    def _ensure_single_flag_assignment(self, param, assignment):
        """Raise if assignment contains multiple flags."""
        if assignment.other_flags:
            self.rpars.setHaltingLevel(1)
            raise ParameterUnknownFlagError(param, assignment.flags_str)

    def _ensure_single_value_assignment(self, param, assignment):
        """Raise if assignment contains multiple values."""
        if assignment.other_values:
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param,
                                               (len(assignment.values), 1))

    def _ensure_no_flags_assignment(self, param, assignment):
        """Raise if assignment contains any flag."""
        if assignment.flags:
            self.rpars.setHaltingLevel(1)
            raise ParameterUnknownFlagError(param, assignment.flag)  # highlight first flag

    def _ensure_simple_assignment(self, param, assignment):
        """Raise if assignment is not simple (i.e., one value, no flags)."""
        self._ensure_no_flags_assignment(param, assignment)
        self._ensure_single_value_assignment(param, assignment)

    def _ensure_single_flag_and_value_assignment(self, param, assignment):
        """Raise if assignment does not have a single flag--value pair."""
        self._ensure_single_flag_assignment(param, assignment)
        self._ensure_single_value_assignment(param, assignment)

    # ----------- Methods to interpret individual parameters ----------
    def interpret_average_beams(self, assignment):                              # TODO: not nice to have multiple types for this. May make sense to make it its own class with some helpers. Otherwise one has to do a bunch of isinstance checks when using it
        """Assign parameter AVERAGE_BEAMS."""
        param = 'AVERAGE_BEAMS'
        self._ensure_no_flags_assignment(param, assignment)
        right_side = assignment.values_str.lower().strip()
        # trivial cases
        if right_side in ['off', 'none', 'false', 'f']:
            self.rpars.AVERAGE_BEAMS = False
            return
        if right_side.lower() in ['all', 'perpendicular', 'perp']:
            self.rpars.AVERAGE_BEAMS = (0., 0.)
            return

        # Otherwise, try to parse the value
        angles = self._parse_incidence_angles(assignment, param, right_side)
        self.rpars.AVERAGE_BEAMS = angles['THETA'], angles['PHI']

    def interpret_beam_incidence(self, assignment):
        """Assign parameter BEAM_INCIDENCE."""
        param = 'BEAM_INCIDENCE'
        self._ensure_no_flags_assignment(param, assignment)
        right_side = assignment.values_str.lower().strip()
        angles = self._parse_incidence_angles(assignment, param, right_side)
        self.rpars.THETA, self.rpars.PHI = angles['THETA'], angles['PHI']

    def interpret_intpol_deg(self, assignment):
        """Assign parameter INTPOL_DEG."""
        param = 'INTPOL_DEG'
        self._ensure_simple_assignment(param, assignment)
        intpol_deg = assignment.value
        if intpol_deg in self.rpars.get_limits(param):
            self.rpars.INTPOL_DEG = int(intpol_deg)
            return
        message = 'Only degree 3 and 5 interpolation supported at the moment.'
        raise ParameterError(parameter=param, message=message)

    def interpret_bulk_repeat(self, assignment):
        """Assign parameter BULK_REPEAT."""
        param = 'BULK_REPEAT'
        # make sure that the slab is defined, otherwise bulk repeat is moot
        if not self.slab:
            raise ParameterError(parameter=param,
                                message="No slab defined for bulk repeat.")
        s = assignment.values_str.lower()
        # TODO: clean up and de-branch below
        if "[" not in s:
            if "(" not in s:
                try:
                    self.rpars.BULK_REPEAT = abs(float(assignment.value))
                except ValueError:
                    raise ParameterFloatConversionError(parameter=param)
            else:
                # regex to match e.g. c(2.0) or z(2.0)
                m = re.match(r'\s*(c|z)\(\s*(?P<val>[0-9.]+)\s*\)', s)
                if not m:
                    raise ParameterParseError(parameter=param)
                else:
                    try:
                        v = abs(float(m.group("val")))
                    except Exception:
                        raise ParameterFloatConversionError(parameter=param)
                    else:
                        if "z" in s:
                            self.rpars.BULK_REPEAT = v
                        else:  # c
                            self.rpars.BULK_REPEAT = self.slab.ucell[2, 2] * v
        else:  # vector
            vec = readVector(s, self.slab.ucell)
            if vec is None:
                raise ParameterParseError(parameter=param)
            self.rpars.BULK_REPEAT = vec

    def interpret_domain(self, assignment):
        """Set the domain path and name."""
        param = 'DOMAIN'
        # Check name
        name = assignment.flag
        names = [n for n, _ in self.rpars.DOMAINS]
        if name in names:  # Already defined
            error_message = f'Multiple sources defined for domain {name}.'
            self.rpars.setHaltingLevel(3)
            raise ParameterError(parameter=param, message=error_message)
        if not name:  # Get unique name                                         # TODO: used in several other places
            i = 1
            while str(i) in names:
                i += 1
            name = str(i)

        # Check path
        right_side = assignment.values_str.strip()
        if Path(right_side).exists():
            path = right_side
        elif Path(right_side).with_suffix('.zip').is_file():
            path = right_side + '.zip'
        else:
            error_message = (f'Value for DOMAIN {name} could not be '
                             'interpreted as either a path or a .zip file.')
            self.rpars.setHaltingLevel(3)
            raise ParameterError(parameter=param, message=error_message)
        self.rpars.DOMAINS.append((name, path))

    def interpret_domain_step(self, assignment):
        """Assign parameter DOMAIN_STEP."""
        param = 'DOMAIN_STEP'
        self._ensure_simple_assignment(param, assignment)
        try:
            i = int(assignment.value)
        except ValueError:
            self.rpars.setHaltingLevel(1)
            raise ParameterIntConversionError(parameter=param)
        if not (1 <= i <= 100):
            raise ParameterRangeError(parameter=param,
                                        given_value=i,
                                        allowed_range=(1, 100))
        domain_step = self.interpret_numerical_parameter(
            param, assignment,
            bounds=NumericBounds(type_=int, range_=(1, 100)),
            return_only=True
            )

        # pylint: disable=compare-to-zero
        # Seems clearer this way than having "if 100 % domain_step"
        if 100 % domain_step != 0:
            j = domain_step - 1
            while 100 % j != 0:
                j -= 1
            message = (f'100 is not divisible by given value {domain_step}. '
                       f'Consider using {j} instead.')
            self.rpars.setHaltingLevel(1)
            raise ParameterError(parameter=param,message=message)

        self.rpars.DOMAIN_STEP = domain_step

    def interpret_element_mix(self, assignment):                                # TODO: don't we check to avoid conflicts for ELEMENT_MIX and ELEMENT_RENAME? We should perhaps have a call to a checker after all parameters are read in?
        """Assign parameter ELEMENT_MIX."""
        param = 'ELEMENT_MIX'
        self._ensure_single_flag_assignment(param, assignment)
        ptl = [el.lower() for el in PERIODIC_TABLE]                             # TODO: nicer to use the leedbase element getter function and catch exceptions
        found = False
        for el in assignment.values:
            if el.lower() not in ptl:
                message = f"Element {el} not found in periodic table."
                self.rpars.setHaltingLevel(2)
                raise ParameterError(parameter=param, message=message)
        if not found:                                                           # TODO: this seems very odd; I'll leave it in for now because I'm unsure what the intended behaviour was
            self.rpars.ELEMENT_MIX[assignment.flag.capitalize()] = [
                el.capitalize() for el in assignment.values
                ]

    def interpret_element_rename(self, assignment):
        """Assign parameter ELEMENT_RENAME."""
        param = 'ELEMENT_RENAME'
        self._ensure_single_flag_and_value_assignment(param, assignment)
        ptl = [el.lower() for el in PERIODIC_TABLE]                             # TODO: nicer to use the leedbase element getter function and catch exceptions
        el = assignment.value.lower()
        if el not in ptl:
            message = f"Element {el} not found in periodic table."
            self.rpars.setHaltingLevel(2)
            raise ParameterError(parameter=param, message=message)
        else:
            self.rpars.ELEMENT_RENAME[assignment.flag.capitalize()] = (
                assignment.value.capitalize()
                )

    def interpret_filament_wf(self, assignment):
        """Assign parameter FILAMENT_WF."""
        param = 'FILAMENT_WF'
        self._ensure_simple_assignment(param, assignment)
        # check common filaments (e.g W), otherwise interpret as float
        known_filaments = self.rpars.get_default(param)
        try:
            self.rpars.FILAMENT_WF = known_filaments[assignment.value.lower()]
        except KeyError:
            self.interpret_numerical_parameter(param, assignment)

    def interpret_fortran_comp(self, assignment, skip_check=False):             # TODO: would be nicer to have a namedtuple or dataclass or similar. It could then have .pre, .post, .mpi, etc...
        """Assign parameter FORTRAN_COMP."""
        param = 'FORTRAN_COMP'

        flags, compiler_str = assignment.flags, assignment.value
        # complain if more than one flag is given
        if len(flags) > 1:
            message = (f"Only one flag allowed for {param} per line. "
                       f"Got {assignment.flags}.")
            raise ParameterError(param, message)
        if (not flags and compiler_str.lower() in ["ifort", "gfortran"]):
            # get default compiler flags
            self.rpars.getFortranComp(comp=compiler_str.lower(), skip_check=skip_check)
        elif (flags and assignment.flag.lower() == "mpi"
                    and compiler_str.lower() in ["mpifort", "mpiifort"]):
            # get default mpi compiler flags
            self.rpars.getFortranMpiComp(comp=compiler_str.lower(),
                                         skip_check=skip_check)
        else:  # set custom compiler flags
            delim = compiler_str[0]     # should be quotation marks                     # TODO: is this needed now that we use f-strings?
            if delim not in ["'", '"']:
                message = ("No valid shorthand and not delimited by "
                                "quotation marks.")
                raise ParameterError(param, message)
            else:
                setTo = assignment.values_str.split(delim)[1]
            if not flags:
                self.rpars.FORTRAN_COMP[0] = setTo
            elif flags[0].lower() == "post":
                self.rpars.FORTRAN_COMP[1] = setTo
            elif flags[0].lower() == "mpi":
                self.rpars.FORTRAN_COMP_MPI[0] = setTo
            elif flags[0].lower() == "mpipost":
                self.rpars.FORTRAN_COMP_MPI[1] = setTo
            else:
                raise ParameterUnknownFlagError(parameter=param,
                                                        flag=flags[0])

    def interpret_log_level(self, assignment):
        """Assign parameter LOG_LEVEL."""
        param = 'LOG_LEVEL'
        self._ensure_simple_assignment(param, assignment)

        # Try to interpret as bool. This is the way
        # the previous LOG_DEBUG used to work
        try:
            log_debug = self.interpret_bool_parameter(param, assignment,
                                                      return_only=True)
        except ParameterError:
            pass
        else:
            self.rpars.LOG_LEVEL = logging.DEBUG if log_debug else logging.INFO
            return

        # Otherwise interpret as int
        try:
            log_level = self.interpret_numerical_parameter(
                param, assignment,
                bounds=NumericBounds(type_=int, range_=(0, 50)),
                return_only=True
                )
        except ParameterError as exc:
            raise ParameterValueError(param) from exc
        self.rpars.LOG_LEVEL = log_level

    def interpret_v0_real(self, assignment):
        """Assign parameter V0_REAL."""
        param = 'V0_REAL'
        v0r_type = assignment.value.lower()
        if v0r_type == 'rundgren':
            rundgren_constants = assignment.other_values
            if len(rundgren_constants) != 4:
                message = ('Rundgren-type function expects four constants '
                           'separated by whitespace.')
                raise ParameterParseError(param, message)
            try:
                setTo = [float(rundgren_constants[i]) for i in range(4)]
            except (ValueError, IndexError) as err:
                message = ("Could not parse constants for Rundgren-type "
                            "function.")
                self.rpars.setHaltingLevel(1)
                raise ParameterError(param, message=message) from err
            return

        # Pass a specific function to FORTRAN, but replace
        # 'EE' with 'EEV+workfn' to get the energies right
        v0_real = re.sub('(?i)EE', 'EEV+workfn', assignment.values_str)
        self.rpars.V0_REAL = v0_real.rstrip()

    def interpret_symmetry_bulk(self, assignment):
        """Assign parameter SYMMETRY_BULK."""
        param = 'SYMMETRY_BULK'
        recombined_list = []                                                    # TODO: use ast
        values = assignment.all_values
        while values:
            v = values.pop(0).lower()
            if '[' in v and ']' in values[0]:
                v += ' ' + values.pop(0).lower()
            recombined_list.append(v)
        for v in recombined_list:
            if v.startswith('r'):                                               # Problem: how about "rcm[1 0]"?
                try:
                    i = int(v[1])
                except (ValueError, IndexError):
                    self.rpars.setHaltingLevel(2)
                    raise ParameterValueError(param, v)
                if 'rotation' not in self.rpars.SYMMETRY_BULK:
                    self.rpars.SYMMETRY_BULK['rotation'] = []
                if i not in self.rpars.SYMMETRY_BULK['rotation']:
                    self.rpars.SYMMETRY_BULK['rotation'].append(i)
            elif v.startswith('m'):
                if '[' not in v or ']' not in v:
                    message = f'Error reading value {v!r}: no direction recognized.'
                    self.rpars.setHaltingLevel(2)
                    raise ParameterParseError(param, message)
                str_vals = v.split('[')[1].split(']')[0].split()
                if len(str_vals) != 2:
                    self.rpars.setHaltingLevel(2)
                    raise ParameterNumberOfInputsError(
                        param,
                        found_and_expected=(len(str_vals), 2)
                        )
                try:
                    int_vals = tuple(int(v) for v in str_vals)
                except (ValueError, IndexError) as err:
                    self.rpars.setHaltingLevel(2)
                    raise ParameterValueError(param, v) from err
                if int_vals[0] < 0:
                    int_vals = (-int_vals[0], -int_vals[1])
                if 'mirror' not in self.rpars.SYMMETRY_BULK:
                    self.rpars.SYMMETRY_BULK['mirror'] = []
                if int_vals in self.rpars.SYMMETRY_BULK['mir']:
                    self.rpars.SYMMETRY_BULK['mirror'].append(int_vals)
            else:
                try:
                    self.rpars.SYMMETRY_BULK['group'] = v                       # This can never error out
                except ValueError as err:
                    self.rpars.setHaltingLevel(2)
                    raise ParameterValueError(param, v) from err

    def interpret_symmetry_fix(self, assignment):                               # TODO: use symmetry groups from elsewhere once symmetry and guilib are merged
        param = 'SYMMETRY_FIX'
        group = assignment.value.lower()
        if group.startswith('t'):
            return  # same as default, determine symmetry automatically
        if group.startswith('f'):
            self.rpars.SYMMETRY_FIX = 'p1'
            return
        if group in self.grouplist and group in ('cm', 'pmg'):
            message = f'For group {group} direction needs to be specified.'
            self.rpars.setHaltingLevel(1)
            raise ParameterParseError(param, message)
        if group in self.grouplist:
            self.rpars.SYMMETRY_FIX = group
            return
        if group.startswith(('pm', 'pg', 'cm', 'rcm', 'pmg')):
            # regex to read
            rgx = re.compile(
                r'\s*(?P<group>(pm|pg|cm|rcm|pmg))\s*'
                + r'\[\s*(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
            m = rgx.match(assignment.values_str.lower())
            if not m:
                self.rpars.setHaltingLevel(2)
                raise ParameterParseError(param)
            i1 = i2 = -2
            group = m.group('group')
            try:
                i1 = int(m.group('i1'))
                i2 = int(m.group('i2'))
            except ValueError as err:
                self.rpars.setHaltingLevel(2)
                raise ParameterParseError(param) from err
            if (group in ['pm', 'pg', 'cm', 'rcm', 'pmg']
                    and i1 in range(-1, 3) and i2 in range(-1, 3)):
                self.rpars.SYMMETRY_FIX = m.group(0)
            else:
                self.rpars.setHaltingLevel(2)
                raise ParameterParseError(param)
        else:
            self.rpars.setHaltingLevel(2)
            raise ParameterParseError(param)

    def interpret_theo_energies(self, assignment):
        """Assign parameter THEO_ENERGIES."""
        param = 'THEO_ENERGIES'
        energies = assignment.values
        if not assignment.other_values:
            energy = assignment.value
            # single value input - only one energy requested
            try:
                f = float(energy)
            except ValueError as err:
                self.rpars.setHaltingLevel(1)
                raise ParameterFloatConversionError(param, energy) from err
            else:
                if f > 0:
                    self.rpars.THEO_ENERGIES = [f, f, 1]
                    return
                else:
                    self.rpars.setHaltingLevel(2)
                    raise ParameterParseError(param, energy)
            return
        if len(energies) != 3:
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param, (len(energies), 3))

        fl = []
        defined = 0
        for val in energies:
            if val == '_':
                fl.append(-1)
                continue  # s in values loop
            try:
                f = float(val)
            except ValueError as err:
                self.rpars.setHaltingLevel(1)
                raise ParameterFloatConversionError(param, val) from err
            else:
                if f > 0:
                    fl.append(f)
                    defined += 1
                else:
                    message = f'{param} values have to be positive.'
                    self.rpars.setHaltingLevel(1)
                    raise ParameterError(param, message)

        if len(fl) != 3:
            raise ParameterNumberOfInputsError(param)
        if defined < 3:
            self.rpars.THEO_ENERGIES = fl
            return
        if (fl[0] > 0 and fl[1] > fl[0] and fl[2] > 0):
            if (fl[1] - fl[0]) % fl[2] != 0:
                # if the max is not hit by the steps exactly,
                #   correct max up to make it so
                fl[0] -= fl[2] - (fl[1] - fl[0]) % fl[2]
                if fl[0] <= 0:
                    fl[0] = fl[0] % fl[2]
                    if fl[0] == 0:
                        fl[0] += fl[2]
                logger.info('THEO_ENERGIES parameter: '
                            '(Eto-Efrom)%Estep != 0, Efrom was '
                            f'corrected to {fl[0]}')
            self.rpars.THEO_ENERGIES = fl
        else:
            self.rpars.setHaltingLevel(1)
            raise ParameterParseError(param)

    def interpret_search_beams(self, assignment):
        param = 'SEARCH_BEAMS'
        if assignment.other_values:
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param)
        value = assignment.value.lower()
        if value.startswith(('0', 'a')):
            self.rpars.SEARCH_BEAMS = 0
        elif value.startswith(('1', 'i')):
            self.rpars.SEARCH_BEAMS = 1
        elif value.startswith(('2', 'f')):
            self.rpars.SEARCH_BEAMS = 2
        else:
            raise ParameterValueError(param, value)

    def interpret_search_convergence(self, assignment, is_updating=False):
        param = 'SEARCH_CONVERGENCE'
        search_convergence_known = self.rpars.search_convergence_known

        if len(assignment.values) > 2:
            raise ParameterNumberOfInputsError(param)

        if (not assignment.flags and len(assignment.values) == 1 and
            assignment.value.lower() == 'off' and
            not is_updating):                                                   # TODO: this is the behaviour of updatePARAMETERS. Was skipping this intended there?
            self.rpars.GAUSSIAN_WIDTH_SCALING = 1.
            return

        if not assignment.flags:
            raise ParameterNeedsFlagError(param)
        elif assignment.flag.lower() not in ['dgen', 'gaussian']:
            raise ParameterUnknownFlagError(param, f'{assignment.flag!r}')

        try:
            numeric = [float(value) for value in assignment.values]
        except ValueError as err:
            raise ParameterFloatConversionError(param) from err

        _errors = []
        if assignment.flag.lower() == 'gaussian':
            if len(numeric) == 1:
                numeric.append(0.5)
            gauss_width, scaling = numeric
            should_update = (not is_updating
                             or gauss_width != self.rpars.searchConvInit['gaussian'])
            if gauss_width is not None and gauss_width > 0 and should_update:
                self.rpars.GAUSSIAN_WIDTH = gauss_width
                if is_updating:
                    self.rpars.searchConvInit['gaussian'] = gauss_width
            elif should_update:
                message = 'gaussian width should be a positive number '
                _errors.append(
                    ParameterError(param, message=message)
                    )
            if scaling is not None and 0 < scaling <= 1:
                self.rpars.GAUSSIAN_WIDTH_SCALING = scaling
            elif scaling is not None:
                _errors.append(
                    ParameterRangeError(param, scaling, (0, 1))
                    )
        elif assignment.flag.lower() == 'dgen':
            if len(assignment.flags) == 1:
                target = 'dec'
            elif assignment.flags[1].lower() in ['dec', 'best', 'all']:
                target = assignment.flags[1].lower()
            else:
                raise ParameterUnknownFlagError(param, f'{assignment.flags[1]!r}')
            max_dgen, scaling = numeric
            should_update = (not is_updating
                             or max_dgen != self.rpars.searchConvInit['dgen'][target])
            if max_dgen is not None and max_dgen > 0 and should_update:
                if not search_convergence_known:  # clear default values
                    self.rpars.SEARCH_MAX_DGEN = {'all': 0, 'best': 0, 'dec': 0}
                    search_convergence_known = True
                self.rpars.SEARCH_MAX_DGEN[target] = max_dgen
                if is_updating:
                    self.rpars.searchConvInit['dgen'][target] = max_dgen
            elif should_update:
                message = 'dgen should be a positive number '
                _errors.append(
                    ParameterError(param, message=message)
                    )
            if scaling is not None and scaling >= 1:
                self.rpars.SEARCH_MAX_DGEN_SCALING[target] = scaling
            elif scaling:
                message = 'scaling value cannot be smaller than 1.'
                _errors.append(
                    ParameterError(param, message=message)
                    )

        if _errors:
            raise _errors[0]
        self.rpars.search_convergence_known = search_convergence_known

    def interpret_search_cull(self, assignment):
        param = 'SEARCH_CULL'
        cull_value = assignment.value
        try:
            cull_float = float(cull_value)
        except ValueError as err:
            self.rpars.setHaltingLevel(1)
            raise ParameterFloatConversionError(param, cull_value) from err
        if cull_float >= 1:
            if cull_float - int(cull_float) < 1e-6:
                self.rpars.SEARCH_CULL = int(cull_float)
            else:
                message = f'{param} value has to be integer if greater than 1.'
                self.rpars.setHaltingLevel(1)
                raise ParameterError(param, message)
        elif cull_float >= 0:
            self.rpars.SEARCH_CULL = cull_float
        else:
            message = f'{param} value has to be non-negative.'
            self.rpars.setHaltingLevel(1)
            raise ParameterError(param, message)
        if assignment.other_values:
            next_value = assignment.other_values[0].lower()
            if next_value in ['clone', 'genetic', 'random']:
                self.rpars.SEARCH_CULL_TYPE = next_value
            else:
                self.rpars.setHaltingLevel(1)
                raise ParameterValueError(param, next_value)
            if len(assignment.other_values) > 1:
                self.rpars.setHaltingLevel(1)
                raise ParameterNumberOfInputsError(param)

    def interpret_search_population(self, assignment):
        param = 'SEARCH_POPULATION'
        self.interpret_numerical_parameter(
            param, assignment,
            bounds=NumericBounds(type_=int, range_=(1, None))
            )
        if self.rpars.SEARCH_POPULATION < 16:
            logger.warning('SEARCH_POPULATION is very small. A '
                        'minimum value of 16 is recommended.')

    def interpret_search_start(self, assignment):
        param = 'SEARCH_START'
        # there should only be one values
        if assignment.other_values:
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param)

        # there are only a few options, so we can just check them all
        search_start = assignment.value
        if search_start.startswith('rand'):
            self.rpars.SEARCH_START = 'random'
        elif search_start.startswith('center'):
            self.rpars.SEARCH_START = 'centered'
        elif search_start.startswith('control'):
            self.rpars.SEARCH_START = 'control'
        elif search_start.startswith('crandom'):
            self.rpars.SEARCH_START = 'crandom'
        else:
            self.rpars.setHaltingLevel(1)
            raise ParameterUnknownFlagError(param, search_start)

    def interpret_site_def(self, assignment):                                  # TODO: clean up this mess, also write tests
        param = 'SITE_DEF'
        newdict = {}
        sublists = splitSublists(assignment.values, ',')
        for sl in sublists:
            atnums = []
            for i in range(1, len(sl)):
                ir = readIntRange(sl[i])
                if len(ir) > 0:
                    atnums.extend(ir)
                elif 'top(' in sl[i]:
                    if self.slab is None:
                        self.rpars.setHaltingLevel(3)
                        raise ParameterError(
                            param,
                            ('SITE_DEF parameter contains a top() '
                             'function, but no slab was passed.')
                        )
                    n = int(sl[i].split('(')[1].split(')')[0])
                    csatlist = sorted(self.slab.atlist,
                                        key=lambda atom: atom.pos[2])
                    while n > 0:
                        at = csatlist.pop()
                        if at.el == assignment.flags[0]:
                            atnums.append(at.oriN)
                            n -= 1
                else:
                    self.rpars.setHaltingLevel(3)
                    raise ParameterError(param, 'Problem with SITE_DEF input format')
            newdict[sl[0]] = atnums
        self.rpars.SITE_DEF[assignment.flags[0]] = newdict

    def interpret_superlattice(self, assignment):
        param = 'SUPERLATTICE'
        if 'M' not in assignment.flags:
            self._ensure_no_flags_assignment(param, assignment)
            read_result = self._read_woods_notation(param, assignment)
            setattr(self.rpars, param, read_result)
            self.rpars.superlattice_defined = True
        else:
            matrix = self._read_matrix_notation(param, assignment.values)
            # write to param
            setattr(self.rpars, param, matrix)
            self.rpars.superlattice_defined = True

    def interpret_symmetry_cell_transform(self, assignment):
        param = 'SYMMETRY_CELL_TRANSFORM'
        if 'M' not in assignment.flags:
            self._ensure_no_flags_assignment(param, assignment)
            read_result = self._read_woods_notation(param, assignment)
            setattr(self.rpars, param, read_result)
        else:
            matrix = self._read_matrix_notation(param, assignment.values)
            # write to param
            setattr(self.rpars, param, matrix)

    def interpret_symmetry_eps(self, assignment):
        param = 'SYMMETRY_EPS'
        self._ensure_no_flags_assignment(param, assignment)
        if len(assignment.values) not in (1, 2):
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param)

        # warning specific to SYMMETRY_EPS
        warning_str = ('PARAMETERS file: SYMMETRY_EPS:\n'
                       'Given value {}is greater than one Ångström. This is a '
                       'very loose constraint and might lead to incorrect '
                       'symmetry detection. Be sure to check the output.')

        # interpret first value as SYMMETRY_EPS
        bounds = NumericBounds(type_=float, range_=(1e-100, None))
        self.interpret_numerical_parameter(param, assignment, bounds=bounds)
        if self.rpars.SYMMETRY_EPS > 1.0:
            logger.warning(warning_str.format(''))
        # interpret possible second value as SYMMETRY_EPS_Z
        if assignment.other_values:
            z_assignment = Assignment(assignment.other_values, param)
            self.interpret_numerical_parameter('SYMMETRY_EPS_Z',
                                               z_assignment,
                                               bounds=bounds)
            if self.rpars.SYMMETRY_EPS_Z > 1.0:
                logger.warning(warning_str.format('for z '))

    def interpret_tensor_output(self, assignment):
        param = 'TENSOR_OUTPUT'
        nl = recombineListElements(assignment.values, '*')
        for s in nl:
            s = re.sub(r'[Tt](rue|RUE)?', '1', s)
            s = re.sub(r'[Ff](alse|ALSE)?', '0', s)
            try:
                v = int(s)
                if v not in (0, 1):
                    message = ('Problem with TENSOR_OUTPUT input format: '
                               f'Found value {v}, expected 0 or 1.')
                    self.rpars.setHaltingLevel(1)
                    raise ParameterParseError(param, message)
                else:
                    self.rpars.TENSOR_OUTPUT.append(v)
            except ValueError:
                if '*' in s:
                    sl = s.split('*')
                    try:
                        r = int(sl[0])
                        v = int(sl[1])
                    except ValueError as err:
                        raise ParameterIntConversionError(param, sl) from err
                    if v not in (0, 1):
                        message = ('Problem with TENSOR_OUTPUT input format: '
                                   f'Found value {v}, expected 0 or 1.')
                        self.rpars.setHaltingLevel(1)
                        raise ParameterParseError(param, message)
                    else:
                        self.rpars.TENSOR_OUTPUT.extend([v]*r)
                else:
                    self.rpars.setHaltingLevel(1)
                    raise ParameterValueError(param, s)

    def interpret_run(self, assignment):                                       # TODO: important param, write tests
        param = 'RUN'
        run_list = []
        for section_str in assignment.values:
            try:
                run_list.extend(Section.sequence_from_string(section_str))
            except ValueError as err:
                self.rpars.setHaltingLevel(2)
                raise ParameterValueError(param, section_str) from err
        if Section.DOMAINS in run_list:
            logger.info('Found domain search.')
        if run_list:
            # insert initialization section if not present
            if run_list[0] is not Section.INITIALIZATION:
                run_list.insert(0, Section.INITIALIZATION)
            self.rpars.RUN = [s.value for s in run_list]                        # TODO: replace with "rl" to keep Section objects
        else:
            self.rpars.setHaltingLevel(3)
            raise ParameterError(
                param,
                'RUN was defined, but no values were read.'
                )

    def interpret_iv_shift_range(self, assignment):
        param = 'IV_SHIFT_RANGE'
        if len(assignment.values) not in (2, 3):
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(parameter=param)
        try:
            fl = [float(s) for s in assignment.values]
        except ValueError:
            raise ParameterFloatConversionError(parameter=param)
        if fl[1] < fl[0]:
            message = 'IV_SHIFT_RANGE end energy has to >= start energy.'
            self.rpars.setHaltingLevel(1)
            raise ParameterError(param, message)
            return

        for i in range(0, 2):
            self.rpars.IV_SHIFT_RANGE[i] = fl[i]
        if len(fl) == 3 and fl[2] <= 0:
            message = 'IV_SHIFT_RANGE step has to be positive.'
            self.rpars.setHaltingLevel(1)
            raise ParameterError(parameter=param, message=message)
        self.rpars.IV_SHIFT_RANGE[2] = fl[2]

    def interpret_layer_cuts(self, assignment):
        param = 'LAYER_CUTS'
        layer_cuts = list(assignment.values)
        # some simple filtering here, but leave as list of strings
        if all(c in assignment.values_str for c in '<>'):
            self.rpars.setHaltingLevel(1)
            raise ParameterParseError(param,
                                    'Cannot parse list with both "<" and ">".')
        if any(c in assignment.values_str for c in '<>'):
            layer_cuts = []
            for s in assignment.values:
                s = s.replace('<', ' < ')
                s = s.replace('>', ' > ')
                layer_cuts.extend(s.split())

        rgx = re.compile(r'\s*(dz|dc)\s*\(\s*(?P<cutoff>[0-9.]+)\s*\)')
        for (i, s) in enumerate(layer_cuts):
            if 'dz' in s.lower() or 'dc' in s.lower():
                m = rgx.match(assignment.values_str.lower())
                if m:
                    try:
                        float(m.group('cutoff'))
                        layer_cuts[i] = m.group(0)
                    except Exception:                                           # TODO: catch better; only 1 statement in try.
                        self.rpars.setHaltingLevel(1)
                        raise ParameterParseError(param,
                                                  f'Could not parse function {s}')
            elif s != '<' and s != '>':
                try:
                    float(s)
                except Exception:
                    self.rpars.setHaltingLevel(1)
                    raise ParameterParseError(param)
        self.rpars.LAYER_CUTS = layer_cuts

    def interpret_lmax(self, assignment):
        param = 'LMAX'
        _min, _max = self.rpars.get_limits(param)

        values = re.sub(r'[:-]', ' ', assignment.values_str).split()
        try:
            lmax_list = [int(v) for v in values]
        except ValueError:
            self.rpars.setHaltingLevel(1)
            raise ParameterIntConversionError(param)
        if len(lmax_list) > 2:
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param)
        if len(lmax_list) == 1:
            if not _min < lmax_list[0] <= _max:
                _, lmax_list[0], _ = sorted((_min, lmax_list[0], _max))
                raise ParameterError(
                    param,
                    f'LMAX must be between {_min} and {_max}.'
                )
            self.rpars.LMAX = [lmax_list[0], lmax_list[0]]
        elif len(lmax_list) == 2:
            if lmax_list[1] < lmax_list[0]:
                lmax_list.reverse()
            if lmax_list[0] < _min:
                raise ParameterError(
                    param,
                    'LMAX lower bound must be positive.'
                )
            if lmax_list[1] > _max:
                raise ParameterError(
                    param,
                    f'LMAX values >{_max} are currently not supported.'
                    )
            self.rpars.LMAX = lmax_list

    def interpret_optimize(self, assignment):
        param = 'OPTIMIZE'
        if not assignment.flag:
            message = 'Parameter to optimize not defined.'
            self.rpars.setHaltingLevel(3)
            raise ParameterError(param, message)
        which = assignment.flag.lower()
        if which not in _OPTIMIZE_OPTIONS:
            self.rpars.setHaltingLevel(3)
            raise ParameterUnknownFlagError(param, f'{which!r}')
        self.rpars.OPTIMIZE['which'] = which
        if not assignment.other_values:
            try:
                self.rpars.OPTIMIZE['step'] = float(assignment.value)
            except ValueError:
                pass   # will be caught below
            else:
                return
        sublists = splitSublists(assignment.values, ',')
        for sl in sublists:
            if len(sl) != 2:
                message = 'Expected "flag value" pairs, found ' + ' '.join(sl)
                self.rpars.setHaltingLevel(2)
                raise ParameterError(param, message)
            flag = sl[0].lower()
            if flag not in ['step', 'convergence',
                            'minpoints', 'maxpoints', 'maxstep']:
                self.rpars.setHaltingLevel(2)
                raise ParameterUnknownFlagError(param, f'{flag!r}')
            partype = {'step': float, 'convergence': float,
                        'minpoints': int, 'maxpoints': int,
                        'maxstep': float}
            value_error = ('PARAMETERS file: OPTIMIZE: Value '
                            f'{sl[1]} is not valid for flag {sl[0]}. '
                            'Value will be ignored.')
            try:
                self.rpars.OPTIMIZE[flag] = partype[flag](sl[1])
            except ValueError as err:
                self.rpars.setHaltingLevel(1)
                raise ParameterError(param, value_error) from err

    def interpret_parabola_fit(self, assignment):
        param = 'PARABOLA_FIT'
        if assignment.value.lower() == 'off':
            self.rpars.PARABOLA_FIT['type'] = 'none'
            return
        sublists = splitSublists(assignment.values, ',')
        for sl in sublists:
            flag = sl[0].lower()
            if flag.lower() == 'localise':
                flag = 'localize'
            value_error = (f'PARAMETERS file: PARABOLA_FIT: Value {sl[1]} '
                            f'is not valid for flag {sl[0]}. Value will be '
                            'ignored.')
            if flag == 'type':
                if sl[1].lower() in ('linear', 'linearregression', 'lasso',
                                        'ridge', 'elasticnet', 'none'):
                    self.rpars.PARABOLA_FIT['type'] = sl[1]
                else:
                    self.rpars.setHaltingLevel(1)
                    raise ParameterError(param, value_error)
            elif flag in ('alpha', 'mincurv', 'localize'):
                try:
                    f = float(sl[1])
                except ValueError:
                    f = -1
                if f >= 0:
                    self.rpars.PARABOLA_FIT[flag] = f
                else:
                    self.rpars.setHaltingLevel(1)
                    raise ParameterError(param, value_error)

    def interpret_phaseshift_eps(self, assignment):
        param = 'PHASESHIFT_EPS'
        ps_eps_value = assignment.value
        try:
            ps_eps = float(ps_eps_value)
        except ValueError:
            # check if one of default values (e.g. 'fine')
            s = ps_eps_value.lower()[0]
            ps_eps_default_dict = self.rpars.get_default(param)
            ps_eps = ps_eps_default_dict.get(s, None)
            if ps_eps is None:
                self.rpars.setHaltingLevel(1)
                raise ParameterFloatConversionError(param)
        if 0 < ps_eps < 1:
            self.rpars.PHASESHIFT_EPS = ps_eps
        else:
            self.rpars.setHaltingLevel(1)
            raise ParameterRangeError(param,
                                      given_value=ps_eps, allowed_range=(0,1))

    def interpret_plot_iv(self, assignment):
        param = 'PLOT_IV'
        # there should be a flag
        if not assignment.flags:
            self.rpars.setHaltingLevel(1)
            raise ParameterNeedsFlagError(param)

        flag = assignment.flags[0].lower()
        value = assignment.value.lower()
        if flag not in ('color', 'colour', 'colors', 'colours', 'perpage',
                        'border', 'borders', 'axes', 'legend', 'legends',
                        'layout', 'overbar', 'overline', 'plot'):
            self.rpars.setHaltingLevel(1)
            raise ParameterUnknownFlagError(param, f'{flag!r}')
        if flag == 'plot':
            # should it plot?
            if value in ('true', 't', '1'):
                self.rpars.PLOT_IV['plot'] = True
            elif value in ('false', 'none', 'f', '0'):
                self.rpars.PLOT_IV['plot'] = False
        if flag in ('border', 'borders', 'axes'):
            if value in ('all', 'none'):
                self.rpars.PLOT_IV['axes'] = value
            elif value in ('less', 'lb'):
                self.rpars.PLOT_IV['axes'] = 'lb'
            elif value in ('bottom', 'b'):
                self.rpars.PLOT_IV['axes'] = 'b'
            else:
                self.rpars.setHaltingLevel(1)
                raise ParameterParseError(param)
        elif flag in ('color', 'colour', 'colors', 'colours'):
            self.rpars.PLOT_IV['colors'] = assignment.values
        elif flag in ('legend', 'legends'):
            if value in ('all', 'first', 'none'):
                self.rpars.PLOT_IV['legend'] = value
            elif value in ('topright', 'tr'):
                self.rpars.PLOT_IV['legend'] = 'tr'
            else:
                self.rpars.setHaltingLevel(1)
                raise ParameterParseError(param)
        elif flag in ('overbar', 'overline'):
            if value.startswith('t'):
                self.rpars.PLOT_IV['overbar'] = True
            elif value.startswith('f'):
                self.rpars.PLOT_IV['overbar'] = False
            else:
                message = f'Value for flag {flag} not recognized'
                self.rpars.setHaltingLevel(1)
                raise ParameterParseError(param, message)
        elif flag in ('perpage', 'layout'):
            if not assignment.other_values:
                try:
                    i = int(value)
                except (ValueError, IndexError):
                    self.rpars.setHaltingLevel(1)
                    raise ParameterIntConversionError(param, value)
                if i <= 0:
                    message = 'perpage value has to be positive integer.'
                    self.rpars.setHaltingLevel(1)
                    raise ParameterParseError(param, message)
                self.rpars.PLOT_IV['perpage'] = i
            elif len(assignment.values) >= 2:
                try:
                    il = [int(v) for v in assignment.values[:2]]
                except ValueError:
                    self.rpars.setHaltingLevel(1)
                    raise ParameterIntConversionError(param,
                                                      assignment.all_[:2])
                if any(i <= 0 for i in il):
                    message = 'perpage values have to be positive integers.'
                    raise ParameterParseError(param, message)
                self.rpars.PLOT_IV['perpage'] = tuple(il)

    def interpret_vibr_amp_scale(self, assignment):
        # this parameter is not interpreted right now
        self.rpars.VIBR_AMP_SCALE.extend(assignment.values_str.split(','))

    # Helper methods
    def _parse_incidence_angles(self, assignment, param, right_side):
        bounds = {
            'THETA': NumericBounds(type_=float, range_=(-90, 90),
                                   out_of_range_event='fail'),
            'PHI': NumericBounds(type_=float, range_=(0, 360),
                                 out_of_range_event='modulo')
            }
        d = {'THETA': 0, 'PHI': 0}
        if ',' in right_side:
            sublists = splitSublists(assignment.values, ',')
            for sl in sublists:
                for name in ['THETA', 'PHI']:
                    if sl[0].upper() == name:
                        d[name] = self.interpret_numerical_parameter(
                            f'{param} {name}',
                            Assignment(sl[1], param),
                            bounds=bounds[name],
                            return_only=True
                            )
                        break
                else:
                    self.rpars.setHaltingLevel(1)
                    raise ParameterUnknownFlagError(parameter=param,
                                                    flag=sl[0])
        else:
            if len(assignment.values) != 2:
                self.rpars.setHaltingLevel(1)
                raise ParameterNumberOfInputsError(
                    parameter=param,
                    found_and_expected=(len(assignment.values), 2)
                    )
            for value, name in zip(assignment.values, ('THETA', 'PHI')):
                d[name] = self.interpret_numerical_parameter(
                    f'{param} {name}',
                    Assignment(value, param),
                    bounds=bounds[name],
                    return_only=True
                    )

        if any(v is None for v in d.values()):
            self.rpars.setHaltingLevel(1)
            raise ParameterParseError(parameter=param)
        # check for negative theta and adjust phi
        if d['THETA'] < 0:
            d['THETA'] = abs(d['THETA'])
            d['PHI'] = (d['PHI'] + 180) % 360
        return d

    def _read_woods_notation(self, param, assignment):
        if self.slab is None:
            message = (f'{param} parameter appears to be in Wood, '
                       'notation but no slab was passed. Cannot '
                       'calculate bulk unit cell!')
            self.rpars.setHaltingLevel(2)
            raise ParameterError(param, message)
        read_result = readWoodsNotation(assignment.values_str, self.slab.ucell)
        return read_result

    def _read_matrix_notation(self, param, values):
        sublists = splitSublists(values, ',')
        if len(sublists) != 2:
            message = 'Number of lines in matrix is not equal to 2.'
            self.rpars.setHaltingLevel(2)
            raise ParameterParseError(param, message)
        nl = []
        for sl in sublists:
            if len(sl) == 2:
                try:
                    nl.append([float(s) for s in sl])
                except ValueError:
                    self.rpars.setHaltingLevel(1)
                    raise ParameterFloatConversionError(param, sl)
            else:
                self.rpars.setHaltingLevel(2)
                message = 'Number of columns in matrix is not equal to 2.'
                raise ParameterParseError(param, message)
        matrix = np.array(nl, dtype=float)
        return matrix


# Dynamically produce methods for the 'simple parameters' listed above
# Notice that we should do this only once right here, not for each
# instance separately. Calling the methods again is possible, and
# would simply override the bindings.
ParameterInterpreter.make_boolean_interpreter_methods()
ParameterInterpreter.make_numerical_interpreter_methods()
