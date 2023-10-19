# -*- coding: utf-8 -*-
"""Module _interpret of viperleed.tleedmlib.files.parameters.

Created on Tue Aug 18 16:56:39 2020

@author: Florian Kraushofer (@fkraushofer)
@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)

Initial version by @fkraushofer in 2020, major rewrite by @amimre
and @michele-riva in June 2023. This module used to be part of
parameters.py. Refactored in October 2023.

Functions and classes for interpreting the contents previously
read from a PARAMETERS file.
"""

import ast
from collections.abc import Sequence
from functools import partialmethod
import logging
from pathlib import Path
import re

import numpy as np

from viperleed.tleedmlib import periodic_table
from viperleed.tleedmlib.base import readIntRange, readVector
from viperleed.tleedmlib.base import recombineListElements, splitSublists
from viperleed.tleedmlib.files.woods_notation import readWoodsNotation
from viperleed.tleedmlib.sections._sections import TLEEDMSection as Section

from .errors import ParameterError, ParameterBooleanConversionError
from .errors import ParameterFloatConversionError, ParameterIntConversionError
from .errors import ParameterNeedsFlagError, ParameterNotRecognizedError
from .errors import ParameterNumberOfInputsError, ParameterParseError
from .errors import ParameterRangeError, ParameterUnknownFlagError
from .errors import ParameterValueError
from ._known_parameters import KNOWN_PARAMS
from ._utils import Assignment, NumericBounds, POSITIVE_FLOAT, POSITIVE_INT


logger = logging.getLogger('tleedm.files.parameters')


# Bool parameters for which to create interpret...() methods automatically.
# key is parameter name, value is tuple of keyword arguments that should be
# passed to interpret_bool_parameter, in order.
_SIMPLE_BOOL_PARAMS = {
    'LOG_SEARCH' : (),
    'PHASESHIFTS_CALC_OLD' : (),
    'PHASESHIFTS_OUT_OLD' : (),
    'R_FACTOR_LEGACY' : (),
    'STOP': (),
    'SUPPRESS_EXECUTION' : (),
    'SYMMETRIZE_INPUT' : (),
    'SYMMETRY_FIND_ORI' : (),
    'TL_IGNORE_CHECKSUM' : (),
    'LAYER_STACK_VERTICAL' : ({False: 'c', True: 'z'},),
    }


# Numerical parameters for which to create interpret...() methods
# automatically. Key is parameter name, value is a NumericBounds
_SIMPLE_NUMERICAL_PARAMS = {
    # Positive-only integers
    'BULKDOUBLING_MAX' : POSITIVE_INT,
    'N_CORES' : POSITIVE_INT,
    'SEARCH_MAX_GEN' : POSITIVE_INT,
    'TENSOR_INDEX' : POSITIVE_INT,
    # Positive-only floats
    'T_DEBYE' : POSITIVE_FLOAT,
    'T_EXPERIMENT' : POSITIVE_FLOAT,
    'V0_IMAG' : POSITIVE_FLOAT,
    'TL_VERSION' : POSITIVE_FLOAT,
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


# parameters that can be optimized in FD optimization
_OPTIMIZE_OPTIONS = {'theta', 'phi', 'v0i', 'a', 'b', 'c', 'ab', 'abc',}


def interpretPARAMETERS(rpars, slab=None, silent=False):
    """Interpret rpars.readParams to actual values.

    Parameters
    ----------
    rpars : Rparams
        Object storing parameters for current run. Created previously
        by parameters.read, and should already contain raw string data.
    slab : Slab, optional
        Slab object with elements and atomic position data. If not
        passed, some parameters will not be interpreted.
    silent : bool, optional
        If True, less output will be printed. The default is False.

    Raises
    ------
    ParameterError
        If any parameter interpretation fails.
    """
    interpreter = ParameterInterpreter(rpars)
    interpreter.interpret(slab=slab, silent=silent)


_BELOW_DEBUG = 2


class ParameterInterpreter:                                                     # TODO: order alphabetically, self-tear-down
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
            # Check if we are doing a domain calculation
            _is_domain_calc = 4 in self.rpars.RUN or self.rpars.domainParams
            if _is_domain_calc and param in self.domains_ignore_params:
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
            (param_name, assignment)
            for param_name in self.param_names
            for assignment in self.rpars.readParams[param_name]
            )
        yield from flat_params

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
        ordered_params = 'LOG_LEVEL', 'RUN'
        self.param_names = [p for p in ordered_params
                            if p in self.rpars.readParams]
        self.param_names.extend(
            p for p in KNOWN_PARAMS
            if p in self.rpars.readParams and p not in self.param_names
            )

    # ---------- Methods for interpreting simple parameters -----------
    def interpret_bool_parameter(self, assignment,
                                 allowed_values=None,
                                 param=None, return_only=False,
                                 no_flags=False):
        """Set a parameter to a boolean value.

        Parameters
        ----------
        assignment : Assignment
            The assignment of the parameter, read from PARAMETERS.
        allowed_values : dict, optional
            Additional string values which should be interpreted as
            False or True. E.g. by default, 'f' and 'false' are False,
            't' and 'true' are True. Keys should be True and False,
            values are Sequences of acceptable strings.
        param : str, optional
            The name of the parameter in PARAMETERS, if it differs
            from assignment.parameter. Also used as attribute name
            for self.rpars. Default is None.
        return_only: bool, optional
            If True, only return the value of the parameter, but do
            not set it. Default is False.
        no_flags : bool, optional
            If True, complain if assignment has flags. Default is
            False.

        Returns
        -------
        value : bool
            The value of the parameter.

        Raises
        ------
        ValueError
            If allowed_values contains keys other than True and False.
        ValueError
            If allowed_values contains identical aliases for True
            and False.
        ParameterUnknownFlagError
            When assignment has flags but no_flag is True.
        ParameterBooleanConversionError
            If the string value of the parameter does not
            have an acceptable boolean correspondent.
        """
        if allowed_values is None:
            allowed_values = {}
        if no_flags:
            self._ensure_no_flags_assignment(param, assignment)

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
            raise ParameterBooleanConversionError(assignment.parameter,
                                                  assignment.value) from None

        param = param or assignment.parameter
        if not return_only:
            setattr(self.rpars, param.upper(), value)
        return value

    def interpret_numerical_parameter(self, assignment,
                                      param=None, return_only=False,
                                      bounds=NumericBounds(),                   # TODO: ideally one could default to actually using the limits known from self.rpars.get_limits!
                                      no_flags=False):
        """Set a parameter to a numeric (int or float) value.

        Parameters
        ----------
        assignment : Assignment
            The assignment of the parameter, read from PARAMETERS.
        param : str, optional
            The name of the parameter in PARAMETERS, if it differs
            from assignment.parameter. Also used as attribute name
            for self.rpars. Default is None.
        return_only: bool, optional
            If True, only return the value of the parameter, but
            do not set it. Default is False.
        bounds : NumericBounds, optional
            Acceptable limits for value, and what to do in case an
            out-of-bounds value is given. Default is an unbounded
            float.
        no_flags : bool, optional
            If True, complain if assignment has flags. Default is
            False.

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
        ParameterUnknownFlagError
            When assignment has flags but no_flag is True.
        """
        if no_flags:
            self._ensure_no_flags_assignment(param, assignment)
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
            raise exc(parameter=assignment.parameter) from None

        in_range = all(bounds.is_in_range(value))
        if not in_range and bounds.fail:
            self.rpars.setHaltingLevel(1)
            out_of_range = bounds.format_out_of_range(value)
            raise ParameterRangeError(assignment.parameter,
                                      message=out_of_range)

        if not in_range:
            in_range_value = bounds.make_in_range(value)
            out_of_range = bounds.format_out_of_range(value)
            logger.warning(f'PARAMETERS file: {assignment.parameter}: '
                           f'{out_of_range}. '
                           f'Value will be set to {in_range_value}.')
            value = in_range_value

        param = param or assignment.parameter
        if not return_only:
            setattr(self.rpars, param.upper(), value)
        return value

    @classmethod
    def _make_methods(cls, wrapped, kwargs_names, new_methods_info):
        """Dynamically add methods for this class.

        Parameters
        ----------
        wrapped : callable
            The callable that will be wrapped to create methods.
            The a call to new_method(*args, **kwargs) becomes a
            wrapped(*args, **wrapped_kwargs, **kwargs) call,
            where wrapped_kwargs is generated here using
            kwargs_names and new_methods_info.
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
            kwargs['param'] = param
            kwargs['no_flags'] = True  # no flags for simple parameters
            method = partialmethod(wrapped, **kwargs)
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
    def _ensure_single_flag_assignment(self, param, assignment, message=None):
        """Raise if assignment contains multiple flags."""
        if message is None:
            message = assignment.flags_str
        if assignment.other_flags:
            self.rpars.setHaltingLevel(1)
            raise ParameterUnknownFlagError(param, message)

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
            raise ParameterUnknownFlagError(param, assignment.flag)             # TODO: too-many-flags?

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

    def interpret_bulk_repeat(self, assignment):
        """Assign parameter BULK_REPEAT."""
        param = 'BULK_REPEAT'
        if not self.slab:  # BULK_REPEAT is moot without a slab
            raise ParameterError(parameter=param,
                                 message='No slab defined for bulk repeat.')

        bulk_repeat_str = assignment.values_str.lower()
        # (1) Vector
        if '[' in bulk_repeat_str:
            vec = readVector(bulk_repeat_str, self.slab.ucell)
            if vec is None:
                raise ParameterParseError(parameter=param)
            self.rpars.BULK_REPEAT = vec
            return

        # (2) Z distance
        if '(' not in bulk_repeat_str:
            try:
                self.rpars.BULK_REPEAT = abs(float(assignment.value))
            except ValueError:
                raise ParameterFloatConversionError(parameter=param) from None
            return

        # (3) C or Z distance. Should match, e.g., c(2.0) or z(2.0)
        match = re.match(r'\s*(c|z)\(\s*(?P<val>[0-9.]+)\s*\)',
                         bulk_repeat_str)
        if not match:
            raise ParameterParseError(parameter=param)
        try:
            val = abs(float(match.group('val')))
        except ValueError:
            raise ParameterFloatConversionError(parameter=param) from None

        val = self.slab.ucell[2, 2] * val if 'c' in bulk_repeat_str else val
        self.rpars.BULK_REPEAT = val

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
        domain_step = self.interpret_numerical_parameter(
            assignment,
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
        self._ensure_chemical_elements(param, assignment.values)

        element = assignment.flag.capitalize()
        self.rpars.ELEMENT_MIX[element] = [el.capitalize()
                                           for el in assignment.values]

    def interpret_element_rename(self, assignment):
        """Assign parameter ELEMENT_RENAME."""
        param = 'ELEMENT_RENAME'
        self._ensure_single_flag_and_value_assignment(param, assignment)
        self._ensure_chemical_elements(param, assignment.values)

        self.rpars.ELEMENT_RENAME[assignment.flag.capitalize()] = (
            assignment.value.capitalize()
            )

    def interpret_filament_wf(self, assignment):
        """Assign parameter FILAMENT_WF."""
        param = 'FILAMENT_WF'
        self._ensure_simple_assignment(param, assignment)
        # Check common filaments (e.g., W), otherwise assume a float
        known_filaments = self.rpars.get_default(param)
        try:
            self.rpars.FILAMENT_WF = known_filaments[assignment.value.lower()]
        except KeyError:
            self.interpret_numerical_parameter(assignment)                      # TODO: bounds? We probably want a POSITIVE_FLOAT

    def interpret_fortran_comp(self, assignment, skip_check=False):             # TODO: would be nicer to have a namedtuple or dataclass or similar. It could then have .pre, .post, .mpi, etc...
        """Assign parameter FORTRAN_COMP."""
        param = 'FORTRAN_COMP'
        message = (f'Only one flag allowed for {param} per line. '
                   f'Got {assignment.flags}.')
        self._ensure_single_flag_assignment(param, assignment, message=message)
        flag, compiler_str = assignment.flag.lower(), assignment.values_str

        # (1) Default (i.e., non-MPI) compiler flags
        if not flag and compiler_str.lower() in ['ifort', 'gfortran']:
            self.rpars.getFortranComp(comp=compiler_str.lower(),
                                      skip_check=skip_check)
            return

        # (2) MPI compiler flags
        if flag == 'mpi' and compiler_str.lower() in ['mpifort', 'mpiifort']:
            self.rpars.getFortranMpiComp(comp=compiler_str.lower(),
                                         skip_check=skip_check)
            return

        # (3) Custom compiler flags or full compilation string.
        #     Both need quotation marks                                         # TODO: is this needed now that we use f-strings?
        if not compiler_str.startswith(("'", '"')):
            message = ('No valid shorthand and not delimited '
                       'by quotation marks.')
            raise ParameterError(param, message)

        delim = assignment.values_str[0]
        _, compiler_str, _ = assignment.values_str.split(delim, maxsplit=2)
        if not flag:
            self.rpars.FORTRAN_COMP[0] = compiler_str
        elif flag == 'post':
            self.rpars.FORTRAN_COMP[1] = compiler_str
        elif flag == 'mpi':
            self.rpars.FORTRAN_COMP_MPI[0] = compiler_str
        elif flag == 'mpipost':
            self.rpars.FORTRAN_COMP_MPI[1] = compiler_str
        else:
            raise ParameterUnknownFlagError(parameter=param,
                                            flag=assignment.flag)

    def interpret_intpol_deg(self, assignment):
        """Assign parameter INTPOL_DEG."""
        param = 'INTPOL_DEG'
        self._ensure_simple_assignment(param, assignment)
        intpol_deg = assignment.value
        if intpol_deg in self.rpars.get_limits(param):
            self.rpars.INTPOL_DEG = int(intpol_deg)
            return
        self.rpars.setHaltingLevel(1)
        message = 'Only degree 3 and 5 interpolation supported at the moment.'
        raise ParameterError(parameter=param, message=message)

    def interpret_iv_shift_range(self, assignment):                             # TODO: would be very convenient to have a simple EnergyRange (namedtuple or dataclass) to use for this and THEO_ENERGIES. Then we could have .start, .stop, .step instead of indices.
        """Assign parameter IV_SHIFT_RANGE."""
        param = 'IV_SHIFT_RANGE'
        if len(assignment.values) not in (2, 3):
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(parameter=param)
        iv_range = self._parse_energy_range(param, assignment,
                                            assignment.values,
                                            accept_underscore=True)
        # Interpret underscores as defaults
        _no_value = self.rpars.no_value
        _defaults = self.rpars.get_default(param)
        for i, (bound, default) in enumerate(zip(iv_range, _defaults)):
            if bound is _no_value:
                iv_range[i] = default

        if len(iv_range) == 3:
            start, stop, step = iv_range
        else:
            start, stop, step = *iv_range, _no_value

        if step is not _no_value and (stop - start) * step < 0:
            message = (f'Inconsistent {param} step. Cannot shift from '
                       f'{start:.2f} to {stop:.2f} with {step=:.2f}')
            self.rpars.setHaltingLevel(1)
            raise ParameterError(param, message)

        if stop < start:
            start, stop = stop, start
            step = step if step is _no_value else -step

        self.rpars.IV_SHIFT_RANGE = [start, stop, step]

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
            raise ParameterIntConversionError(param) from None
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

    def interpret_log_level(self, assignment):
        """Assign parameter LOG_LEVEL."""
        param = 'LOG_LEVEL'
        self._ensure_simple_assignment(param, assignment)

        # Try to interpret as bool. This is the way
        # the previous LOG_DEBUG used to work
        try:
            log_debug = self.interpret_bool_parameter(assignment,
                                                      return_only=True)
        except ParameterError:
            pass
        else:
            self.rpars.LOG_LEVEL = logging.DEBUG if log_debug else logging.INFO
            return

        # Try to interpret as string (e.g. 'verbose', 'vv')
        _defaults = self.rpars.get_default(param)
        if assignment.values_str.lower() in _defaults:
            self.rpars.LOG_LEVEL = _defaults[assignment.values_str.lower()]
            return

        # Otherwise interpret as int
        try:
            log_level = self.interpret_numerical_parameter(
                assignment,
                bounds=NumericBounds(type_=int, range_=(0, 50)),
                return_only=True
                )
        except ParameterError as exc:
            raise ParameterValueError(param) from exc
        self.rpars.LOG_LEVEL = log_level

    def interpret_optimize(self, assignment):
        param = 'OPTIMIZE'
        if not assignment.flag:
            message = 'Parameter to optimize not defined.'
            self.rpars.setHaltingLevel(3)
            raise ParameterNeedsFlagError(param, message)
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
            value_error = (f'Value {sl[1]} is not valid for flag {sl[0]}. '
                           'Value will be ignored')
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
                raise ParameterFloatConversionError(param) from None
        if 0 < ps_eps < 1:
            self.rpars.PHASESHIFT_EPS = ps_eps
        else:
            self.rpars.setHaltingLevel(1)
            raise ParameterRangeError(param, given_value=ps_eps,
                                      allowed_range=(0, 1))

    def interpret_plot_iv(self, assignment):
        param = 'PLOT_IV'
        # there should be a flag
        if not assignment.flags:
            self.rpars.setHaltingLevel(1)
            raise ParameterNeedsFlagError(param)

        flag = assignment.flag.lower()
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
                    raise ParameterIntConversionError(param, value) from None
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

    def interpret_run(self, assignment):
        """Assign parameter RUN, inserting an initialization if needed."""
        param = 'RUN'
        segments = []
        for section_str in assignment.values:
            try:
                segments.extend(Section.sequence_from_string(section_str))
            except ValueError as exc:
                self.rpars.setHaltingLevel(2)
                raise ParameterValueError(param, section_str) from exc
        if Section.DOMAINS in segments:
            logger.info('Found domain search.')
        if not segments:
            self.rpars.setHaltingLevel(3)
            raise ParameterError(param,
                                 'RUN was defined, but no values were read.')
        # Insert initialization section if not present
        if segments[0] is not Section.INITIALIZATION:
            segments.insert(0, Section.INITIALIZATION)
        self.rpars.RUN = [s.value for s in segments]                            # TODO: replace with "segments" to keep Section objects

    def interpret_search_beams(self, assignment):
        """Assign parameter SEARCH_BEAMS."""
        param = 'SEARCH_BEAMS'
        self._ensure_simple_assignment(param, assignment)
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
            not is_updating):                                                   # TODO: this is the behaviour of parameters.update(). Was skipping this intended there?
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

    def interpret_search_cull(self, assignment):                                # TODO: custom class to merge the value and type, probably a float subclass (using __new__)
        """Assign parameter SEARCH_CULL and SEARCH_CULL_TYPE."""
        param = 'SEARCH_CULL'
        rpars = self.rpars
        cull_value = assignment.value
        try:
            cull_float = float(cull_value)
        except ValueError as exc:
            rpars.setHaltingLevel(1)
            raise ParameterFloatConversionError(param, cull_value) from exc
        cull_int = int(cull_float)
        if cull_float >= 1 and abs(cull_float - cull_int) < 1e-6:
            rpars.SEARCH_CULL = int(cull_float)
        elif cull_float >= 1:
            message = 'Values greater than one must be integers'
            rpars.setHaltingLevel(1)
            raise ParameterValueError(param, message)
        elif cull_float >= 0:
            rpars.SEARCH_CULL = cull_float
        else:
            message = f'{param} value must be non-negative'
            rpars.setHaltingLevel(1)
            raise ParameterValueError(param, message)

        if not assignment.other_values:
            rpars.SEARCH_CULL_TYPE = rpars.get_default('SEARCH_CULL_TYPE')
            return
        if len(assignment.other_values) > 1:
            rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param)
        cull_type = assignment.other_values[0].lower()
        if cull_type in ['clone', 'genetic', 'random']:
            rpars.SEARCH_CULL_TYPE = cull_type
        else:
            rpars.setHaltingLevel(1)
            raise ParameterValueError(param, cull_type)

    def interpret_search_population(self, assignment):
        param = 'SEARCH_POPULATION'
        self.interpret_numerical_parameter(
            assignment,
            bounds=NumericBounds(type_=int, range_=(1, None))
            )
        if self.rpars.SEARCH_POPULATION < 16:
            logger.warning(f'{param} is very small. A minimum '
                           'value of 16 is recommended.')

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
                            atnums.append(at.num)
                            n -= 1
                else:
                    self.rpars.setHaltingLevel(3)
                    raise ParameterError(param, 'Problem with SITE_DEF input format')
            newdict[sl[0]] = atnums
        self.rpars.SITE_DEF[assignment.flags[0]] = newdict

    def interpret_superlattice(self, assignment):
        """Assign parameter SUPERLATTICE (Wood's or matrix notation)."""
        param = 'SUPERLATTICE'
        self._interpret_wood_or_matrix(param, assignment)
        self.rpars.superlattice_defined = True

    def interpret_symmetry_bulk(self, assignment):
        """Assign parameter SYMMETRY_BULK."""
        param = 'SYMMETRY_BULK'

        # We accept mirrors with syntax "m[a b]", rotations with
        # syntax "rN" and plane groups. Here we only check that
        # the syntax is OK. We defer the complaints about invalid
        # plane group, directions, or rotation orders, to when we
        # know more about the slab
        self.rpars.SYMMETRY_BULK = {  # 'group' added below
            'mirror': set(),
            'rotation': set()
            }
        unrecognized = assignment.values_str

        _mirror_re = re.compile(r'(\s+m\[\s*(-?[012])\s*,?\s*(-?[012])\])',
                                re.IGNORECASE)
        for token, *directions in _mirror_re.findall(unrecognized):
            # All matches of mirrors are acceptable
            unrecognized = unrecognized.replace(token, '')
            direction = tuple(int(v) for v in directions)
            if direction[0] < 0:
                direction = -direction[0], -direction[1]
            self.rpars.SYMMETRY_BULK['mirror'].add(direction)

        _rotation_re = re.compile(r'(\s+r([2346]))', re.IGNORECASE)
        for token, order in _rotation_re.findall(unrecognized):
            unrecognized = unrecognized.replace(token, '')
            self.rpars.SYMMETRY_BULK['rotation'].add(int(order))

        _group_re = re.compile(                                                 # TODO: For now borrowed from guilib. Eventually will try to instantiate a PlaneGroup
            r'(\s*(\w+)\s*(?:\[\s*-?[012]\s*-?[012]\s*\])?)',
            re.IGNORECASE
            )
        for token, group in _group_re.findall(unrecognized):
            if group not in self.grouplist:
                continue
            unrecognized = unrecognized.replace(token, '')
            if 'group' in self.rpars.SYMMETRY_BULK:
                message = 'Only one symmetry group can be given.'
                raise ParameterValueError(param, message=message)
            self.rpars.SYMMETRY_BULK['group'] = token.strip().lower()

        if 'group' not in self.rpars.SYMMETRY_BULK:
            message = 'Need to specify exactly one symmetry group.'
            raise ParameterValueError(param, message=message)

        if unrecognized:
            self.rpars.setHaltingLevel(2)
            message = f'Could not recognize {unrecognized!r}: '
            if 'm' in unrecognized:
                message += ('Syntax for mirrors is "m[n1 n2]"; '
                            'n1 and n2 must be 0, +-1, or +-2; ')
            if 'r' in unrecognized:
                message += ('Syntax for rotations is rN; N must be '
                            '2, 3, 4, or 6')
            raise ParameterValueError(param, message=message)

    def interpret_symmetry_cell_transform(self, assignment):
        """Assign parameter SYMMETRY_CELL_TRANSFORM (Wood's or matrix)."""
        param = 'SYMMETRY_CELL_TRANSFORM'
        self._interpret_wood_or_matrix(param, assignment)

    def interpret_symmetry_eps(self, assignment):
        """Assign parameters SYMMETRY_EPS/SYMMETRY_EPS_Z."""
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
        self.interpret_numerical_parameter(assignment, bounds=bounds)
        if self.rpars.SYMMETRY_EPS > 1.0:
            # pylint: disable-next=logging-format-interpolation
            logger.warning(warning_str.format(''))
        # interpret possible second value as SYMMETRY_EPS_Z
        if not assignment.other_values:
            self.rpars.SYMMETRY_EPS_Z = self.rpars.SYMMETRY_EPS                 # TODO: not nice, as changes to _EPS that may occur are not reflected on _Z. Could do this with a @property of Rparams or a dedicated class.
            return

        z_assignment = Assignment(assignment.other_values, param)
        self.interpret_numerical_parameter(z_assignment,
                                           param='SYMMETRY_EPS_Z',
                                           bounds=bounds)
        if self.rpars.SYMMETRY_EPS_Z > 1.0:
            # pylint: disable-next=logging-format-interpolation
            logger.warning(warning_str.format('for z '))

    def interpret_symmetry_fix(self, assignment):                               # TODO: use symmetry groups from elsewhere once symmetry and guilib are merged
        param = 'SYMMETRY_FIX'
        group = assignment.values_str.lower()
        if group.startswith('t'):  # determine symmetry automatically
            self.rpars.SYMMETRY_FIX = self.rpars.get_default(param)
            return
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

    def interpret_tensor_output(self, assignment):
        """Assign parameter TENSOR_OUTPUT."""
        param = 'TENSOR_OUTPUT'
        tensor_output = []
        # Acceptable format for each token is (repeats*)value,
        # with the bracketed part as optional and 'value' either
        # zero or one (or a bool that can be translated to those)
        rgx = re.compile(r'\s*((?P<repeats>\d+)\s*[*])?\s*(?P<value>[01])\s*$')
        tokens = recombineListElements(assignment.values, '*')
        for token in tokens:
            token_01 = re.sub(r'[Tt](rue|RUE)?', '1', token)
            token_01 = re.sub(r'[Ff](alse|ALSE)?', '0', token_01)

            _match = rgx.match(token_01)
            if not _match:
                self.rpars.setHaltingLevel(1)
                raise ParameterParseError(param,
                                          f'Invalid format in: {token!r}')
            repeats = int(_match['repeats'] or '1')
            value = int(_match['value'])
            tensor_output.extend((value,)*repeats)
        self.rpars.TENSOR_OUTPUT = tensor_output

    def interpret_theo_energies(self, assignment):
        """Assign parameter THEO_ENERGIES. Correct start if inappropriate."""
        param = 'THEO_ENERGIES'
        energies = assignment.values

        # (1) Single value input: treat as default if it is a single
        # underscore, otherwise as a request for a single energy
        if not assignment.other_values and assignment.value == '_':
            self.rpars.THEO_ENERGIES = self.rpars.get_default(param)
            return
        if not assignment.other_values:
            energies = (assignment.value, assignment.value, '1')

        # (2) Three values. Any can be an "_" meaning "default".
        # Internally, we store those as an Rparams.no_value. The
        # others should be positive floats.
        theo_energies = self._parse_energy_range(param, assignment, energies,
                                                 accept_underscore=True)
        non_defaults = [v for v, s in zip(theo_energies, energies) if s != '_']
        if not all(e > 0 for e in non_defaults):
            message = f'{param} values have to be positive.'
            self.rpars.setHaltingLevel(1)
            raise ParameterRangeError(param, message)

        if len(theo_energies) != 3:
            raise ParameterNumberOfInputsError(param)
        if len(non_defaults) < 3:
            self.rpars.THEO_ENERGIES = theo_energies
            return

        start, stop, step = theo_energies
        if stop < start:
            message = f'maximum {param} value should be >= than the minimum.'
            self.rpars.setHaltingLevel(1)
            raise ParameterValueError(param)

        # If the max is not hit by the steps exactly, correct
        # the lower bound down so that this is the case
        # pylint: disable=compare-to-zero  # Clearer this way
        if (stop - start) % step == 0:
            self.rpars.THEO_ENERGIES = theo_energies
            return
        # pylint: enable=compare-to-zero

        start -= step - (stop - start) % step
        if start < 0:
            start = start % step
        if not start:
            start = step
        logger.info('THEO_ENERGIES parameter: (Eto - Efrom) % Estep != 0, '
                    f'Efrom was corrected to {start}')
        self.rpars.THEO_ENERGIES = [start, stop, step]

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
                self.rpars.V0_REAL = [float(c) for c in rundgren_constants]
            except ValueError as exc:
                message = (f'Could not parse constants {rundgren_constants} '
                           'for Rundgren-type function.')
                self.rpars.setHaltingLevel(1)
                raise ParameterError(param, message=message) from exc
            return

        # Pass a specific function to FORTRAN, but replace
        # 'EE' with 'EEV+workfn' to get the energies right
        v0_real = re.sub('(?i)EE', 'EEV+workfn', assignment.values_str)
        self.rpars.V0_REAL = v0_real.rstrip()

    def interpret_vibr_amp_scale(self, assignment):
        """Store raw values: Requires sites for interpreting."""
        self.rpars.VIBR_AMP_SCALE.extend(assignment.values_str.split(','))

    def _interpret_wood_or_matrix(self, param, assignment):
        """Store a Wood- or matrix-notation parameter value."""
        if 'M' not in assignment.flags:
            self._ensure_no_flags_assignment(param, assignment)
            wood = self._read_woods_notation(param, assignment)
            setattr(self.rpars, param, wood)
        else:
            matrix = self._read_matrix_notation(param, assignment.values)
            setattr(self.rpars, param, matrix)

    # ------------------------ Helper methods -------------------------
    def _ensure_chemical_elements(self, param, elements):
        """Raise unless all entries are valid chemical elements."""
        invalid = []
        for element in elements:
            try:
                _ = periodic_table.get_atomic_number(element)
            except ValueError:
                invalid.append(element)

        if invalid:
            message = f'Element(s) {invalid} not found in periodic table.'
            self.rpars.setHaltingLevel(2)
            raise ParameterError(parameter=param, message=message)

    def _parse_energy_range(self, param, assignment,
                            energies, accept_underscore=True):
        """Return a tuple of floats for energies.

        Parameters
        ----------
        param : str
            The parameter to be assigned. Used only for error reporting
        assignment : Assignment
            The assignment for the parameter, Used only for error
            reporting.
        energies : Sequence
            The energies that should be parsed. Elements are
            strings. If accept_underscore is True, they may
            also contain single underscore characters, which
            will be replaced with an Rparams.NO_VALUE.
        accept_underscore : bool, optional
            Whether the energy range should accept underscores in
            the fields to signify "default value"

        Returns
        -------
        float_energies : list
            The energies as a list of floats

        Raises
        ------
        ParameterFloatConversionError
            If conversion to float fails for some of the inputs
        ParameterParseError
            If other parsing errors occur (typically a malformed
            input that cannot be converted to a simple tuple of
            of numbers).
        ParameterRangeError
            In case there are at least 3 entries and the last
            entry (interpreted as a step size) is zero.
        """
        energies_str = ','.join(energies)
        if accept_underscore:
            # We have to replace '_' with something that AST can
            # handle. 'None' seems easy enough. We replace it
            # again further down when converting the rest to float
            energies_str = energies_str.replace('_', 'None')
        try:
            enrange = [float(v) if v is not None else self.rpars.no_value
                       for v in ast.literal_eval(energies_str)]
        except (SyntaxError, ValueError, TypeError,
                MemoryError, RecursionError) as exc:
            self.rpars.setHaltingLevel(1)
            if 'float' in exc.args[0]:
                new_exc = ParameterFloatConversionError
            else:  # Some other weird input
                new_exc = ParameterParseError
            raise new_exc(param, assignment.values_str) from exc

        # Make sure step, the last one, is NOT ZERO
        if len(enrange) < 3:
            return enrange
        try:
            1 / enrange[-1]
        except ZeroDivisionError:
            message = 'Step cannot be zero'
            raise ParameterRangeError(param, message=message) from None
        except TypeError:  # It's a NO_VALUE
            pass
        return enrange

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
                            Assignment(sl[1], param),
                            param=f'{param} {name}',
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
                    Assignment(value, param),
                    param=f'{param} {name}',
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
        """Return a Woods notation from an Assignment, if slab exists."""
        if self.slab is None:
            message = (f'{param} parameter appears to be in Wood '
                       'notation but no slab was passed. Cannot '
                       'calculate bulk unit cell!')
            self.rpars.setHaltingLevel(2)
            raise ParameterError(param, message)
        return readWoodsNotation(assignment.values_str, self.slab.ucell)

    def _read_matrix_notation(self, param, values):
        """Try interpreting values as a 2x2 matrix."""
        matrix = splitSublists(values, ',')
        if len(matrix) != 2:
            message = 'Number of lines in matrix is not equal to 2'
            self.rpars.setHaltingLevel(2)
            raise ParameterParseError(param, message)
        if any(len(row) != 2 for row in matrix):
            message = 'Number of columns in matrix is not equal to 2'
            self.rpars.setHaltingLevel(2)
            raise ParameterParseError(param, message)
        try:
            return np.array(matrix).astype(float)
        except ValueError:
            self.rpars.setHaltingLevel(1)
            raise ParameterFloatConversionError(param, matrix) from None


# Dynamically produce methods for the 'simple parameters' listed above
# Notice that we should do this only once right here, not for each
# instance separately. Calling the methods again is possible, and
# would simply override the bindings.
ParameterInterpreter.make_boolean_interpreter_methods()
ParameterInterpreter.make_numerical_interpreter_methods()
