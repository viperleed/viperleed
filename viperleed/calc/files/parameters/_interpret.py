"""Module _interpret of viperleed.calc.files.parameters.

Initial version by @fkraushofer in 2020, major rewrite by @amimre
and @michele-riva in June 2023. This module used to be part of
parameters.py. Refactored in October 2023.

Functions and classes for interpreting the contents previously
read from a PARAMETERS file.
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
import copy
from functools import partialmethod, wraps
import logging
from pathlib import Path
import re

import numpy as np
try:
    from matplotlib.colors import is_color_like  # For PLOT_IV
except ImportError:
    _CAN_PLOT = False
else:
    _CAN_PLOT = True

from viperleed import __version__
from viperleed.calc.lib import periodic_table
from viperleed.calc.lib.base import readIntRange, readVector
from viperleed.calc.lib.base import recombineListElements, splitSublists
from viperleed.calc.classes.rparams import EnergyRange
from viperleed.calc.classes.rparams import IVShiftRange
from viperleed.calc.classes.rparams import LayerCuts
from viperleed.calc.classes.rparams import SymmetryEps
from viperleed.calc.classes.rparams import TheoEnergies
from viperleed.calc.lib.woods_notation import readWoodsNotation
from viperleed.calc.sections.calc_section import CalcSection as Section

from .errors import ParameterError
from .errors import ParameterBooleanConversionError
from .errors import ParameterFloatConversionError
from .errors import ParameterHasNoValueError
from .errors import ParameterIntConversionError
from .errors import ParameterNeedsFlagError
from .errors import ParameterNeedsSlabError
from .errors import ParameterNotRecognizedError
from .errors import ParameterNumberOfInputsError
from .errors import ParameterParseError
from .errors import ParameterRangeError
from .errors import ParameterUnknownFlagError
from .errors import ParameterValueError
from .errors import SuperfluousParameterError
from ._checker import ParametersChecker
from ._known_parameters import KNOWN_PARAMS, is_deprecated, warn_if_deprecated
from ._utils import Assignment, NumericBounds, POSITIVE_FLOAT, POSITIVE_INT


_LOGGER_NAME, _ = __name__.rsplit('.', maxsplit=1)
_LOGGER = logging.getLogger(_LOGGER_NAME)


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
    'BULKDOUBLING_EPS' : NumericBounds(range_=(1e-4, None),
                                       out_of_range_event='coerce'),
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


def interpret(rpars, slab=None, silent=False):
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


def _with_dynamic_methods(cls):
    """Attach dynamically generated methods to the ParameterInterpreter cls."""

    def _interpret_deprecated(interpreter_method, param):
        """Interpret a deprecated PARAMETER."""
        @wraps(interpreter_method)
        def _wrapper_method(self, assignment, **kwargs):
            if warn_if_deprecated(param, __version__):
                self.rpars.setHaltingLevel(3)
            interpreter_method(self, assignment, **kwargs)
        return _wrapper_method

    # Dynamically produce methods for the 'simple parameters' listed
    # above. Notice that we should do this only once for the class,
    # and not for each instance separately. Calling the methods again
    # is possible, and would simply override the bindings
    cls.make_boolean_interpreter_methods()
    cls.make_numerical_interpreter_methods()

    # Now wrap all the interpreter methods for deprecated parameters
    for method_name, method in cls.__dict__.items():
        if not method_name.startswith('interpret_'):
            continue
        *_, param = method_name.split('interpret_')
        if not is_deprecated(param, __version__):
            continue
        setattr(cls, method_name, _interpret_deprecated(method, param))
    return cls


_BELOW_DEBUG = 2


# Disable 'too-many-public-methods' because of the mechanics of a
# ParameterInterpreter: a dispatcher for interpreting PARAMETERS.
# Methods are public, but are very rarely used. Most of the times
# they are used via the .interpret() interface.
@_with_dynamic_methods
class ParameterInterpreter:  # pylint: disable=too-many-public-methods
    """Class to interpret parameters from the PARAMETERS file.

    To add a new parameter, add a method with the name
    'interpret_<param>' and with a single assignment argument.
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
        self.slab = None
        self._param_names = []  # In precedence order

        # Some flags
        self._search_conv_read = False

    def interpret(self, slab, silent=False):
        """Interpret all known parameters using slab."""
        self.slab = slab
        self._search_conv_read = False

        _backup_log_level = _LOGGER.level
        if silent:
            _LOGGER.setLevel(logging.ERROR)

        self._update_param_order()
        for param, assignment in self._get_param_assignments():
            self._complain_if_invalid_param_in_domain_calc(param)
            self._interpret_param(param, assignment)
            _LOGGER.log(_BELOW_DEBUG,
                        f'Successfully interpreted parameter {param}')

        # Make sure there is no inconsistent combination of parameters
        checker = ParametersChecker()
        checker.check_parameter_conflicts(self.rpars)

        # Finally set the log level back to what it was
        _LOGGER.setLevel(_backup_log_level)

    # ----------------  Helper methods for interpret() ----------------
    def _complain_if_invalid_param_in_domain_calc(self, param):
        """Raise if an invalid param is given in a DOMAINS calculation."""
        _is_domain_calc = 4 in self.rpars.RUN or self.rpars.domainParams
        if not _is_domain_calc or param not in self.domains_ignore_params:
            return
        err_ = (
            f'Found parameter {param} in the main PARAMETERS file of a '
            'DOMAIN calculation. This parameter can only be defined in '
            'the PARAMETERS files of the individual structural domains'
            )
        raise SuperfluousParameterError(param, message=err_)

    def _get_param_assignments(self):
        """Yield parameters and assignments for each PARAMETER read."""
        flat_params = (
            (param_name, assignment)
            for param_name in self._param_names
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
        self._param_names = [p for p in ordered_params
                             if p in self.rpars.readParams]
        self._param_names.extend(
            p for p in KNOWN_PARAMS
            if p in self.rpars.readParams and p not in self._param_names
            )

    # ---------- Methods for interpreting simple parameters -----------
    # Disable pylint warning since, while 6 is a bit many, we can't
    # really do much better than this. Merging some in a container
    # does not seem clearer.
    # pylint: disable-next=too-many-arguments
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
            self._ensure_no_flags_assignment(assignment, param)

        _bool_synonyms = copy.deepcopy(self.bool_synonyms)
        for option, values in allowed_values.items():
            try:
                _bool_synonyms[option].update(v.lower() for v in values)
            except KeyError:
                raise ValueError(f'Unexpected option={option} '
                                 'in allowed_values') from None

        # Make sure there is no intersection between the two sets
        if set(_bool_synonyms[True]) & set(_bool_synonyms[False]):
            raise ValueError('The sets of allowed values for '
                             'True and False must not overlap')

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

    # Disable pylint warning since, while 6 is a bit many, we can't
    # really do much better than this. Merging some in a container
    # does not seem clearer.
    # pylint: disable-next=too-many-arguments
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
            self._ensure_no_flags_assignment(assignment, param)
        type_ = bounds.type_
        exc_cls = (ParameterFloatConversionError if type_ is float
                   else ParameterIntConversionError)
        exc = exc_cls(assignment.parameter, assignment.value)

        try:
            # First convert to float. Necessary for, e.g., 1e6 as int
            float_value = float(assignment.value)
        except ValueError:
            self.rpars.setHaltingLevel(1)
            raise exc from None

        value = type_(float_value)
        if type_ is int and not np.isclose(value, float_value):
            raise exc

        in_range = all(bounds.is_in_range(value))
        if not in_range and bounds.fail:
            self.rpars.setHaltingLevel(1)
            out_of_range = bounds.format_out_of_range(value)
            raise ParameterRangeError(assignment.parameter,
                                      message=out_of_range)

        if not in_range:
            in_range_value = bounds.make_in_range(value)
            out_of_range = bounds.format_out_of_range(value)
            _LOGGER.warning(f'PARAMETERS file: {assignment.parameter}: '
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
            The call to new_method(*args, **kwargs) becomes a
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
    def _ensure_single_flag_assignment(self, assignment, param=None,
                                       message=None,
                                       must_have_exaclty_one=True):
        """Raise if assignment contains multiple or, optionally, no flags."""
        if param is None:
            param = assignment.parameter
        if must_have_exaclty_one and not assignment.flag:
            self.rpars.setHaltingLevel(1)
            raise ParameterNeedsFlagError(param, message)
        if message is None:
            message = assignment.flags_str
        if assignment.other_flags:
            self.rpars.setHaltingLevel(1)
            raise ParameterUnknownFlagError(param, message)

    def _ensure_single_value_assignment(self, assignment, param=None):
        """Raise if assignment contains multiple values."""
        n_values = len(assignment.values)
        if param is None:
            param = assignment.parameter
        if not n_values:
            self.rpars.setHaltingLevel(1)
            raise ParameterHasNoValueError(param)
        if n_values != 1:
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param, (n_values, 1))

    def _ensure_no_flags_assignment(self, assignment, param=None):
        """Raise if assignment contains any flag."""
        if param is None:
            param = assignment.parameter
        if assignment.flags:
            self.rpars.setHaltingLevel(1)
            raise ParameterUnknownFlagError(param, assignment.flag)

    def _ensure_simple_assignment(self, assignment, param=None):
        """Raise if assignment is not simple (i.e., one value, no flags)."""
        self._ensure_no_flags_assignment(assignment, param)
        self._ensure_single_value_assignment(assignment, param)

    def _ensure_single_flag_and_value_assignment(self, assignment, param=None):
        """Raise if assignment does not have a single flag--value pair."""
        self._ensure_single_flag_assignment(assignment, param)
        self._ensure_single_value_assignment(assignment, param)

    # ----------- Methods to interpret individual parameters ----------
    def interpret_average_beams(self, assignment):                              # TODO: Issue #135
        """Assign parameter AVERAGE_BEAMS."""
        param = 'AVERAGE_BEAMS'
        self._ensure_no_flags_assignment(assignment)
        right_side = assignment.values_str.lower().strip()
        # trivial cases
        if right_side in {'off', 'none', 'false', 'f'}:
            self.rpars.AVERAGE_BEAMS = False
            return
        if right_side.lower() in {'all', 'perpendicular', 'perp'}:
            self.rpars.AVERAGE_BEAMS = (0., 0.)
            return

        # Otherwise, try to parse the value
        angles = self._parse_incidence_angles(param, assignment)
        self.rpars.AVERAGE_BEAMS = angles['THETA'], angles['PHI']

    def interpret_beam_incidence(self, assignment):
        """Assign parameter BEAM_INCIDENCE."""
        param = 'BEAM_INCIDENCE'
        self._ensure_no_flags_assignment(assignment)
        angles = self._parse_incidence_angles(param, assignment)
        self.rpars.THETA, self.rpars.PHI = angles['THETA'], angles['PHI']

    def interpret_bulk_repeat(self, assignment):
        """Assign parameter BULK_REPEAT."""
        param = 'BULK_REPEAT'
        if not self.slab:  # BULK_REPEAT is moot without a slab
            raise ParameterNeedsSlabError(param)

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

        val = self.slab.c_vector[2] * val if 'c' in bulk_repeat_str else val
        self.rpars.BULK_REPEAT = val

    def interpret_domain(self, assignment):
        """Set the domain path and name."""
        param = 'DOMAIN'
        # Check name
        name = assignment.flag
        if name in self.rpars.DOMAINS:  # Already defined
            error_message = f'Multiple sources defined for domain {name}'
            self.rpars.setHaltingLevel(3)
            raise ParameterValueError(param, message=error_message)
        if not name:  # Get unique name                                         # TODO: used in several other places
            i = 1
            while str(i) in self.rpars.DOMAINS:
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
                             'interpreted as either a path or a .zip file')
            self.rpars.setHaltingLevel(3)
            raise ParameterValueError(param, message=error_message)
        self.rpars.DOMAINS[name] = path

    def interpret_domain_step(self, assignment):
        """Assign parameter DOMAIN_STEP."""
        param = 'DOMAIN_STEP'
        self._ensure_simple_assignment(assignment)
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
                       f'Consider using {j} instead')
            self.rpars.setHaltingLevel(1)
            raise ParameterValueError(param, message=message)
        self.rpars.DOMAIN_STEP = domain_step

    def interpret_element_mix(self, assignment):
        """Assign parameter ELEMENT_MIX."""
        param = 'ELEMENT_MIX'
        self._ensure_single_flag_assignment(assignment)
        element = self._get_valid_slab_element_from_flag(param, assignment)

        self._ensure_chemical_elements(param, assignment.values)
        mix = [el.capitalize() for el in assignment.values]
        if len(mix) == 1:
            rename = mix[0]
            self.rpars.setHaltingLevel(3)
            message = (f'Only one mixed element found: {rename}. If you '
                       f'intended to rename {element} to {rename}, you should '
                       f'use "ELEMENT_RENAME {element} = {rename}"')
            raise ParameterNumberOfInputsError(param, message=message)
        self.rpars.ELEMENT_MIX[element] = mix

    def interpret_element_rename(self, assignment):
        """Assign parameter ELEMENT_RENAME."""
        param = 'ELEMENT_RENAME'
        self._ensure_single_flag_and_value_assignment(assignment)
        element = self._get_valid_slab_element_from_flag(param, assignment)
        self._ensure_chemical_elements(param, assignment.values)
        self.rpars.ELEMENT_RENAME[element] = assignment.value.capitalize()

    def interpret_filament_wf(self, assignment):
        """Assign parameter FILAMENT_WF."""
        param = 'FILAMENT_WF'
        self._ensure_simple_assignment(assignment)
        # Check common filaments (e.g., W), otherwise assume a float
        known_filaments = self.rpars.get_default(param)
        try:
            self.rpars.FILAMENT_WF = known_filaments[assignment.value.lower()]
        except KeyError:
            self.interpret_numerical_parameter(assignment)

    def interpret_fortran_comp(self, assignment):                               # TODO: Issue #136
        """Assign parameter FORTRAN_COMP."""
        param = 'FORTRAN_COMP'
        message = (f'Only one flag allowed for {param} per line. '
                   f'Got {assignment.flags}')
        self._ensure_single_flag_assignment(assignment, message=message,
                                            must_have_exaclty_one=False)
        flag, compiler_str = assignment.flag.lower(), assignment.values_str

        # (1) Default (i.e., non-MPI) compiler flags
        if not flag and compiler_str.lower() in {'ifort', 'gfortran'}:
            self.rpars.FORTRAN_COMP[0] = compiler_str.lower()
            return

        # (2) MPI compiler flags
        if flag == 'mpi' and compiler_str.lower() in {'mpifort', 'mpiifort'}:
            self.rpars.FORTRAN_COMP_MPI[0] = compiler_str.lower()
            return

        # (3) Custom compiler flags or full compilation string.
        #     Remove optional quotes
        delim = assignment.values_str[0]
        if delim in {'"', "'"}:
            compiler_str = assignment.values_str[1:]
            if not compiler_str.endswith(delim):
                raise ParameterValueError(param,
                                          message='No closing quotation mark')
            compiler_str = compiler_str[:-1]

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
        self._ensure_simple_assignment(assignment)
        intpol_deg = assignment.value
        if intpol_deg in self.rpars.get_limits(param):
            self.rpars.INTPOL_DEG = int(intpol_deg)
            return
        self.rpars.setHaltingLevel(1)
        message = 'Only degree 3 and 5 interpolation supported at the moment'
        raise ParameterValueError(param, message=message)

    def interpret_iv_shift_range(self, assignment):
        """Assign parameter IV_SHIFT_RANGE."""
        param = 'IV_SHIFT_RANGE'
        energies = self._parse_energy_range(assignment, assignment.values)
        if len(energies) not in (2, 3):
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(parameter=param)

        try:
            iv_range = IVShiftRange(*energies)
        except (ValueError, RuntimeError) as exc:
            # RuntimeError: was fixed, and we tried to unfix it
            self.rpars.setHaltingLevel(1)
            raise ParameterValueError(param, message=str(exc)) from exc

        # Make sure that the user gave us a reasonable
        # step. Inform them if we had to change stuff
        if iv_range.defined and not iv_range.is_equivalent(energies):
            _LOGGER.info(f'IV_SHIFT_RANGE bounds were not integer multiples '
                         f'of {iv_range.step}. Bounds have been adjusted to '
                         f'{iv_range.min} {iv_range.max}.')

        # Interpret underscores as defaults
        iv_range.set_undefined_values(self.rpars.get_default(param))
        self.rpars.IV_SHIFT_RANGE = iv_range

    def interpret_layer_cuts(self, assignment):
        """Assign parameter LAYER_CUTS."""
        param = 'LAYER_CUTS'
        try:
            cuts = LayerCuts.from_string(assignment.values_str)
        except ValueError as exc:
            # String has invalid syntax
            raise ParameterParseError(param, message=str(exc)) from None
        if not cuts:
            raise ParameterValueError(param,
                                      message='At least one layer cut needed')
        self.rpars.LAYER_CUTS = cuts

    @staticmethod
    def _tokenize_int_range(numeric_string):
        r"""Return integers representing range bounds in `numeric_string`.

        A token of the form '-\d+' is considered as a negative number
        unless it is preceded by another '\d+' without delimiters, or
        only by '-' delimiters. This means that
        >>> _tokenize_int_range('1-15')
        [1, 15]
        >>> _tokenize_int_range('-1---15')
        [-1, 15]
        >>> _tokenize_int_range('-1-- -15')
        [-1, -15]
        >>> _tokenize_int_range('1 -15')
        [1, -15]

        Acceptable delimiters are spaces, '-', or ':'.

        Parameters
        ----------
        numeric_string : str
            The string to be analysed.

        Returns
        -------
        integers : list
            Integers in `numeric_string`.

        Raises
        ------
        ValueError
            If any of the entries is neither a number or a delimiter
        ValueError
            If, after stripping spaces, `numeric_string` ends with a
            delimiter instead of a number
        """
        string = re.sub(r'\s', ' ', numeric_string.strip())
        delimiter_chars = ' :-'
        delimiters = tuple(delimiter_chars)
        # Remove all adjacent identical delimiters
        for delimiter in delimiter_chars:
            string = re.sub(fr'{delimiter}*{delimiter}', delimiter, string)
        if not string:
            return []
        _digits_or_delim_re = re.compile(fr'[\d{delimiter_chars}]+')
        _digits_re = re.compile(r'\d+')
        if not _digits_or_delim_re.fullmatch(string):
            raise ValueError('contains non-integers or invalid delimiters')
        if string.endswith(delimiters):
            raise ValueError('ends with a delimiter')
        integers = []
        for match_digits in _digits_re.finditer(string):
            token, start = match_digits.group(), match_digits.start()
            delimiter = string[start-1:start]
            if not delimiter:  # First token in string
                integers.append(int(token))
                continue
            if delimiter != '-':
                integers.append(int(token))
                continue
            # Figure out if the '-' is a delimiter or it means negative
            previous = string[:max(start-1, 0)]
            if not previous or previous.endswith(delimiters):
                # Negative number
                integers.append(int(delimiter + token))
            else:
                integers.append(int(token))
        return integers

    def interpret_lmax(self, assignment):
        """Assign parameter LMAX."""
        param = 'LMAX'
        try:
            values = self._tokenize_int_range(assignment.values_str)
        except ValueError as exc:
            raise ParameterValueError(
                param,
                message=f'Invalid {assignment.values_str!r}: {exc}'
                ) from None

        if not values:
            self.rpars.setHaltingLevel(1)
            raise ParameterHasNoValueError(param)
        if len(values) > 2:
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param)
        try:
            self.rpars.LMAX = self.rpars.LMAX.from_value(values)
        except ValueError as exc:  # out-of-range
            raise ParameterRangeError(param, message=str(exc))

    def interpret_log_level(self, assignment):
        """Assign parameter LOG_LEVEL."""
        param = 'LOG_LEVEL'
        self._ensure_simple_assignment(assignment)

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
            raise ParameterValueError(param, assignment.values_str) from exc
        self.rpars.LOG_LEVEL = log_level

    def interpret_optimize(self, assignment):
        """Assign parameter OPTIMIZE."""
        param = 'OPTIMIZE'
        if not assignment.flag:
            message = 'Parameter to optimize not defined'
            self.rpars.setHaltingLevel(3)
            raise ParameterNeedsFlagError(param, message)
        which = assignment.flag.lower()
        if which not in _OPTIMIZE_OPTIONS:
            self.rpars.setHaltingLevel(3)
            raise ParameterUnknownFlagError(param, f'{which!r}')
        if Section.FD_OPTIMIZATION.value not in self.rpars.RUN:                 # TODO: remove .value when using TLEEDMSection in RUN
            err_ = ('RUN does not include a full-dynamic-optimization section '
                    f'(RUN = {Section.FD_OPTIMIZATION.value}). It makes no '
                    'sense to provide an OPTIMIZE parameter')
            raise SuperfluousParameterError(param, message=err_)
        self.rpars.OPTIMIZE['which'] = which
        if not assignment.other_values:
            try:
                self.rpars.OPTIMIZE['step'] = float(assignment.value)
            except ValueError:
                self.rpars.setHaltingLevel(1)
                raise ParameterFloatConversionError(param,
                                                    assignment.value) from None
            return
        for flag_value_pair in assignment.values_str.split(','):
            self._interpret_optimize_flag_value_pair(param, flag_value_pair)

    def _interpret_optimize_flag_value_pair(self, param, flag_value_pair):
        """Interpret one 'flag value' pair for OPTIMIZE."""
        flag, value = self._get_flag_value_from_pair(param, flag_value_pair)
        value_error = f'Value {value!r} is invalid for flag {flag!r}'
        partype = {'step': float, 'convergence': float,
                   'minpoints': int, 'maxpoints': int,
                   'maxstep': float}
        try:
            numeric = partype[flag](value)
        except ValueError as exc:
            self.rpars.setHaltingLevel(1)
            raise ParameterValueError(param, message=value_error) from exc
        except KeyError:
            self.rpars.setHaltingLevel(2)
            raise ParameterUnknownFlagError(param, f'{flag!r}') from None
        self.rpars.OPTIMIZE[flag] = numeric

    def interpret_parabola_fit(self, assignment):                               # TODO: Issue #137
        """Assign parameter PARABOLA_FIT."""
        if assignment.values_str.lower() == 'off':
            self.rpars.PARABOLA_FIT['type'] = 'none'
            return
        for flag_value_pair in assignment.values_str.split(','):
            self._interpret_parabola_fit_flag_value_pair(flag_value_pair)

    def _interpret_parabola_fit_flag_value_pair(self, flag_value_pair):
        """Interpret one 'flag value' pair for PARABOLA_FIT."""
        param = 'PARABOLA_FIT'
        flag, value = self._get_flag_value_from_pair(param, flag_value_pair)
        if flag == 'localise':
            flag = 'localize'
        value_error = (f'Invalid value {value} for flag {flag}. '
                       'Value will be ignored')
        if flag not in ('type', 'alpha', 'mincurv', 'localize'):
            self.rpars.setHaltingLevel(1)
            raise ParameterValueError(param, f'Unknown flag={flag!r}')
        if flag == 'type' and value not in ('linear', 'linearregression',
                                            'lasso', 'ridge', 'elasticnet',
                                            'none'):
            self.rpars.setHaltingLevel(1)
            raise ParameterValueError(param, value_error)
        if flag == 'type':
            self.rpars.PARABOLA_FIT[flag] = value
            return
        try:
            value_float = float(value)
        except ValueError:
            self.rpars.setHaltingLevel(1)
            raise ParameterFloatConversionError(param, value) from None
        if value_float >= 0:
            self.rpars.PARABOLA_FIT[flag] = value_float
            return
        self.rpars.setHaltingLevel(1)
        raise ParameterRangeError(param,
                                  message=f'{flag} value must be non-negative')

    def interpret_phaseshift_eps(self, assignment):
        """Assign parameter PHASESHIFT_EPS."""
        param = 'PHASESHIFT_EPS'
        self._ensure_simple_assignment(assignment)
        try:
            ps_eps = float(assignment.value)
        except ValueError:
            # check if one of default values (e.g. 'fine')
            preset = assignment.value.lower()[0]
            _defaults = self.rpars.get_default(param)
            ps_eps = _defaults.get(preset, None)
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
        """Assign parameter PLOT_IV."""
        if not _CAN_PLOT:
            # Cannot interpret this parameter
            return
        param = 'PLOT_IV'
        self._ensure_single_flag_assignment(assignment)
        flag_aliases = {
            'plot': 'plot',
            'axes': 'axes', 'border': 'axes', 'borders': 'axes',
            'colors': 'colors', 'colours': 'colors',
            'color': 'colors', 'colour': 'colors',
            'font': 'font_size', 'fontsize': 'font_size',
            'font_size': 'font_size',
            'legend': 'legend', 'legends': 'legend',
            'lw': 'line_width', 'linewidth': 'line_width',
            'line_width': 'line_width', 'line': 'line_width',
            'overbar': 'overbar', 'overline': 'overbar',
            'perpage': 'perpage', 'layout': 'perpage',
            }
        flag = assignment.flag.lower()
        try:
            flag = flag_aliases[flag]
        except KeyError:
            self.rpars.setHaltingLevel(1)
            raise ParameterUnknownFlagError(param, f'{flag!r}') from None
        setter = getattr(self, f'_interpret_plot_iv__{flag}')
        setter(assignment)

    def _interpret_plot_iv__axes(self, assignment):
        """Assign PLOT_IV['axes']."""
        self._ensure_single_value_assignment(assignment)
        synonyms = {'all': 'all',
                    'none': 'none',
                    'lb': 'lb', 'less': 'lb',
                    'b': 'b', 'bottom': 'b',}
        try:
            self.rpars.PLOT_IV['axes'] = synonyms[assignment.value.lower()]
        except KeyError:
            self.rpars.setHaltingLevel(1)
            raise ParameterParseError(assignment.parameter) from None

    def _interpret_plot_iv__colors(self, assignment):
        """Assign PLOT_IV['colors']."""
        colors = assignment.values
        if not colors:
            self.rpars.setHaltingLevel(1)
            raise ParameterHasNoValueError(assignment.parameter)
        if not all(is_color_like(c) for c in colors):
            err_ = ('Cannot interpret at least one of '
                    f'{assignment.values_str!r} as color(s)')
            raise ParameterValueError(assignment.parameter, message=err_)
        self.rpars.PLOT_IV['colors'] = colors

    def _interpret_plot_iv__font_size(self, assignment):
        """Assign PLOT_IV['font_size']."""
        self._ensure_single_value_assignment(assignment)
        try:
            font_size = float(assignment.value)
        except ValueError:
            self.rpars.setHaltingLevel(1)
            raise ParameterFloatConversionError(assignment.parameter,
                                                assignment.value) from None
        if font_size <= 0:
            self.rpars.setHaltingLevel(1)
            raise ParameterRangeError(assignment.parameter,
                                      message='Font size must be positive.')
        self.rpars.PLOT_IV['font_size'] = font_size

    def _interpret_plot_iv__legend(self, assignment):
        """Assign PLOT_IV['legend']."""
        self._ensure_single_value_assignment(assignment)
        value = assignment.value.lower()
        if value in {'all', 'first', 'none'}:
            self.rpars.PLOT_IV['legend'] = value
        elif value in {'topright', 'tr'}:
            self.rpars.PLOT_IV['legend'] = 'tr'
        else:
            self.rpars.setHaltingLevel(1)
            raise ParameterParseError(assignment.parameter)

    def _interpret_plot_iv__line_width(self, assignment):
        """Assign PLOT_IV['line_width']."""
        self._ensure_single_value_assignment(assignment)
        try:
            line_width = float(assignment.value)
        except ValueError:
            self.rpars.setHaltingLevel(1)
            raise ParameterFloatConversionError(assignment.parameter,
                                               assignment.value) from None
        if line_width <= 0:
            self.rpars.setHaltingLevel(1)
            raise ParameterRangeError(assignment.parameter,
                                      message='Line width must be > 0.')
        self.rpars.PLOT_IV['line_width'] = line_width

    def _interpret_plot_iv__overbar(self, assignment):
        """Assign PLOT_IV['overbar']."""
        self._ensure_single_value_assignment(assignment)
        value = assignment.value.lower()
        if value.startswith('t'):
            self.rpars.PLOT_IV['overbar'] = True
        elif value.startswith('f'):
            self.rpars.PLOT_IV['overbar'] = False
        else:
            message = f'Value for flag {assignment.flag!r} not recognized'
            self.rpars.setHaltingLevel(1)
            raise ParameterParseError(assignment.parameter, message)

    def _interpret_plot_iv__perpage(self, assignment):
        """Assign PLOT_IV['perpage']."""
        # Can be one or two integers
        if not assignment.values:
            self.rpars.setHaltingLevel(1)
            raise ParameterHasNoValueError(assignment.parameter)
        if len(assignment.values) > 2:
            self.rpars.setHaltingLevel(1)
            message = 'perpage accepts one or two positive integers'
            raise ParameterNumberOfInputsError(assignment.parameter,
                                               message=message)
        try:
            values = tuple(int(v) for v in assignment.values)
        except ValueError:
            self.rpars.setHaltingLevel(1)
            raise ParameterIntConversionError(assignment.parameter,
                                              assignment.values) from None
        if any(v <= 0 for v in values):
            _plural = 's' if len(values)>1 else ''
            message = f'perpage value{_plural} must be positive'
            raise ParameterRangeError(assignment.parameter, message=message)
        self.rpars.PLOT_IV['perpage'] = values if len(values)>1 else values[0]

    def _interpret_plot_iv__plot(self, assignment):
        """Assign PLOT_IV['plot']."""
        self._ensure_single_value_assignment(assignment)
        extra_synonyms = {False: ('none',)}
        self.rpars.PLOT_IV['plot'] = self.interpret_bool_parameter(
            assignment, extra_synonyms, return_only=True
            )

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
            _LOGGER.info('Found domain search.')
        if not segments:
            self.rpars.setHaltingLevel(3)
            message = f'{param} was defined, but no values were read'
            raise ParameterHasNoValueError(param, message)
        # Insert initialization section if not present
        if segments[0] is not Section.INITIALIZATION:
            segments.insert(0, Section.INITIALIZATION)
        self.rpars.RUN = [s.value for s in segments]                            # TODO: replace with "segments" to keep Section objects

    def interpret_search_beams(self, assignment):
        """Assign parameter SEARCH_BEAMS."""
        param = 'SEARCH_BEAMS'
        self._ensure_simple_assignment(assignment)
        value = assignment.value.lower()
        if value.startswith(('0', 'a')):
            self.rpars.SEARCH_BEAMS = 0
        elif value.startswith(('1', 'i')):
            self.rpars.SEARCH_BEAMS = 1
        elif value.startswith(('2', 'f')):
            self.rpars.SEARCH_BEAMS = 2
        else:
            raise ParameterValueError(param, value)

    def interpret_search_convergence(self, assignment, is_updating=False):      # TODO: Issue #139
        """Interpret SEARCH_CONVERGENCE parameter.

        This method also sets rpars.searchConvInit if `is_updating`
        is True-thy and the values to be assigned have changed.

        Parameters
        ----------
        assignment : Assignment
            The assignment line, containing details of the
            flags and values to be interpreted.
        is_updating : bool, optional
            Whether this method is being called as part of a
            parameters.update (True) or a parameters.interpret
            (False). Default is False.

        Raises
        ------
        ParameterHasNoValueError
            If assignment carries no values.
        ParameterNumberOfInputsError
            If there are more than two values for this parameter.
        ParameterNeedsFlagError
            If more values are present, but no flag was given.
        ParameterUnknownFlagError
            If the first flag is neither 'dgen' nor 'gaussian'.
        ParameterUnknownFlagError
            If the first flag is 'gaussian' and there are more
            flags.
        ParameterUnknownFlagError
            If first flag is 'dgen' and there is more than one
            additional flag, or if the additional flag is not
            one of the acceptable ones ('dec', 'best', 'all').
        ParameterFloatConversionError
            If numeric values cannot be converted to float.
        ParameterRangeError
            If Gaussian width, its scaling factor, dgen, or its
            scaling factor are out of bounds
        """
        param = 'SEARCH_CONVERGENCE'
        if not assignment.values:
            raise ParameterHasNoValueError(param)
        if len(assignment.values) > 2:
            raise ParameterNumberOfInputsError(param)

        if (not assignment.flags
                and not assignment.other_values
                and assignment.value.lower() == 'off'):
            self.rpars.GAUSSIAN_WIDTH_SCALING = 1.
            return

        if not assignment.flags:
            raise ParameterNeedsFlagError(param)
        if assignment.flag.lower() not in ('dgen', 'gaussian'):
            raise ParameterUnknownFlagError(param, f'{assignment.flag!r}')

        numeric = [None, None]
        for i, value in enumerate(assignment.values):
            try:
                numeric[i] = float(value)
            except ValueError as err:
                raise ParameterFloatConversionError(param) from err

        if assignment.flag.lower() == 'gaussian':
            self._ensure_single_flag_assignment(assignment)
            self._interpret_search_convergence_gaussian(
                param, *numeric, is_updating
                )
        else:
            self._interpret_search_convergence_dgen(
                param, assignment.other_flags, *numeric, is_updating
                )

    # disable: The only option would be to pack max_dgen and scaling
    # into a single object, but I'm not sure this would make it
    # more understandable.
    # pylint: disable-next=too-many-arguments
    def _interpret_search_convergence_dgen(self, param, flags, max_dgen,
                                           scaling, is_updating):
        """Assign the SEARCH_MAX_DGEN/_SCLING parameters."""
        if len(flags) > 1:
            raise ParameterUnknownFlagError(param, ' '.join(flags[1:]))
        target = flags[0] if flags else 'dec'
        if target not in ['dec', 'best', 'all']:
            raise ParameterUnknownFlagError(param, f'{target!r}')
        should_store_new_value = (
            not is_updating
            or max_dgen != self.rpars.searchConvInit['dgen'][target]
            )
        valid = max_dgen is not None and max_dgen > 0
        if should_store_new_value and not valid:
            message = 'dgen should be a positive number'
            raise ParameterRangeError(param, message=message)

        if should_store_new_value and not self._search_conv_read:
            # Clear default values
            self.rpars.reset_default('SEARCH_MAX_DGEN')
            self._search_conv_read = True
        if should_store_new_value:
            self.rpars.SEARCH_MAX_DGEN[target] = max_dgen
        if should_store_new_value and is_updating:
            self.rpars.searchConvInit['dgen'][target] = max_dgen

        if scaling is not None and scaling >= 1:
            self.rpars.SEARCH_MAX_DGEN_SCALING[target] = scaling
        elif scaling:
            message = 'dgen scaling value cannot be smaller than 1'
            raise ParameterRangeError(param, message=message)

    def _interpret_search_convergence_gaussian(self, param, gauss_width,
                                               scaling, is_updating):
        """Assign the GAUSSIAN_WIDTH/_SCALING parameters."""
        should_store_new_value = (
            not is_updating
            or gauss_width != self.rpars.searchConvInit['gaussian']
            )
        valid = gauss_width is not None and gauss_width > 0
        if should_store_new_value and not valid:
            message = 'gaussian width should be a positive number'
            raise ParameterRangeError(param, message=message)

        if should_store_new_value:
            self.rpars.GAUSSIAN_WIDTH = gauss_width
        if should_store_new_value and is_updating:
            self.rpars.searchConvInit['gaussian'] = gauss_width

        if scaling is not None and 0 < scaling <= 1:
            self.rpars.GAUSSIAN_WIDTH_SCALING = scaling
        elif scaling is not None:
            raise ParameterRangeError(param, scaling, (0, 1))

    def interpret_search_cull(self, assignment):
        """Assign parameter SEARCH_CULL."""
        param = 'SEARCH_CULL'
        rpars = self.rpars
        if len(assignment.values) not in (1, 2):
            rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param)
        try:
            cull = rpars.SEARCH_CULL.from_value(assignment.values)
        except TypeError as exc:   # Numeric value is not a number
            rpars.setHaltingLevel(1)
            raise ParameterFloatConversionError(param,
                                                message=str(exc)) from exc
        except ValueError as exc:  # Numeric or type are wrong
            rpars.setHaltingLevel(1)
            raise ParameterValueError(param, message=str(exc)) from exc

        if not cull.has_type:
            _default = rpars.get_default('SEARCH_CULL')
            cull.type_ = _default.type_
        rpars.SEARCH_CULL = cull

    def interpret_search_population(self, assignment):
        param = 'SEARCH_POPULATION'
        self.interpret_numerical_parameter(
            assignment,
            bounds=NumericBounds(type_=int, range_=(1, None))
            )
        if self.rpars.SEARCH_POPULATION < 16:
            _LOGGER.warning(f'{param} is very small. A minimum '
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

    def interpret_site_def(self, assignment):                                   # TODO: Issue #143
        """Assign the SITE_DEF for one POSCAR element."""
        param = 'SITE_DEF'
        if assignment.values_str.count('top(') > 1:
            self.rpars.setHaltingLevel(3)
            raise ParameterValueError(
                param,
                message='only a single top() allowed per SITE_DEF line'
                )
        self._ensure_single_flag_assignment(assignment)
        site_element = self._get_valid_slab_element_from_flag(param,
                                                              assignment)
        sorted_atoms = []
        if 'top(' in assignment.values_str:
            sorted_atoms = sorted(
                (at for at in self.slab if at.el == site_element),
                key=lambda atom: atom.pos[2],
                reverse=True
                )

        site_def_dict = {}
        for flag_and_values in assignment.values_str.strip().split(','):
            site_label, site_specs = self._get_valid_sitedef_spec(
                param, site_element, flag_and_values
                )
            try:
                atnums = self._get_atom_numbers_for_site(site_specs,
                                                         sorted_atoms)
            except ValueError as exc:  # Invalid syntax
                self.rpars.setHaltingLevel(3)
                err_ = f'Invalid syntax in {assignment.values_str!r}: {exc!r}'
                raise ParameterParseError(param, message=err_) from None
            try:
                self._ensure_sitedef_makes_sense(site_element,
                                                 site_label,
                                                 atnums)
            except ValueError as exc:
                self.rpars.setHaltingLevel(3)
                raise ParameterValueError(assignment.parameter,
                                          message=str(exc)) from None
            site_def_dict[site_label] = atnums
        self.rpars.SITE_DEF[site_element] = site_def_dict

    def _get_valid_sitedef_spec(self, param, site_element, label_and_specs):
        """Return a site label and atom selection specifications, or raise."""
        try:
            site_label, *site_specs = (s.strip()
                                       for s in label_and_specs.split())
        except ValueError:
            self.rpars.setHaltingLevel(2)
            err_ = ('Expected "site_label atom_selection_patterns ". '
                    f'Found "{label_and_specs}"')
            raise ParameterNumberOfInputsError(param, message=err_) from None
        if not site_specs:
            self.rpars.setHaltingLevel(2)
            err_ = ('No atom selection pattern found for site '
                    f'label {site_label}, element {site_element}.')
            raise ParameterParseError(param, message=err_)
        return site_label, site_specs

    @staticmethod
    def _get_atom_numbers_for_site(site_specs, sorted_atoms):
        """Return atom numbers given a selection specification.

        Parameters
        ----------
        site_specs : Sequence
            Items are strings. Each corresponds to a specification for
            selecting some atom numbers. Acceptable specifications are:
            - single integer
            - range of integers (e.g., '1-5', '10:28')
            - 'top(N)'
        sorted_atoms : Sequence
            Atom objects, with element already correct for the specific
            site, sorted by decreasing fractional position. Used only
            if one of the `site_specs` is 'top(N)'.

        Returns
        -------
        atnums : list
            Atom numbers that are selected by the given `site_specs`.

        Raises
        ------
        ValueError
            If any of the `site_specs` has invalid syntax.
        """
        atnums = []
        for site_spec in site_specs:
            extra_atom_numbers = readIntRange(site_spec)
            if extra_atom_numbers:
                atnums.extend(extra_atom_numbers)
                continue
            if 'top(' in site_spec:
                n_top = int(site_spec.split('(')[1].split(')')[0])
                atnums.extend(at.num for at in sorted_atoms[:n_top])
                continue
            raise ValueError(site_spec)
        return atnums

    def _ensure_sitedef_makes_sense(self, poscar_el, site_label, atom_nrs):
        """Raise ValueError if the atom_nrs selected for a site are wrong."""
        invalid_atom_nums = []
        invalid_elements = []
        for num in atom_nrs:
            try:
                atom = self.slab.atlist.get(num)
            except KeyError:
                invalid_atom_nums.append(num)
                continue
            if atom.el != poscar_el:
                invalid_elements.append(atom)
        if not invalid_atom_nums and not invalid_elements:
            return
        error = f'Invalid selection for {site_label}:'
        if invalid_atom_nums:
            error += f'\nAtom number(s) {invalid_atom_nums} not found in slab'
        if invalid_elements:
            error += (f'\nAtom(s) {[str(at) for at in invalid_elements]} '
                      f'do not have POSCAR element {poscar_el}')
        raise ValueError(error)

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
        unrecognized = self._interpret_symmetry_bulk_glides(unrecognized)
        unrecognized = self._interpret_symmetry_bulk_screws(unrecognized)
        unrecognized = self._interpret_symmetry_bulk_group(unrecognized)

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

    def _interpret_symmetry_bulk_glides(self, unrecognized):
        """Look up `unrecognized` for glide-plane-direction specifications."""
        _mirror_re = re.compile(r'(\s+m\[\s*(-?[012])\s*,?\s*(-?[012])\])',
                                re.IGNORECASE)
        for token, *directions in _mirror_re.findall(unrecognized):
            # All matches of mirrors are acceptable
            unrecognized = unrecognized.replace(token, '')
            direction = tuple(int(v) for v in directions)
            if (direction[0] < 0
                    or not direction[0] and direction[1] < 0):
                direction = -direction[0], -direction[1]
            self.rpars.SYMMETRY_BULK['mirror'].add(direction)
        return unrecognized

    def _interpret_symmetry_bulk_group(self, unrecognized):
        """Look up `unrecognized` for a symmetry group specification."""
        param = 'SYMMETRY_BULK'
        _group_re = re.compile(                                                 # TODO: For now borrowed from guilib. Eventually will try to instantiate a PlaneGroup
            r'(\s*(\w+)\s*(?:\[\s*-?[012]\s*-?[012]\s*\])?)',
            re.IGNORECASE
            )
        for token, group in _group_re.findall(unrecognized):
            if group not in self.grouplist:
                continue
            unrecognized = unrecognized.replace(token, '')
            if 'group' in self.rpars.SYMMETRY_BULK:
                message = 'Only one symmetry group can be given'
                raise ParameterValueError(param, message=message)
            self.rpars.SYMMETRY_BULK['group'] = token.strip().lower()

        if 'group' not in self.rpars.SYMMETRY_BULK:
            message = 'Need to specify exactly one symmetry group'
            raise ParameterValueError(param, message=message)
        return unrecognized

    def _interpret_symmetry_bulk_screws(self, unrecognized):
        """Look up `unrecognized` for screw-axis-order specifications."""
        _rotation_re = re.compile(r'(\s+r([2346]))', re.IGNORECASE)
        for token, order in _rotation_re.findall(unrecognized):
            unrecognized = unrecognized.replace(token, '')
            self.rpars.SYMMETRY_BULK['rotation'].add(int(order))
        return unrecognized

    def interpret_symmetry_cell_transform(self, assignment):
        """Assign parameter SYMMETRY_CELL_TRANSFORM (Wood's or matrix)."""
        param = 'SYMMETRY_CELL_TRANSFORM'
        self._interpret_wood_or_matrix(param, assignment)

    def interpret_symmetry_eps(self, assignment):
        """Assign parameters SYMMETRY_EPS and SYMMETRY_EPS.z."""
        param = 'SYMMETRY_EPS'
        self._ensure_no_flags_assignment(assignment)
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
        eps_value = self.interpret_numerical_parameter(assignment,
                                                       bounds=bounds,
                                                       return_only=True)
        if eps_value > 1.0:
            # pylint: disable-next=logging-format-interpolation
            _LOGGER.warning(warning_str.format(''))
        # interpret possible second value as SYMMETRY_EPS.z
        if not assignment.other_values:
            self.rpars.SYMMETRY_EPS = SymmetryEps(eps_value)
            return

        z_assignment = Assignment(assignment.other_values, param)
        eps_z_value = self.interpret_numerical_parameter(z_assignment,
                                                         bounds=bounds,
                                                         return_only=True)
        self.rpars.SYMMETRY_EPS = SymmetryEps(eps_value, eps_z_value)
        if self.rpars.SYMMETRY_EPS.z > 1.0:
            # pylint: disable-next=logging-format-interpolation
            _LOGGER.warning(warning_str.format('for z '))

    def interpret_symmetry_fix(self, assignment):                               # TODO: use symmetry groups from elsewhere once symmetry and guilib are merged
        param = 'SYMMETRY_FIX'
        group = assignment.values_str.lower()
        if group.startswith('t'):  # determine symmetry automatically
            self.rpars.reset_default(param)
            return
        if group.startswith('f'):
            self.rpars.SYMMETRY_FIX = 'p1'
            return
        if group in self.grouplist and group in ('cm', 'pmg'):
            message = f'For group {group} direction needs to be specified'
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
            self.rpars.reset_default(param)
            return
        if not assignment.other_values:
            energies = (assignment.value, assignment.value, '1')

        # (2) Three values. Any can be an "_" meaning "default".
        # Internally, we store those as an rparams.NO_VALUE. The
        # others should be positive floats.
        theo_energies = self._parse_energy_range(assignment, energies)
        if len(theo_energies) != 3:
            raise ParameterNumberOfInputsError(param)

        try:
            self.rpars.THEO_ENERGIES = TheoEnergies(*theo_energies)
        except ValueError as exc:  # No TypeError, as all are floats
            self.rpars.setHaltingLevel(1)
            raise ParameterValueError(param, message=str(exc)) from None

        # Notify in case start was adjusted to fit step
        if not self.rpars.THEO_ENERGIES.is_equivalent(theo_energies):
            _LOGGER.info(
                f'{param} parameter: (Eto - Efrom) % Estep != 0, Efrom '
                f'was corrected to {self.rpars.THEO_ENERGIES.start}'
                )

    def interpret_v0_real(self, assignment):
        """Assign parameter V0_REAL."""
        param = 'V0_REAL'
        v0r_type = assignment.value.lower()
        if v0r_type == 'rundgren':
            rundgren_constants = assignment.other_values
            if len(rundgren_constants) != 4:
                message = ('Rundgren-type function expects four '
                           'constants separated by whitespace')
                raise ParameterNumberOfInputsError(param, message=message)
            try:
                self.rpars.V0_REAL = [float(c) for c in rundgren_constants]
            except ValueError as exc:
                message = (f'Could not parse constants {rundgren_constants} '
                           'for Rundgren-type function')
                self.rpars.setHaltingLevel(1)
                raise ParameterFloatConversionError(param,
                                                    message=message) from exc
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
            self._ensure_no_flags_assignment(assignment)
            wood = self._read_woods_notation(param, assignment)
            setattr(self.rpars, param, wood)
        else:
            matrix = self._read_matrix_notation(param, assignment.values)
            setattr(self.rpars, param, matrix)

    # ------------------------ Helper methods -------------------------
    def _ensure_chemical_elements(self, param, elements):
        """Raise unless all entries are valid chemical elements."""
        if not elements:
            raise ParameterHasNoValueError(param)
        invalid = []
        for element in elements:
            try:
                _ = periodic_table.get_atomic_number(element)
            except ValueError:
                invalid.append(element)

        if invalid:
            message = f'Element(s) {invalid} not found in periodic table'
            self.rpars.setHaltingLevel(2)
            raise ParameterValueError(param, message=message)

    def _ensure_has_non_empty_slab(self, param):
        """Complain if self.slab is inappropriate for SITE_DEF."""
        if self.slab is None:
            self.rpars.setHaltingLevel(3)
            raise ParameterNeedsSlabError(param)
        if not self.slab.atlist:
            self.rpars.setHaltingLevel(3)
            raise ParameterError(param, 'Slab contains no atoms')

    def _ensure_valid_slab_element(self, param, element):
        """Raise unless element is one of the elements of self.slab."""
        self._ensure_has_non_empty_slab(param)
        known_elements = self.slab.elements
        if not known_elements:
            self.rpars.setHaltingLevel(3)
            raise ParameterError(param, 'Slab has no elements')
        if element not in known_elements:
            msg = (f'{element!r} is not a valid POSCAR element. '
                   'Known elements: ' + ', '.join(known_elements))
            self.rpars.setHaltingLevel(3)
            raise ParameterUnknownFlagError(param, message=msg)

    def _get_flag_value_from_pair(self, param, flag_value_pair):
        """Return lowercase flag and value from a pair, or complain."""
        try:
            flag, value = (s.strip()
                           for s in flag_value_pair.lower().strip().split())
        except ValueError:
            err_ = f'Expected "flag value" pairs, found "{flag_value_pair}"'
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(param, message=err_) from None
        return flag, value

    def _get_valid_slab_element_from_flag(self, param, assignment):
        """Return a valid POSCAR element from assignment; raise otherwise."""
        element = assignment.flag.capitalize()
        self._ensure_valid_slab_element(param, element)
        return element

    def _parse_energy_range(self, assignment, energies):
        """Return a list of floats for energies.

        Parameters
        ----------
        assignment : Assignment
            The assignment for the parameter. Used only for error
            reporting.
        energies : Sequence
            The energies that should be parsed. Elements are strings.
            They may also contain single underscore characters, which
            will be replaced with an Rparams.NO_VALUE.

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
        """
        try:
            return EnergyRange.parse_string_sequence(energies)
        except (TypeError, ValueError) as exc:
            self.rpars.setHaltingLevel(1)
            new_exc = (ParameterParseError if isinstance(exc, ValueError)
                       else ParameterFloatConversionError)
            raise new_exc(assignment.parameter,
                          message=assignment.values_str) from exc

    def _parse_incidence_angles(self, param, assignment):                       # TODO: Issue #147
        """Return a dictionary with incidence angles from an assignment."""
        bounds = {'THETA': NumericBounds(type_=float, range_=(-90, 90),
                                         out_of_range_event='fail'),
                  'PHI': NumericBounds(type_=float, range_=(0, 360),
                                       out_of_range_event='modulo')}
        if ',' not in assignment.values_str and len(assignment.values) != 2:
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(
                parameter=param,
                found_and_expected=(len(assignment.values), 2)
                )
        right_side = assignment.values_str.upper()
        if ',' not in assignment.values_str:
            # pylint: disable=consider-using-f-string
            # disable: Better than repeating 'assignment.values' twice
            right_side = 'THETA {}, PHI {}'.format(*assignment.values)
        self._check_n_incidence_angles_ok(param, right_side)

        angles_specs = (spec.strip().split()
                        for spec in right_side.strip().split(','))
        angles = {}
        for name, *values in angles_specs:
            if name not in ('THETA', 'PHI'):
                self.rpars.setHaltingLevel(1)
                raise ParameterUnknownFlagError(param, name)
            if len(values) != 1:
                self.rpars.setHaltingLevel(1)
                raise ParameterNumberOfInputsError(
                    parameter=param,
                    message=f'Found {len(values)} values for angle {name}'
                    )
            angles[name] = self.interpret_numerical_parameter(
                Assignment(values[0], param),
                param=f'{param} {name}',
                bounds=bounds[name],
                return_only=True
                )

        # Check for negative theta and adjust phi
        if angles['THETA'] < 0:
            angles['THETA'] = abs(angles['THETA'])
            angles['PHI'] = (angles['PHI'] + 180) % 360
        return angles

    def _check_n_incidence_angles_ok(self, param, right_side):
        """Complain if right_side contains the wrong number of angles."""
        nr_angle_specs = len(right_side.split(','))
        if nr_angle_specs != 2:
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(
                parameter=param, found_and_expected=(nr_angle_specs, 2)
                )
        if (right_side.count('THETA'), right_side.count('PHI')) != (1, 1):
            self.rpars.setHaltingLevel(1)
            raise ParameterNumberOfInputsError(
                parameter=param,
                message='Expected exactly one entry for THETA and one for PHI'
                )

    def _read_woods_notation(self, param, assignment):
        """Return a Woods notation from an Assignment, if slab exists."""
        if self.slab is None:
            message = (f'{param} parameter appears to be in Wood '
                       'notation but no slab was passed. Cannot '
                       'calculate bulk unit cell!')
            self.rpars.setHaltingLevel(2)
            raise ParameterNeedsSlabError(param, message)
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
