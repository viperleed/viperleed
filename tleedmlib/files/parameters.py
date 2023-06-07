# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 16:56:39 2020

@author: Florian Kraushofer, Michele Riva, Alexander Imre

Functions for reading from and writing to the PARAMETERS file
"""

import logging
from pathlib import Path
import re
import shutil

import numpy as np

from viperleed.tleedmlib import leedbase
from viperleed.tleedmlib.base import (strip_comments, splitSublists,
                                      readVector, readIntRange,
                                      recombineListElements)
from viperleed.tleedmlib.classes import rparams
from viperleed.tleedmlib.sections._sections import TLEEDMSection as Section


logger = logging.getLogger("tleedm.files.parameters")

# TODO: fill dict of parameter limits here (e.g. LMAX etc.)
# TODO: change module level globals to ALL_CAPS everywhere

# list of allowed parameters
_KNOWN_PARAMS = [                                                               # TODO: IntEnum?
    'ATTENUATION_EPS', 'AVERAGE_BEAMS', 'BEAM_INCIDENCE', 'BULKDOUBLING_EPS',
    'BULKDOUBLING_MAX', 'BULK_LIKE_BELOW', 'BULK_REPEAT', 'DOMAIN',
    'DOMAIN_STEP', 'ELEMENT_MIX',
    'ELEMENT_RENAME', 'FILAMENT_WF', 'FORTRAN_COMP', 'HALTING', 'INTPOL_DEG',
    'IV_SHIFT_RANGE', 'LAYER_CUTS', 'LAYER_STACK_VERTICAL', 'LMAX',
    'LOG_DEBUG', 'LOG_SEARCH', 'N_BULK_LAYERS', 'N_CORES', 'OPTIMIZE',
    'PARABOLA_FIT', 'PHASESHIFT_EPS', 'PHASESHIFTS_CALC_OLD',
    'PHASESHIFTS_OUT_OLD', 'PLOT_IV', 'RUN', 'R_FACTOR_LEGACY', 'R_FACTOR_SMOOTH',
    'R_FACTOR_TYPE', 'S_OVL', 'SCREEN_APERTURE', 'SEARCH_BEAMS', 'SEARCH_CONVERGENCE',
    'SEARCH_CULL', 'SEARCH_MAX_GEN', 'SEARCH_POPULATION', 'SEARCH_START',
    'SITE_DEF', 'SUPERLATTICE', 'SUPPRESS_EXECUTION', 'SYMMETRIZE_INPUT',
    'SYMMETRY_BULK',
    'SYMMETRY_CELL_TRANSFORM', 'SYMMETRY_EPS', 'SYMMETRY_FIND_ORI',
    'SYMMETRY_FIX', 'TENSOR_INDEX', 'TENSOR_OUTPUT', 'THEO_ENERGIES',
    'TL_VERSION', 'TL_IGNORE_CHECKSUM',
    'T_DEBYE', 'T_EXPERIMENT', 'V0_IMAG', 'V0_REAL',
    'V0_Z_ONSET', 'VIBR_AMP_SCALE', 'ZIP_COMPRESSION_LEVEL',
    ]

# _PARAM_ALIAS keys should be all lowercase, with no underscores
_PARAM_ALIAS = {
    'bulklike': 'BULK_LIKE_BELOW',
    'bulksymmetry': 'SYMMETRY_BULK',
    'compiler': 'FORTRAN_COMP',
    'fortrancompile': 'FORTRAN_COMP',
    'fortrancompiler': 'FORTRAN_COMP',
    'fdoptimize': 'OPTIMIZE',
    'fdoptimization': 'OPTIMIZE',
    'plotrfactor': 'PLOT_IV',
    'plotrfactors': 'PLOT_IV',
    'ignorechecksum': 'TL_IGNORE_CHECKSUM',
    'ivplot': 'PLOT_IV',
    'overlap': 'S_OVL',
    'mtoverlap': 'S_OVL',
    'compression_level': 'ZIP_COMPRESSION_LEVEL',
    'compression': 'ZIP_COMPRESSION_LEVEL',
    }


PARAM_LIMITS = {
    'LMAX': (1, 18),
    }

for p in _KNOWN_PARAMS:
    _PARAM_ALIAS[p.lower().replace("_", "")] = p


class ParameterError(Exception):
    '''Base class for errors raised during PARAMETERS interpretation'''

    def __init__(self, parameter, message):
        message = f"PARAMETERS file: parameter {str(parameter)}:\n"
        if message:
            message += message + " "  # add space before "input will be ignored"
        super().__init__(message)

class ParameterNotRecognizedError(ParameterError):
    '''Raised when a parameter is not recognized'''

    def __init__(self, parameter):
        super().__init__(parameter, "Parameter not recognized.")

class ParameterUnexpectedInputError(ParameterError):
    '''Raised when unexpected input is encountered'''

    def __init__(self, parameter):
        super().__init__(parameter,
                         "Encountered unexpected input."
                         "Check parameter syntax.")

class ParameterBooleanConversionError(ParameterError):
    '''Raised when a boolean conversion fails'''

    def __init__(self, parameter):
        super().__init__(parameter,
                         "Failed to convert input to boolean."
                         "Check parameter syntax.")

class ParameterFloatConversionError(ParameterError):
    '''Raised when a float conversion fails'''

    def __init__(self, parameter):
        super().__init__(parameter,
                         "Failed to convert input to float(s)."
                         "Check parameter syntax.")

class ParameterIntConversionError(ParameterError):
    '''Raised when an int conversion fails'''

    def __init__(self, parameter):
        super().__init__(parameter,
                         "Failed to convert input to integer(s)."
                         "Check parameter syntax.")

class ParameterParseError(ParameterError):
    '''Raised when parsing fails'''

    def __init__(self, parameter, supp_message=None):
        if supp_message:
            super().__init__(parameter,
                             f"Could not parse input. {supp_message}")
        else:
            super().__init__(parameter,
                            "Could not parse input."
                            "Check parameter syntax.")

class ParameterNumberOfInputsError(ParameterError):
    '''Raised when the number of inputs is unexpected'''

    def __init__(self, parameter, found_and_expected=None):
        if found_and_expected:
            super().__init__(parameter,
                             f"Expected {found_and_expected[1]} inputs, "
                             "but found {found_and_expected[0]}.")
        else:
            super().__init__(parameter,
                            "Unexpected number of inputs.")

class ParameterRangeError(ParameterError):
    '''Raised when the value is not in the allowed range'''

    def __init__(self, parameter, given_value, allowed_range):
        message = (f"Value {given_value} is outside allowed range "
                   f"({allowed_range[0]} <= {parameter} <= "
                   f"{allowed_range[1]}).")
        super().__init__(parameter, message)

class ParameterUnknownFlagError(ParameterError):
    '''Raised when an unknown flag is encountered'''

    def __init__(self, parameter, flag):
        super().__init__(parameter,
                         f"Unknown flag '{flag}' encountered.")

class ParameterNeedsFlagError(ParameterError):
    '''Raised when a flag is needed but not given'''

    def __init__(self, parameter):
        super().__init__(parameter,
                         "Parameter requires a flag.")


def _interpret_SEARCH_CONVERGENCE(rpars, flags, values,
                                  search_convergence_known=True,
                                  is_updating=False):
    """Interpret SEARCH_CONVERGENCE. Raise ParametersSyntaxError."""
    param = 'SEARCH_CONVERGENCE'
    if (not flags and len(values) == 1 and values[0].lower() == 'off'
            and not is_updating):                                                # TODO: this is the behaviour of updatePARAMETERS. Was skipping this intended there?
        rpars.GAUSSIAN_WIDTH_SCALING = 1.
        return

    if not flags:
        raise ParameterNeedsFlagError(param)
    elif flags[0].lower() not in ['dgen', 'gaussian']:
        raise ParameterUnknownFlagError(param, f"{flags[0]!r}")

    numeric = [None, None]
    for i, value in enumerate(values[:2]):
        try:
            numeric[i] = float(value)
        except ValueError:
            raise ParameterFloatConversionError(param)

    _errors = []
    if flags[0].lower() == 'gaussian':
        gauss_width, scaling = numeric
        should_update = (not is_updating
                         or gauss_width != rpars.searchConvInit["gaussian"])
        if gauss_width is not None and gauss_width > 0 and should_update:
            rpars.GAUSSIAN_WIDTH = gauss_width
            if is_updating:
                rpars.searchConvInit["gaussian"] = gauss_width
        elif should_update:
            message = "gaussian width should be a positive number "
            _errors.append(
                ParameterCustomError(param, message=message)
                )
        if scaling is not None and 0 < scaling <= 1:
            rpars.GAUSSIAN_WIDTH_SCALING = scaling
        elif scaling is not None:
            _errors.append(
                ParameterRangeError(param, scaling, (0, 1))
                )
    elif flags[0].lower() == 'dgen':
        if len(flags) == 1:
            target = 'dec'
        elif flags[1].lower() in ['dec', 'best', 'all']:
            target = flags[1].lower()
        else:
            raise ParameterUnknownFlagError(param, f"{flags[1]!r}")
        max_dgen, scaling = numeric
        should_update = (not is_updating
                         or max_dgen != rpars.searchConvInit["dgen"][target])
        if max_dgen is not None and max_dgen > 0 and should_update:
            if not search_convergence_known:  # clear default values
                rpars.SEARCH_MAX_DGEN = {"all": 0, "best": 0, "dec": 0}
                search_convergence_known = True
            rpars.SEARCH_MAX_DGEN[target] = max_dgen
            if is_updating:
                rpars.searchConvInit["dgen"][target] = max_dgen
        elif should_update:
            message = "dgen should be a positive number "
            _errors.append(
                ParameterError(param, message=message)
                )
        if scaling is not None and scaling >= 1:
            rpars.SEARCH_MAX_DGEN_SCALING[target] = scaling
        elif scaling:
            message = "scaling value cannot be smaller than 1."
            _errors.append(
                ParameterError(param, message=message)
                )

    if _errors:
        raise _errors[0]
    return search_convergence_known


def updatePARAMETERS(rp, filename='PARAMETERS', update_from=""):
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
        with open(update_from / filename, 'r') as rf:
            lines = rf.readlines()
    except FileNotFoundError as err:
        logger.error("updatePARAMETERS routine: PARAMETERS file not found.")
        raise err

    for line in lines:
        line = strip_comments(line)
        for param in ["SEARCH_KILL", "STOP"]:  # SEARCH_KILL is legacy name
            if line.upper().startswith(param):
                if not re.match(fr"\s*{param}\s*=\s*[Ff](alse)?", line):
                    rp.STOP = True
        if "=" not in line:                                                     # TODO: perhaps we should rather look at the first entry of a .split() and see if, by chance, it is a valid parameter. Then warn and continue reading.
            continue  # ignore all lines that don't have an "=" sign at all
        param, value = line.split('=', maxsplit=1)
        if param:
            # get rid of spaces and check the leftmost entry.
            param, *flags = param.split()
        param_alias = param.lower().replace("_", "")
        if param not in _KNOWN_PARAMS and param_alias in _PARAM_ALIAS:
            param = _PARAM_ALIAS[param_alias]
        if param not in _KNOWN_PARAMS:
            continue
        values = value.rstrip().split()
        if not values:
            continue                                                            # TODO: shouldn't we complain?
        if param == 'SEARCH_CONVERGENCE':
            # Previously we would pass if a line was not interpretable, but
            # new behaviour is to raise an error.
            _interpret_SEARCH_CONVERGENCE(rp, flags, values, is_updating=True)


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
        rf = filename.open('r')
    except FileNotFoundError as err:
        logger.error("PARAMETERS file not found.")
        raise err
    # read PARAMETERS:
    rpars = rparams.Rparams()
    for line in rf:
        line = strip_comments(line)
        for param in ["STOP", "SEARCH_KILL"]:
            if (line.upper().startswith(param)
                    and not re.match(fr"\s*{param}\s*=\s*[Ff](alse)?", line)):
                logger.warning(
                    f'PARAMETERS file: {param} was set at start of '
                    f'program. Modifying PARAMETERS to disable {param}; '
                    're-insert it if you actually want to stop.'
                    )
                modifyPARAMETERS(rpars, param, comment='Disabled at program '
                                 'start', path=filename.parent,
                                 suppress_ori=True)
        if "=" not in line:
            continue     # ignore all lines that don't have an "=" sign at all
        param, value = line.split('=', maxsplit=1)  # parameter at left of "="
        if not param:
            continue
        # get rid of spaces and check the leftmost entry.
        param, *flags = param.split()
        if (param not in _KNOWN_PARAMS and
                param.lower().replace("_", "") in _PARAM_ALIAS):
            param = _PARAM_ALIAS[param.lower().replace("_", "")]
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


def interpretPARAMETERS(rpars, slab=None, silent=False):                        # TODO: replace the implementation of this one with a PARAMETERSInterpreter class; make the defs below into methods, make a method that interprets one of the loop iterations, split all the branches into methods. The _interpret_SEARCH_CONVERGENCE above would become one such method.
    """
    Interprets the string values in an Rparams object, read previously by
    readPARAMETERS, to fill the parameter variables.

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
    ValueError
        Raised if the value for some parameter cannot be used.

    Returns
    -------
    None

    """

    def setBoolParameter(rp, param, value, varname=None,
                         addAllowedValues={False: [], True: []}):
        """
        Generic function for setting a given parameter to a boolean value.

        Parameters
        ----------
        rp : Rparams
            The Rparams object for which the parameter should be set
        param : str
            The name of the parameter in PARAMETERS
        value : str
            The value as found in the PARAMETERS file
        varname : str, optional
            The variable name in the Rparams class, if it differs from 'param'
        addAllowedValues : dict, optional
            Additional string values which should be interpreted as False or
            True. By default, 'f' and 'false' are False, 't' and 'true' are
            True.

        Returns
        ----------
        int
            0 if value was set, 1 otherwise

        """
        s = value.lower()
        allowedValues = {False: ['false', 'f'], True: ['true', 't']}
        for key in allowedValues:
            allowedValues[key].extend(addAllowedValues[key])
            if s in allowedValues[key]:
                v = key
                break
        else:
            raise ParameterBooleanConversionError(parameter=param)
        if not varname:
            varname = param
        setattr(rp, varname, v)
        return 0

    def getNumericalParameter(param, value,
                              type_=float,
                              range_=(None, None),
                              range_exclude=(False, False),
                              outOfRangeEvent=('fail', 'fail')):
        """
        Generic function for interpreting a given parameter as an integer or
        float value.

        Parameters
        ----------
        param : str
            The name of the parameter in PARAMETERS, for error messages only
        value : str
            The value as found in the PARAMETERS file
        type_ : int or float, optional
            Which type the variable must take
        range_ : tuple, values can be int or float or None, optional
            Defines upper and lower limits for the value. If either is set to
            None, only the other one will be interpreted. The default is None
            for both.
        range_exclude : tuple (bool, bool), optional
            Whether the upper and lower limits themselves are allowed values
            or not. The default is False for both (range is inclusive).
        outOfRangeEvent : str 'fail' or 'set', optional
            What to do if the given value lies outside the range. For 'fail',
            the value will be ignored. For 'modulo', the value will be
            brought within the range by '%' operation. For 'set', the value
            will be set to the lowest allowed value. The default is
            ('fail', 'fail').

        Returns
        ----------
        int or float or None
            The value to set, or None if interpretation failed.

        """
        if type_ not in (int, float):
            raise ValueError('type_ must be int or float')
        try:
            try:
                v = type_(value)
            except ValueError:
                if type_ != int:
                    raise
                v = int(float(value))  # for e.g. '1e6' to int
        except ValueError:
            if type is int:
                raise ParameterIntegerConversionError(parameter=param)
            else:  # type is float
                raise ParameterFloatConversionError(parameter=param)
        outOfRange = [False, False]
        if range_[0] is not None and (v < range_[0] or
                                      (range_exclude[0] and v == range_[0])):
            outOfRange[0] = True
        if range_[1] is not None and (v > range_[1] or
                                      (range_exclude[1] and v == range_[1])):
            outOfRange[1] = True
        rangeChars = ({False: "[", True: "]"}, {False: "]", True: "["})
        if all([v is not None for v in range_]):
            outOfRangeStr = (
                f'not in range {rangeChars[0][range_exclude[0]]}{range_[0]}, '
                f'{range_[1]}{rangeChars[1][range_exclude[1]]}'
                )
        elif range_[0] is not None:
            if range_exclude[0]:
                outOfRangeStr = f'less than or equal to {range_[0]}'
            else:
                outOfRangeStr = f'less than {range_[0]}'
        else:
            if range_exclude[1]:
                outOfRangeStr = f'larger than or equal to {range_[1]}'
            else:
                outOfRangeStr = f'larger than {range_[1]}'
        if any(outOfRange[i] and outOfRangeEvent[i] == 'fail'
               for i in range(0, 2)):
            message = f"Value {v} is {outOfRangeStr}."
            raise ParameterError(param, message=message)

        for i in range(0, 2):
            if not outOfRange[i]:
                continue
            if outOfRangeEvent == 'modulo':
                if not range_[1]:
                    raise ValueError('Cannot use outOfRangeEvent modulo '
                                     'if upper limit is undefined.')
                if not range_[0]:
                    range_[0] = 0
                setTo = (((v - range_[0]) % (range_[1] - range_[0]))
                         + range_[0])
            elif range_exclude[i]:
                mult = 1
                if i == 1:
                    mult = -1
                if type_ == float:
                    setTo = range_[i] + mult*1e-4
                else:
                    setTo = range_[i] + mult
            else:
                setTo = range_[i]
            # TODO: raise?
            logger.warning(f'PARAMETERS file: {param}: Value {v} is '
                           f'{outOfRangeStr}. Value will be set to {setTo}.')
            v = setTo
            break
        return v

    def setNumericalParameter(rp, param, value, varname=None, type_=float,
                              range_=(None, None), haltingOnFail=1,
                              range_exclude=(False, False),
                              outOfRangeEvent=('fail', 'fail')):
        """
        Generic function for setting a given parameter to an integer or float
        value. Interpretation is done via getNumericalParameter.

        Parameters
        ----------
        rp : Rparams
            The Rparams object for which the parameter should be set
        param : str
            The name of the parameter in PARAMETERS
        value : str
            The value as found in the PARAMETERS file
        varname : str, optional
            The variable name in the Rparams class, if it differs from 'param'
        type_ : int or float, optional
            Which type the variable must take
        range_ : tuple, values can be int or float or None, optional
            Defines upper and lower limits for the value. If either is set to
            None, only the other one will be interpreted. The default is None
            for both.
        range_exclude : tuple (bool, bool), optional
            Whether the upper and lower limits themselves are allowed values
            or not. The default is False for both (range is inclusive).
        outOfRangeEvent : str 'fail' or 'set', optional
            What to do if the given value lies outside the range. For 'fail',
            the value will be ignored. For 'modulo', the value will be
            brought within the range by '%' operation. For 'set', the value
            will be set to the lowest allowed value. The default is
            ('fail', 'fail').
        haltingOnFail : int, optional
            What halting level should be set if the 'fail' event is triggered.

        Returns
        ----------
        int
            0 if value was set, 1 otherwise

        """

        v = getNumericalParameter(param=param, value=value, type_=type_,
                                  range_=range_, haltingOnFail=haltingOnFail,
                                  range_exclude=range_exclude,
                                  outOfRangeEvent=outOfRangeEvent)
        if v is None:
            return 1
        if not varname:
            varname = param
        setattr(rp, varname, v)
        return 0

    # general definitions
    grouplist = [
        "p1", "p2", "pm", "pg", "cm", "rcm", "pmm", "pmg", "pgg", "cmm",
        "rcmm", "p4", "p4m", "p4g", "p3", "p3m1", "p31m", "p6", "p6m"]

    loglevel = logger.level
    if silent:
        logger.setLevel(logging.ERROR)
    # track some parameters while reading
    searchConvRead = False
    # define order that parameters should be read in
    orderedParams = ["LOG_DEBUG", "RUN"]
    param_names = [p for p in orderedParams if p in rpars.readParams]
    param_names.extend(p for p in _KNOWN_PARAMS
                       if p in rpars.readParams and p not in param_names)

    # Flatten out the (possibly multiple) assignments read from
    # the PARAMETERS file for each parameter
    flat_params = ((p_name, *flags_and_value)
                   for p_name in param_names
                   for flags_and_value in rpars.readParams[p_name])
    domainsIgnoreParams = [
        'BULK_REPEAT', 'ELEMENT_MIX', 'ELEMENT_RENAME',
        'LAYER_CUTS', 'N_BULK_LAYERS', 'SITE_DEF', 'SUPERLATTICE',
        'SYMMETRY_CELL_TRANSFORM', 'TENSOR_INDEX', 'TENSOR_OUTPUT']
    _is_doman_calc = Section.DOMAINS in rpars.RUN or rpars.domainParams
    for param, flags, right_side in flat_params:
        if _is_doman_calc and param in domainsIgnoreParams:
            continue
        values = right_side.split()
        value, *other_values = values                                           # TODO: probably cleaner to use a namedtuple or a simple @dataclass?
        # parameters not interpreted right now
        if param == 'VIBR_AMP_SCALE':
            rpars.VIBR_AMP_SCALE.extend(right_side.split(","))
        # simple bool parameters
        elif param in ['LOG_DEBUG', 'LOG_SEARCH', 'PHASESHIFTS_CALC_OLD',
                       'PHASESHIFTS_OUT_OLD', 'R_FACTOR_LEGACY',
                       'SUPPRESS_EXECUTION', 'SYMMETRIZE_INPUT',
                       'SYMMETRY_FIND_ORI', 'TL_IGNORE_CHECKSUM']:
            setBoolParameter(rpars, param, value)
        # slightly more complicated bools
        elif param == 'LAYER_STACK_VERTICAL':
            setBoolParameter(rpars, param, value,
                             addAllowedValues={False: 'c', True: 'z'})
        # positive-only integers
        elif param in ['BULKDOUBLING_MAX', 'N_CORES', 'SEARCH_MAX_GEN',
                       'TENSOR_INDEX']:
            setNumericalParameter(rpars, param, value, type_=int,
                                  range_=(1, None))
        # positive-only floats
        elif param in ['T_DEBYE', 'T_EXPERIMENT',
                       'V0_IMAG', 'TL_VERSION', 'S_OVL']:
            setNumericalParameter(rpars, param, value, range_=(0, None))
        # simple numericals
        elif param == 'V0_Z_ONSET':
            setNumericalParameter(rpars, param, value)
        elif param == 'ATTENUATION_EPS':
            setNumericalParameter(rpars, param, value, range_=(1e-6, 1),
                                  range_exclude=(False, True))
        elif param == 'BULKDOUBLING_EPS':
            setNumericalParameter(rpars, param, value,
                                  range_=(0.0001, None),
                                  outOfRangeEvent=('set', 'fail'))
        elif param == 'BULK_LIKE_BELOW':
            setNumericalParameter(rpars, param, value, range_=(0, 1),
                                  range_exclude=(True, True))
        elif param == 'HALTING':
            setNumericalParameter(rpars, param, value, type_=int,
                                  range_=(1, 3))
        elif param == 'N_BULK_LAYERS':
            setNumericalParameter(rpars, param, value, type_=int,
                                  range_=(1, 2), haltingOnFail=2)
        elif param == 'R_FACTOR_SMOOTH':
            setNumericalParameter(rpars, param, value, type_=int,
                                  range_=(0, 999))
        elif param == 'R_FACTOR_TYPE':
            setNumericalParameter(rpars, param, value, type_=int,
                                  range_=(1, 2))
        elif param == 'SCREEN_APERTURE':
            setNumericalParameter(rpars, param, value, range_=(0, 180))
        elif param == 'SEARCH_POPULATION':
            r = setNumericalParameter(rpars, param, value, type_=int,
                                      range_=(1, None))
            if r == 0 and rpars.SEARCH_POPULATION < 16:
                logger.warning('SEARCH_POPULATION is very small. A '
                               'minimum value of 16 is recommended.')
        elif param == 'SYMMETRY_EPS':
            r = setNumericalParameter(rpars, param, value,
                                      range_=(1e-100, None),
                                      outOfRangeEvent=("set", "set"))
            if r == 0 and rpars.SYMMETRY_EPS > 1.0:
                logger.warning(
                    'PARAMETERS file: SYMMETRY_EPS: Given '
                    'value is greater than one Ångström. This is a '
                    'very loose constraint and might lead to '
                    'incorrect symmetry detection. Be sure to check '
                    'the output!')
                rpars.setHaltingLevel(1)
            if other_values:
                r = setNumericalParameter(rpars, param, other_values[0],
                                          range_=(1e-100, None),
                                          outOfRangeEvent=("set", "set"),
                                          varname='SYMMETRY_EPS_Z')
                if r == 0 and rpars.SYMMETRY_EPS_Z > 1.0:
                    logger.warning(
                        'PARAMETERS file: SYMMETRY_EPS: Given '
                        'value for z is greater than one Ångström. This is a '
                        'very loose constraint and might lead to '
                        'incorrect symmetry detection. Be sure to check '
                        'the output!')
                    rpars.setHaltingLevel(1)
                if r != 0:
                    rpars.SYMMETRY_EPS_Z = rpars.SYMMETRY_EPS
            else:
                rpars.SYMMETRY_EPS_Z = rpars.SYMMETRY_EPS
        elif param == 'ZIP_COMPRESSION_LEVEL':
            setNumericalParameter(rpars, param, value, type_=int,
                                  range_=(0, 9))
        # non-trivial parameters
        elif param in ['AVERAGE_BEAMS', 'BEAM_INCIDENCE']:  # same syntax
            if param == 'AVERAGE_BEAMS':   # special cases
                right_side = right_side.lower().strip()
                if right_side in ['off', 'none', 'false', 'f']:
                    rpars.AVERAGE_BEAMS = False
                    continue
                if right_side.lower() in ['all', 'perpendicular', 'perp']:
                    rpars.AVERAGE_BEAMS = (0., 0.)
                    continue
            range_ = {'THETA': (-90, 90), 'PHI': (0, 360)}
            outOfRangeEvent = {'THETA': ('fail', 'fail'),
                               'PHI': ('modulo', 'modulo')}
            d = {'THETA': 0, 'PHI': 0}
            if ',' in right_side:
                sublists = splitSublists(values, ',')
                for sl in sublists:
                    for name in ['THETA', 'PHI']:
                        if sl[0].upper() == name:
                            d[name] = getNumericalParameter(
                                f'{param} {name}', sl[1],
                                range_=range_[name],
                                outOfRangeEvent=outOfRangeEvent[name]
                                )
                            break
                    else:
                        raise ParameterUnknownFlagError(parameter=param,
                                                        flag=sl[0])
            else:
                if len(values) != 2:
                    raise ParameterNumberOfInputsError(
                        parameter=param,
                        found_and_expected=(len(values), 2)
                        )
                    rpars.setHaltingLevel(1)
                    continue
                for ind, name in enumerate(['THETA', 'PHI']):
                    d[name] = getNumericalParameter(
                        f'{param} {name}', values[ind],
                        range_=range_[name],
                        outOfRangeEvent=outOfRangeEvent[name]
                        )
            if any(v is None for v in d.values()):
                raise ParameterParseError(parameter=param)
            # check for negative theta and adjust phi
            if d['THETA'] < 0:
                d['THETA'] = abs(d['THETA'])
                d['PHI'] = (d['PHI'] + 180) % 360
            # finally set
            if param == 'AVERAGE_BEAMS':
                rpars.AVERAGE_BEAMS = (d['THETA'], d['PHI'])
            else:
                rpars.THETA = d['THETA']
                rpars.PHI = d['PHI']
        elif param == 'BULK_REPEAT':
            s = right_side.lower()
            if "[" not in s:
                if "(" not in s:
                    try:
                        rpars.BULK_REPEAT = abs(float(value))
                    except ValueError:
                        raise ParameterFloatConversionError(parameter=param)
                else:
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
                                rpars.BULK_REPEAT = v
                            else:  # c
                                rpars.BULK_REPEAT = slab.ucell[2, 2] * v        # TODO: what if slab is None??
            else:  # vector
                vec = readVector(s, slab.ucell)                                 # TODO: what if slab is None??
                if vec is None:
                    raise ParameterParseError(parameter=param)
                rpars.BULK_REPEAT = vec
        elif param == 'DOMAIN':
            # check name
            name = flags[0] if flags else ""
            names = [n for (n, _) in rpars.DOMAINS]
            if name in names:
                # TODO: raise in this case?
                error_message = ("Multiple sources defined for"
                                 f"domain {name}.")
                raise ParameterError(parameter=param,
                                           message=error_message)
            if not name:  # get unique name
                i = 1
                while str(i) in names:
                    i += 1
                name = str(i)
            # check path
            right_side = right_side.strip()
            if Path(right_side).exists():
                path = right_side
            elif Path(right_side).with_suffix(".zip").is_file():
                path = right_side + ".zip"
            else:
                error_message = (f"Value for DOMAIN {name} could not be "
                                 "interpreted as either a path or a .zip file.")
                raise ParameterError(parameter=param,
                                           message=error_message)
            rpars.DOMAINS.append((name, path))
        elif param == 'DOMAIN_STEP':
            try:
                i = int(value)
            except ValueError:
                raise ParameterIntConversionError(parameter=param)
            if not (1 <= i <= 100):
                raise ParameterRangeError(parameter=param,
                                          given_value=i,
                                          allowed_range=(1, 100))
            if 100 % i != 0:
                message = (f"100 is not divisible by given value {i}.")
                raise ParameterError(parameter=param,message=message)
            else:
                rpars.DOMAIN_STEP = i
        elif param == 'ELEMENT_MIX':                                            # TODO: don't we check to avoid conflicts for ELEMENT_MIX and ELEMENT_RENAME?
            ptl = [el.lower() for el in leedbase.PERIODIC_TABLE]                # TODO: nicer to use the leedbase element getter function and catch exceptions
            found = False
            for el in values:
                if el.lower() not in ptl:
                    message = f"Element {el} not found in periodic table."
                    raise ParameterError(parameter=param,
                                               message=message)
            if not found:
                rpars.ELEMENT_MIX[flags[0].capitalize()] = [el.capitalize()
                                                            for el in values]
        elif param == 'ELEMENT_RENAME':
            ptl = [el.lower() for el in tl.leedbase.PERIODIC_TABLE]             # TODO: nicer to use the leedbase element getter function and catch exceptions
            if value.lower() not in ptl:
                message = f"Element {el} not found in periodic table."
                raise ParameterError(parameter=param,
                                            message=message)
            else:
                rpars.ELEMENT_RENAME[flags[0].capitalize()] = (
                    value.capitalize()
                    )
        elif param == 'FILAMENT_WF':
            if value.lower() == 'w':
                rpars.FILAMENT_WF = 4.5
            elif value.lower() == 'lab6':
                rpars.FILAMENT_WF = 2.65
            else:
                setNumericalParameter(rpars, param, value)
        elif param == 'FORTRAN_COMP':
            if (not flags and value.lower() in ["ifort", "gfortran"]
                    and not other_values):
                rpars.getFortranComp(comp=value.lower())
            elif (flags and flags[0].lower() == "mpi"
                  and value.lower() in ["mpifort", "mpiifort"]):
                rpars.getFortranMpiComp(comp=value.lower())
            else:
                delim = value[0]     # should be quotation marks                # TODO: is this needed now that we use f-strings?
                if delim not in ["'", '"']:
                    message = ("No valid shorthand and not delimited by "
                               "quotation marks.")
                    raise ParameterError(parameter=param,
                                               message=message)
                else:
                    setTo = right_side.split(delim)[1]
                if not flags:
                    rpars.FORTRAN_COMP[0] = setTo
                elif flags[0].lower() == "post":
                    rpars.FORTRAN_COMP[1] = setTo
                elif flags[0].lower() == "mpi":
                    rpars.FORTRAN_COMP_MPI[0] = setTo
                elif flags[0].lower() == "mpipost":
                    rpars.FORTRAN_COMP_MPI[1] = setTo
                else:
                    raise ParameterUnknownFlagError(parameter=param,
                                                    flag=plist[1])
        elif param == 'INTPOL_DEG':
            if value in ('3', '5'):
                rpars.INTPOL_DEG = int(value)
            else:
                message = ("Only degree 3 and 5 interpolation supported at the "
                           "moment.")
                raise ParameterError(parameter=param,
                                             message=message)
        elif param == 'IV_SHIFT_RANGE':
            if len(values) not in (2, 3):
                raise ParameterNumberOfInputsError(parameter=param)
            try:
                fl = [float(s) for s in values]
            except ValueError:
                logger.warning('PARAMETERS file: Failed to convert '
                               'IV_SHIFT_RANGE input to floats Input will '
                               'be ignored')
                rpars.setHaltingLevel(1)
                continue
            if fl[1] < fl[0]:
                logger.warning('PARAMETERS file: IV_SHIFT_RANGE end energy '
                               'has to be greater than or equal to start '
                               'energy. Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue

            for i in range(0, 2):
                rpars.IV_SHIFT_RANGE[i] = fl[i]
            if len(fl) == 3 and fl[2] <= 0:
                logger.warning('PARAMETERS file: IV_SHIFT_RANGE step has '
                               'to be positive. Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            rpars.IV_SHIFT_RANGE[2] = fl[2]
        elif param == 'LAYER_CUTS':
            # some simple filtering here, but leave as list of strings
            if all(c in right_side for c in "<>"):
                logger.warning('PARAMETERS file: LAYER_CUTS parameter: '
                               'Cannot parse list containing both "<" and ">" '
                               'signs. Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            elif any(c in right_side for c in "<>"):
                newlist = []
                for s in values:
                    s = s.replace("<", " < ")
                    s = s.replace(">", " > ")
                    newlist.extend(s.split())
                values = newlist
            rgx = re.compile(r'\s*(dz|dc)\s*\(\s*(?P<cutoff>[0-9.]+)\s*\)')
            for (i, s) in enumerate(values):
                if "dz" in s.lower() or "dc" in s.lower():
                    m = rgx.match(right_side.lower())
                    if m:
                        try:
                            float(m.group('cutoff'))
                            values[i] = m.group(0)
                        except Exception:                                       # TODO: catch better; only 1 statement in try.
                            logger.warning(
                                'PARAMETERS file: LAYER_CUTS parameter: Could '
                                'not parse function ' + s + '. Input will be '
                                'ignored.')
                            rpars.setHaltingLevel(1)
                            continue
                elif not (s == "<" or s == ">"):
                    try:
                        float(s)
                    except Exception:
                        raise ParameterParseError(parameter=param)
            rpars.LAYER_CUTS = values
        elif param == 'LMAX':
            _min, _max = rpars.get_limits(param)

            values = re.sub(r'[:-]', ' ', right_side).split()
            try:
                il = [int(v) for v in values]
            except ValueError:
                logger.warning('PARAMETERS file: LMAX parameter: Could not '
                               f'parse {right_side!r} as integer(s). Input'
                               'will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            if len(il) > 2:
                logger.warning(
                    'PARAMETERS file: LMAX parameter: Expected one or two '
                    f'values, found {len(il)}. First two values will be used.'
                    )
                il = il[:2]
            if len(il) == 1:
                if not _min < il[0] <= _max:
                    _, il[0], _ = sorted((_min, il[0], _max))
                    logger.warning(
                        f'PARAMETERS file: LMAX must be between {_min} '
                        f'and {_max}. Value will be set to {il[0]}.'
                        )
                rpars.LMAX = [il[0], il[0]]
            elif len(il) == 2:
                if il[1] < il[0]:
                    il.reverse()
                if il[0] < _min:
                    logger.warning('PARAMETERS file: LMAX lower bound must be '
                                   f'positive. Value will be set to {_min}.')
                    il[0] = _min
                if il[1] > _max:
                    logger.warning('PARAMETERS file: LMAX values greater than '
                                   f'{_max} are currently not supported. Upper'
                                   f' bound will be set to {_max}.')
                    il[1] = _max
                rpars.LMAX = il
        elif param == 'OPTIMIZE':
            if not flags:
                logger.warning('PARAMETERS file: OPTIMIZE: Parameter to '
                               'optimize not defined. Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            which = flags[0].lower()
            if which not in ['theta', 'phi', 'v0i',
                             'a', 'b', 'c', 'ab', 'abc', 's_ovl']:
                logger.warning('PARAMETERS file: OPTIMIZE: Parameter '
                               f'{which!r} not recognized. Input will '
                               'be ignored.')
                rpars.setHaltingLevel(1)
                continue
            rpars.OPTIMIZE['which'] = which
            if not other_values:
                try:
                    rpars.OPTIMIZE['step'] = float(value)
                except ValueError:
                    pass   # will be caught below
                else:
                    continue
            sublists = splitSublists(values, ',')
            for sl in sublists:
                if len(sl) != 2:
                    logger.warning('PARAMETERS file: OPTIMIZE: Expected '
                                   '"flag value" pairs, found ' + " ".join(sl))
                    rpars.setHaltingLevel(1)
                    continue
                flag = sl[0].lower()
                if flag not in ['step', 'convergence',
                                'minpoints', 'maxpoints', 'maxstep']:
                    logger.warning(f'PARAMETERS file: OPTIMIZE: Flag {sl[0]!r}'
                                   ' not recognized.')
                    rpars.setHaltingLevel(1)
                    continue
                partype = {'step': float, 'convergence': float,
                           'minpoints': int, 'maxpoints': int,
                           'maxstep': float}
                value_error = ('PARAMETERS file: OPTIMIZE: Value '
                               f'{sl[1]} is not valid for flag {sl[0]}. '
                               'Value will be ignored.')
                try:
                    rpars.OPTIMIZE[flag] = partype[flag](sl[1])
                except ValueError:
                    logger.warning(value_error)
                    rpars.setHaltingLevel(1)
                    continue
        elif param == 'PARABOLA_FIT':
            if value == 'off':
                rpars.PARABOLA_FIT['type'] = 'none'
                continue
            sublists = splitSublists(values, ',')
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
                        rpars.PARABOLA_FIT['type'] = sl[1]
                    else:
                        logger.warning(value_error)
                        rpars.setHaltingLevel(1)
                elif flag in ('alpha', 'mincurv', 'localize'):
                    try:
                        f = float(sl[1])
                    except ValueError:
                        f = -1
                    if f >= 0:
                        rpars.PARABOLA_FIT[flag] = f
                    else:
                        logger.warning(value_error)
                        rpars.setHaltingLevel(1)
        elif param == 'PHASESHIFT_EPS':
            try:
                f = float(value)
            except ValueError:
                s = value.lower()[0]
                ps_eps_default_dict = rpars.get_default(param)
                f = ps_eps_default_dict.get(s, None)
                if f is None:
                    logger.warning('PARAMETERS file: PHASESHIFT_EPS: '
                                   'Could not convert value to float. '
                                   'Input will be ignored.')
                    rpars.setHaltingLevel(1)
                    continue
            if 0 < f < 1:
                rpars.PHASESHIFT_EPS = f
            else:
                logger.warning(
                    'PARAMETERS file: PHASESHIFT_EPS: Unexpected value '
                    '(should be between 0 and 1). Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'PLOT_IV':
            if not flags:
                logger.warning('PARAMETERS file: PLOT_IV: Found no flag. '
                               'Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            flag = flags[0].lower()
            value = value.lower()
            if flag not in ('color', 'colour', 'colors', 'colours', 'perpage',
                            'border', 'borders', 'axes', 'legend', 'legends',
                            'layout', 'overbar', 'overline', 'plot'):
                logger.warning(f'PARAMETERS file: PLOT_IV: Flag {flag} not '
                               'recognized. Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            if flag == 'plot':
                # should it plot?
                if value in ('true'):
                     rpars.PLOT_IV['plot'] = True
                elif value in ('false', 'none'):
                    rpars.PLOT_IV['plot'] = False
            if flag in ('border', 'borders', 'axes'):
                if value in ('all', 'none'):
                    rpars.PLOT_IV['axes'] = value
                elif value in ('less', 'lb'):
                    rpars.PLOT_IV['axes'] = 'lb'
                elif value in ('bottom', 'b'):
                    rpars.PLOT_IV['axes'] = 'b'
                else:
                    logger.warning(f'PARAMETERS file: PLOT_IV {flag}: Value '
                                   'not recognized. Input will be ignored.')
            elif flag in ('color', 'colour', 'colors', 'colours'):
                rpars.PLOT_IV['colors'] = values
            elif flag in ('legend', 'legends'):
                if value in ('all', 'first', 'none'):
                    rpars.PLOT_IV['legend'] = value
                elif value in ('topright', 'tr'):
                    rpars.PLOT_IV['legend'] = 'tr'
                else:
                    logger.warning(f'PARAMETERS file: PLOT_IV {flag}: Value '
                                   'not recognized. Input will be ignored.')
                    continue
            elif flag in ('overbar', 'overline'):
                if value.startswith("t"):
                    rpars.PLOT_IV['overbar'] = True
                elif value.startswith("f"):
                    rpars.PLOT_IV['overbar'] = False
                else:
                    logger.warning(f'PARAMETERS file: PLOT_IV {flag}: Value '
                                   'not recognized. Input will be ignored.')
                    continue
            elif flag in ('perpage', 'layout'):
                if not other_values:
                    try:
                        i = int(value)
                    except (ValueError, IndexError):
                        logger.warning(
                            'PARAMETERS file: PLOT_IV perpage: Could not '
                            'convert value to integer. Input will be ignored.')
                        continue
                    if i <= 0:
                        logger.warning(
                            'PARAMETERS file: PLOT_IV perpage: Value has '
                            'to be positive integer. Input will be ignored.')
                        continue
                    rpars.PLOT_IV['perpage'] = i
                elif len(values) >= 2:
                    try:
                        il = [int(v) for v in values[:2]]
                    except ValueError:
                        logger.warning(
                            'PARAMETERS file: PLOT_IV perpage: Could not '
                            'convert values to integers. Input will be '
                            'ignored.')
                        continue
                    if any(i <= 0 for i in il):
                        logger.warning(
                            'PARAMETERS file: PLOT_IV perpage: Values '
                            'have to be positive integers. Input will be '
                            'ignored.')
                        continue
                    rpars.PLOT_IV['perpage'] = tuple(il)
        elif param == 'RUN':
            rl = []
            for s in values:
                try:
                    rl.extend(Section.sequence_from_string(s))
                except ValueError as err:
                    logger.warning(
                        'PARAMETERS file: RUN: Could not interpret '
                        f'value {s}, skipping value. Info: {err}'
                        )
                    rpars.setHaltingLevel(2)
            if Section.DOMAINS in rl:
                logger.info('Found domain search.')
            if rl:
                if rl[0] is not Section.INITIALIZATION:
                    rl.insert(0, Section.INITIALIZATION)
                rpars.RUN = rl
            else:
                # could in principle try carrying on with some default, but
                #   this is probably a case in which a hard stop is
                #   appropriate
                logger.warning('PARAMETERS file: RUN was defined, but '
                               'no values were read. Execution will stop.')
                rpars.setHaltingLevel(3)
        elif param == 'SEARCH_BEAMS':
            value = value.lower()
            if value.startswith(('0', 'a')):
                rpars.SEARCH_BEAMS = 0
            elif value.startswith(('1', 'i')):
                rpars.SEARCH_BEAMS = 1
            elif value.startswith(('2', 'f')):
                rpars.SEARCH_BEAMS = 2
            else:
                logger.warning('PARAMETERS file: SEARCH_BEAMS: value not '
                               'recognized. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'SEARCH_CONVERGENCE':
            try:
                searchConvRead = _interpret_SEARCH_CONVERGENCE(
                    rpars, flags, values,
                    search_convergence_known=searchConvRead
                    )
            except ParametersSyntaxError as exc:
                logger.warning(exc)
                rpars.setHaltingLevel(1)
        elif param == 'SEARCH_CULL':
            try:
                f = float(value)
            except ValueError:
                logger.warning('PARAMETERS file: SEARCH_CULL: could not '
                               'interpret value.')
                rpars.setHaltingLevel(1)
                continue
            if f >= 1:
                if f - int(f) < 1e-6:
                    rpars.SEARCH_CULL = int(f)
                else:
                    logger.warning('PARAMETERS file: SEARCH_CULL: Values '
                                   'greater than 1 have to be integers.')
                    rpars.setHaltingLevel(1)
                    continue
            elif f >= 0:
                rpars.SEARCH_CULL = f
            else:
                logger.warning('PARAMETERS file: SEARCH_CULL cannot be '
                               'negative.')
                rpars.setHaltingLevel(1)
                continue
            if other_values:
                next_value = other_values[0].lower()
                if next_value in ["clone", "genetic", "random"]:
                    rpars.SEARCH_CULL_TYPE = next_value
                else:
                    logger.warning('PARAMETERS file: SEARCH_CULL type '
                                   'not recognized.')
                    rpars.setHaltingLevel(1)
        elif param == 'SEARCH_START':
            value = value.lower()
            if value.startswith("rand"):
                rpars.SEARCH_START = "random"
            elif value.startswith("center"):
                rpars.SEARCH_START = "centered"
            elif value.startswith("control"):
                rpars.SEARCH_START = "control"
            elif value.startswith("cr"):
                rpars.SEARCH_START = "crandom"
            else:
                logger.warning('PARAMETERS file: SEARCH_START: flag '
                               'not recognized.')
                rpars.setHaltingLevel(1)
        elif param == 'SITE_DEF':
            newdict = {}
            sublists = splitSublists(values, ',')
            for sl in sublists:
                atnums = []
                for i in range(1, len(sl)):
                    ir = readIntRange(sl[i])
                    if len(ir) > 0:
                        atnums.extend(ir)
                    elif "top(" in sl[i]:
                        if slab is None:
                            logger.warning(
                                'PARAMETERS file: SITE_DEF '
                                'parameter contains a top() function, '
                                'but no slab was passed. The atoms '
                                'will be assigned the default site '
                                'type instead.'
                                )
                            rpars.setHaltingLevel(1)
                            continue
                        n = int(sl[i].split('(')[1].split(')')[0])
                        csatlist = sorted(slab.atlist,
                                          key=lambda atom: atom.pos[2])
                        while n > 0:
                            at = csatlist.pop()
                            if at.el == flags[0]:
                                atnums.append(at.oriN)
                                n -= 1
                    else:
                        logger.error('PARAMETERS file: Problem with '
                                     'SITE_DEF input format')
                        raise ParametersSyntaxError("Invalid syntax",           # TODO: more explicit message?
                                                    parameter=param)
                newdict[sl[0]] = atnums
            rpars.SITE_DEF[flags[0]] = newdict
        elif param in ['SUPERLATTICE', 'SYMMETRY_CELL_TRANSFORM']:
            if 'M' not in flags:
                if slab is None:
                    logger.warning(
                        f'PARAMETERS file: {param} parameter appears to be '
                        'in Wood notation, but no slab was passed; cannot'
                        ' calculate bulk unit cell!'
                        )
                    rpars.setHaltingLevel(2)
                else:
                    setattr(rpars,
                            param,
                            leedbase.readWoodsNotation(right_side, slab.ucell))
                    if param == 'SUPERLATTICE':
                        rpars.superlattice_defined = True
            else:
                sublists = splitSublists(values, ',')
                if not len(sublists) == 2:
                    logger.warning(f'PARAMETERS file: error reading {param} '
                                   'matrix: number of lines is not equal 2.')
                    rpars.setHaltingLevel(2)
                else:
                    write = True
                    nl = []
                    for sl in sublists:
                        if len(sl) == 2:
                            try:
                                nl.append([float(s) for s in sl])
                            except ValueError:
                                logger.warning(
                                    f'PARAMETERS file: error reading {param} '
                                    f'matrix: could not convert {sl} to floats.'
                                    )
                                rpars.setHaltingLevel(2)
                                write = False
                        else:
                            logger.warning(
                                f'PARAMETERS file: error reading {param} '
                                'matrix: number of columns is not equal 2.'
                                )
                            rpars.setHaltingLevel(2)
                            write = False
                    if write:
                        setattr(rpars, param, np.array(nl, dtype=float))
                        if param == 'SUPERLATTICE':
                            rpars.superlattice_defined = True
        elif param == 'SYMMETRY_BULK':
            recombined_list = []                                                # TODO: use ast
            while values:
                v = values.pop(0).lower()
                if "[" in v and "]" in values[0]:
                    v += " " + values.pop(0).lower()
                recombined_list.append(v)
            for v in recombined_list:
                if v.startswith("r"):
                    try:
                        i = int(v[1])
                    except (ValueError, IndexError):
                        logger.warning("PARAMETERS file: error reading "
                                       f"value {v!r} in SYMMETRY_BULK.")
                        rpars.setHaltingLevel(2)
                        continue
                    if 'rotation' not in rpars.SYMMETRY_BULK:
                        rpars.SYMMETRY_BULK['rotation'] = []
                    if i not in rpars.SYMMETRY_BULK['rotation']:
                        rpars.SYMMETRY_BULK['rotation'].append(i)
                elif v.startswith("m"):
                    if "[" not in v or "]" not in v:
                        logger.warning(
                            f"PARAMETERS file: error reading value {v!r} "
                            "in SYMMETRY_BULK: no direction recognized."
                            )
                        rpars.setHaltingLevel(2)
                        continue
                    str_vals = v.split("[")[1].split("]")[0].split()
                    if len(str_vals) != 2:
                        logger.warning(
                            f"PARAMETERS file: error reading value {v!r} "
                            "in SYMMETRY_BULK: expected 2 values."
                            )
                        rpars.setHaltingLevel(2)
                        continue
                    try:
                        int_vals = tuple(int(v) for v in str_vals)
                    except (ValueError, IndexError):
                        logger.warning("PARAMETERS file: error reading "
                                       f"value {v!r} in SYMMETRY_BULK.")
                        rpars.setHaltingLevel(2)
                        continue
                    if int_vals[0] < 0:
                        int_vals = (-int_vals[0], -int_vals[1])
                    if 'mirror' not in rpars.SYMMETRY_BULK:
                        rpars.SYMMETRY_BULK['mirror'] = []
                    if int_vals in rpars.SYMMETRY_BULK['mir']:
                        rpars.SYMMETRY_BULK['mirror'].append(int_vals)
                else:
                    try:
                        rpars.SYMMETRY_BULK['group'] = v
                    except ValueError as err:
                        logger.warning("PARAMETERS file: error reading "
                                       f"value {v!r} in SYMMETRY_BULK: {err}.")
                        rpars.setHaltingLevel(2)
        elif param == 'SYMMETRY_FIX':
            group = value.lower()
            if group.startswith('t'):
                continue  # same as default, determine symmetry automatically
            if group.startswith('f'):
                rpars.SYMMETRY_FIX = 'p1'
                continue
            if group in grouplist and group in ("cm", "pmg"):
                logger.warning('PARAMETERS file: SYMMETRY_FIX: For '
                               f'group {group}, direction needs to be '
                               'specified. Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            if group in grouplist:
                rpars.SYMMETRY_FIX = group
                continue
            if group.startswith(("pm", "pg", "cm", "rcm", "pmg")):
                # regex to read
                rgx = re.compile(
                    r'\s*(?P<group>(pm|pg|cm|rcm|pmg))\s*'
                    + r'\[\s*(?P<i1>[-012]+)\s+(?P<i2>[-012]+)\s*\]')
                m = rgx.match(right_side.lower())
                if not m:
                    logger.warning(
                        'PARAMETERS file: SYMMETRY_FIX: Could '
                        'not parse given value. Input will be ignored.')
                    rpars.setHaltingLevel(1)
                    continue
                i1 = i2 = -2
                group = m.group('group')
                try:
                    i1 = int(m.group('i1'))
                    i2 = int(m.group('i2'))
                except ValueError:
                    logger.warning('PARAMETERS file: SYMMETRY_FIX: '
                                   'Could not parse given value. '
                                   'Input will be ignored.')
                    rpars.setHaltingLevel(1)
                    continue
                if (group in ["pm", "pg", "cm", "rcm", "pmg"]
                        and i1 in range(-1, 3) and i2 in range(-1, 3)):
                    rpars.SYMMETRY_FIX = m.group(0)
                else:
                    logger.warning('PARAMETERS file: SYMMETRY_FIX: '
                                   'Could not parse given value. '
                                   'Input will be ignored.')
                    rpars.setHaltingLevel(1)
            else:
                logger.warning('PARAMETERS file: SYMMETRY_FIX: Could not '
                               'parse given value. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'TENSOR_OUTPUT':
            nl = recombineListElements(values, '*')
            for s in nl:
                s = re.sub(r'[Tt](rue|RUE)?', '1', s)
                s = re.sub(r'[Ff](alse|ALSE)?', '0', s)
                try:
                    v = int(s)
                    if v not in (0, 1):
                        logger.warning(
                            'PARAMETERS file: Problem with '
                            f'TENSOR_OUTPUT input format: Found value {v}'
                            ', expected 0 or 1. Value will be ignored.'
                            )
                        rpars.setHaltingLevel(1)
                    else:
                        rpars.TENSOR_OUTPUT.append(v)
                except ValueError:
                    if '*' in s:
                        sl = s.split('*')
                        try:
                            r = int(sl[0])
                            v = int(sl[1])
                        except ValueError:
                            logger.warning(
                                'PARAMETERS file: Problem with '
                                'TENSOR_OUTPUT input format: could not read '
                                f'value {s}, value will be ignored.')
                            rpars.setHaltingLevel(1)
                        if v not in (0, 1):
                            logger.warning(
                                'PARAMETERS file: Problem with '
                                'TENSOR_OUTPUT input format: Found value '
                                f'{v}, expected 0 or 1. Value will be '
                                'ignored.')
                            rpars.setHaltingLevel(1)
                        else:
                            rpars.TENSOR_OUTPUT.extend([v]*r)
                    else:
                        logger.warning(
                            'PARAMETERS file: Problem with '
                            'TENSOR_OUTPUT input format: could not read '
                            f'value {s}, value will be ignored.')
                        rpars.setHaltingLevel(1)
        elif param == 'THEO_ENERGIES':
            if not other_values:
                # single value input - only one energy requested
                try:
                    f = float(value)
                except ValueError:
                    logger.warning(
                        'PARAMETERS file: Failed to convert '
                        'THEO_ENERGIES input to floats Input will be ignored')
                    rpars.setHaltingLevel(1)
                else:
                    if f > 0:
                        rpars.THEO_ENERGIES = [f, f, 1]
                        continue
                    logger.warning('PARAMETERS file: Unexpected input for '
                                   'THEO_ENERGIES. Input will be ignored.')
                    rpars.setHaltingLevel(1)
                continue
            if len(values) != 3:
                logger.warning('PARAMETERS file: Unexpected input for '
                               'THEO_ENERGIES. Input will be ignored.')
                rpars.setHaltingLevel(1)
                continue
            fl = []
            defined = 0
            for s in values:
                if s == "_":
                    fl.append(-1)
                    continue  # s in values loop
                try:
                    f = float(s)
                except ValueError:
                    logger.warning('PARAMETERS file: Failed to convert '
                                   'THEO_ENERGIES input to floats. '
                                   'Input will be ignored.')
                    rpars.setHaltingLevel(1)
                else:
                    if f > 0:
                        fl.append(f)
                        defined += 1
                    else:
                        logger.warning('PARAMETERS file: THEO_ENERGIES '
                                       'values have to be positive. Input '
                                       'will be ignored.')
                        rpars.setHaltingLevel(1)

            if len(fl) != 3:
                logger.warning('PARAMETERS file: THEO_ENERGIES '
                               'appears to have no value.')
                rpars.setHaltingLevel(1)
                continue
            if defined < 3:
                rpars.THEO_ENERGIES = fl
                continue
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
                rpars.THEO_ENERGIES = fl
            else:
                logger.warning('PARAMETERS file: Unexpected input for '
                               'THEO_ENERGIES. Input will be ignored.')
                rpars.setHaltingLevel(1)
        elif param == 'V0_REAL':
            if value.lower() == 'rundgren':
                try:
                    setTo = [float(other_values[i]) for i in range(4)]
                except (ValueError, IndexError):
                    message = ("Could not parse constants for Rundgren-type "
                               "function.")
                    raise ParameterCustomError(param, custom_message=message)
                    logger.warning(
                    rpars.setHaltingLevel(1)
            else:
                setTo = re.sub("(?i)EE", "EEV+workfn", right_side)
            if isinstance(setTo, str):
                setTo = setTo.rstrip()
            rpars.V0_REAL = setTo
    logger.setLevel(loglevel)


def modifyPARAMETERS(rp, modpar, new="", comment="", path="",
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
    path : str or path-like, optional
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
    file = _path / "PARAMETERS"
    oriname = f"PARAMETERS_ori_{rp.timestamp}"
    ori = _path / oriname
    if oriname not in rp.manifest and not suppress_ori and file.is_file():
        try:
            shutil.copy2(file, ori)
        except Exception:
            logger.error(
                "modifyPARAMETERS: Could not copy PARAMETERS file "
                "to PARAMETERS_ori. Proceeding, original file will be lost.")
        rp.manifest.append(oriname)
    if "PARAMETERS" not in rp.manifest and _path == Path():
        rp.manifest.append("PARAMETERS")
    output = ""
    headerPrinted = False

    try:
        with open(file, "r") as rf:
            plines = rf.readlines()
    except FileNotFoundError:
        plines = []
    except Exception:
        logger.error("Error reading PARAMETERS file.")
        raise
    found = False
    for line in plines:
        if "! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #" in line:
            headerPrinted = True
        valid = False
        param = ""
        stripped_line = strip_comments(line)
        if "=" in stripped_line:
            # Parameter is defined left of "="
            param, *_ = stripped_line.split('=')
            if param:
                valid = True
                param, *_ = param.split()
        else:
            for p in ["STOP", "SEARCH_KILL"]:
                if stripped_line.upper().startswith(p):
                    valid = True
                    param = p
        if valid and param == modpar:
            found = True
            if new and f"{modpar} = {new.strip()}" == stripped_line:
                output += line
            elif new:
                output += f"!{line[:-1]} ! line automatically changed to:\n"
                if not include_left:
                    output += f"{modpar} = "
                output += f"{new} ! {comment}\n"
            else:
                comment = comment or "line commented out automatically"
                output += f"!{line.rstrip():<34} ! {comment}\n"
        else:
            output += line
    if new and not found:
        if not headerPrinted:
            output += """

! ######################################################
! #  THE FOLLOWING LINES WERE GENERATED AUTOMATICALLY  #
! ######################################################
"""
        output += f"\n{modpar} = {new}"
        if comment:
            output += f" ! {comment}"
    try:
        with open(file, "w") as wf:
            wf.write(output)
    except Exception:
        logger.error("modifyPARAMETERS: Failed to write PARAMETERS file.")
        raise
