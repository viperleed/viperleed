"""Module known_parameters of viperleed.calc.files.parameters.

Initial version by @fkraushofer in 2020, major rewrite by @amimre
and @michele-riva in June 2023. This module used to be part of
parameters.py. Refactored in October 2023.

Contains definitions of all known parameters and their supported
aliases. It also defines useful functions for retrieving parameter
names used internally when reading/writing/interpreting a PARAMETERS
file.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-08-18'
__license__ = 'GPLv3+'

from difflib import get_close_matches
import logging

from .errors import ParameterNotRecognizedError


_LOGGER_NAME, _ = __name__.rsplit('.', maxsplit=1)
_LOGGER = logging.getLogger(_LOGGER_NAME)


# Allowed parameters with their 'standard' names in alphabetic order
KNOWN_PARAMS = (
    'ATTENUATION_EPS',
    'AVERAGE_BEAMS',
    'BEAM_INCIDENCE',
    'BULKDOUBLING_EPS',
    'BULKDOUBLING_MAX',
    'BULK_LIKE_BELOW',
    'BULK_REPEAT',
    'DOMAIN',
    'DOMAIN_STEP',
    'ELEMENT_MIX',
    'ELEMENT_RENAME',
    'FILAMENT_WF',
    'FORTRAN_COMP',
    'HALTING',
    'INTPOL_DEG',
    'IV_SHIFT_RANGE',
    'LAYER_CUTS',
    'LAYER_STACK_VERTICAL',
    'LMAX',
    'LOG_LEVEL',
    'LOG_SEARCH',
    'N_BULK_LAYERS',
    'N_CORES',
    'OPTIMIZE',
    'PARABOLA_FIT',
    'PHASESHIFT_EPS',
    'PHASESHIFTS_CALC_OLD',
    'PHASESHIFTS_OUT_OLD',
    'PLOT_IV',
    'RUN',
    'R_FACTOR_LEGACY',
    'R_FACTOR_SMOOTH',
    'R_FACTOR_TYPE',
    'SCREEN_APERTURE',
    'SEARCH_BEAMS',
    'SEARCH_CONVERGENCE',
    'SEARCH_CULL',
    'SEARCH_MAX_GEN',
    'SEARCH_POPULATION',
    'SEARCH_START',
    'SITE_DEF',
    'SUPERLATTICE',
    'SUPPRESS_EXECUTION',
    'STOP',
    'SYMMETRIZE_INPUT',
    'SYMMETRY_BULK',
    'SYMMETRY_CELL_TRANSFORM',
    'SYMMETRY_EPS',
    'SYMMETRY_FIND_ORI',
    'SYMMETRY_FIX',
    'TENSOR_INDEX',
    'TENSOR_OUTPUT',
    'THEO_ENERGIES',
    'TL_VERSION',
    'TL_IGNORE_CHECKSUM',
    'T_DEBYE',
    'T_EXPERIMENT',
    'V0_IMAG',
    'V0_REAL',
    'V0_Z_ONSET',
    'VIBR_AMP_SCALE',
    'ZIP_COMPRESSION_LEVEL',
    )


_DEPRECATED = {  # alias: (since version, till version), both inclusive
    'parabolafit': ('0.11.0', None),
}


def _make_alias(param):
    """Return an alias of param."""
    return param.lower().replace('_', '')


# _PARAM_ALIAS keys should be all lowercase, with no underscores
_PARAM_ALIAS = {_make_alias(p): p for p in KNOWN_PARAMS}
_PARAM_ALIAS.update({    # Sort keys alphabetically!
    'bulklike': 'BULK_LIKE_BELOW',
    'bulksymmetry': 'SYMMETRY_BULK',
    'compiler': 'FORTRAN_COMP',
    'compression': 'ZIP_COMPRESSION_LEVEL',
    'compression_level': 'ZIP_COMPRESSION_LEVEL',
    'fdoptimize': 'OPTIMIZE',
    'fdoptimization': 'OPTIMIZE',
    'fortrancompile': 'FORTRAN_COMP',
    'fortrancompiler': 'FORTRAN_COMP',
    'ignorechecksum': 'TL_IGNORE_CHECKSUM',
    'ivplot': 'PLOT_IV',
    'logdebug' : 'LOG_LEVEL',
    'plotrfactor': 'PLOT_IV',
    'plotrfactors': 'PLOT_IV',
    'searchkill': 'STOP',
    })


def from_alias(known_param_or_alias):
    """Return a parameter name, possibly from its alias.

    Parameters
    ----------
    known_param_or_alias : str
        The name of a parameter or one of its aliases.

    Returns
    -------
    known_param : str
        The name of the known parameter.

    Raises
    ------
    TypeError
        If known_param_or_alias is not a string.
    ParameterNotRecognizedError
        If known_param_or_alias is not a known parameter.
    """
    if known_param_or_alias in KNOWN_PARAMS:
        return known_param_or_alias
    try:
        return _PARAM_ALIAS[_make_alias(known_param_or_alias)]
    except KeyError:
        pass
    except AttributeError as exc:  # No .lower or no .replace
        raise TypeError('from_alias: must be string') from exc
    raise ParameterNotRecognizedError(parameter=known_param_or_alias)


def did_you_mean(unknown_param):
    """Return one close match for `unknown_param` from the known ones."""
    try:
        alias = _make_alias(unknown_param)
    except AttributeError as exc:  # No .lower or no .replace
        raise TypeError('from_alias: must be string') from exc
    close_matches = get_close_matches(alias, _PARAM_ALIAS)
    if not close_matches:
        raise ParameterNotRecognizedError(unknown_param)
    return _PARAM_ALIAS[close_matches[0]]


def is_deprecated(param, version):
    """Return whether param is deprecated for a ViPErLEED version."""
    try:
        since, till = _DEPRECATED.get(_make_alias(param))
    except (TypeError, KeyError):
        return False
    if version < since:
        return False
    if till is not None and version > till:
        return False
    return True


def warn_if_deprecated(param, version):
    """Emit a log warning if `param` is deprecated for `version`.

    Parameters
    ----------
    param : str
        Parameter to be checked for deprecation. It may also be
        one of the aliases.
    version : str
        The ViPErLEED version for which deprecation should be
        checked.

    Returns
    -------
    is_deprecated : bool
        True if `param` is deprecated for `version`.
    """
    if not is_deprecated(param, version):
        return False
    _LOGGER.warning(f'PARAMETER {param} is deprecated in ViPErLEED version '
                    f'{version}. This may cause errors or undefined behavior. '
                    'Use at your own risk')
    return True
