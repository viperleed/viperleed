"""Moule defaults of viperleed.calc.classes.rparams.

Defines the default values of 'simple' user PARAMETERS. Not-so-simple
parameters, which are defined as their own classes in package special,
also take care of their own default values. This module was originally
(2019-06-13) part of the rparams.py module, refactored by Michele Riva
in Oct 2023.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-23'
__license__ = 'GPLv3+'

import logging


# The name we expect for folders containing TensErLEED Fortran code
TENSERLEED_FOLDER_NAME = 'TensErLEED'

# Notice that we cannot use a module-level global object(), as this
# module may be imported a number of times when using multiprocessing
NO_VALUE = None


# Notice that the defaults in here that may be mutated during execution
# are saved as immutable types to prevent inadvertent modification of
# this global, and are rather converted to their mutable equivalent
# in the relevant places. The only difference is dictionaries. Copies
# are used for them.
DEFAULTS = {
    # BASIC PARAMETERS
    'FILAMENT_WF': {
        'lab6': 2.65,  # This is the default if nothing is given
        'w': 4.5,
        },
    # GAUSSIAN_WIDTH and GAUSSIAN_WIDTH_SCALING are set by parameter
    # SEARCH_CONVERGENCE
    'GAUSSIAN_WIDTH' : 0.1,
    'GAUSSIAN_WIDTH_SCALING' : 0.5,
    'LOG_LEVEL' : {
        NO_VALUE: logging.INFO,
        'debug': logging.DEBUG,
        'v' : 5, 'verbose' : 5,
        'vv' : 1, 'vverbose' : 1,
        },
    'OPTIMIZE': {  # settings for fd optimization
        'which': 'none',
        'step': 0.,
        'minpoints': 4,
        'maxpoints': 10,
        'convergence': 0.,
        'maxstep': 0.,
        },
    'PHASESHIFT_EPS': {
        'r': 0.1,
        'n': 0.05,
        'd': 0.02, # default value
        'f': 0.01,
        },
    'RUN': (0, 1, 2, 3),
    'SEARCH_EVAL_TIME': 60,  # time interval between reads of SD.TL,            # TODO: should be dynamic?
    'SEARCH_MAX_DGEN': {'all': 0, 'best': 0, 'dec': 100},
    'SYMMETRY_FIX': '',
    'THETA': 0,   # perpendicular incidence
    'TL_VERSION': None,
    'PHI': 0,     # not needed in case of perpendicular incidence
    'ZIP_COMPRESSION_LEVEL': 2,

    # SPECIAL PARAMETERS
    'IV_SHIFT_RANGE': (-3, 3, NO_VALUE),  # NO_VALUE step: from data
    'LAYER_CUTS': 'dz(1.2)',
    'LMAX': (NO_VALUE, NO_VALUE),
    'SEARCH_CULL': (0.1, 'genetic'),
    'SYMMETRY_EPS': 0.1,  # z always equal to in-plane
    'THEO_ENERGIES': (NO_VALUE, NO_VALUE, NO_VALUE),
    'THEO_ENERGIES - no experiments': (20, 800, 3),
    }
