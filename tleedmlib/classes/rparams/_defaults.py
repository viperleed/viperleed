# -*- coding: utf-8 -*-
"""Moule _defaults of viperleed.tleedmlib.classes.rparams.

Created on 2023-10-23, originally Jun 13 2019

@author: Florian Kraushofer (@fkraushofer)

Defines the default values of 'simple' user PARAMETERS. Not-so-simple
parameters, which are defined as their own classes in package special,
also take care of their own default values. This module was originally
part of the rparams.py module, refactored by Michele Riva in Oct 2023.
"""

import logging


_LOGGER = logging.getLogger('tleedm.rparams')

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
    'LOG_LEVEL' : {
        NO_VALUE: logging.INFO,
        'debug': logging.DEBUG,
        'v' : 5, 'verbose' : 5,
        'vv' : 1, 'vverbose' : 1,
        },
    'PHASESHIFT_EPS': {
        'r': 0.1,
        'n': 0.05,
        'f': 0.01,  # This is the default if nothing is given
        'e': 0.001,
        },
    'RUN': (0, 1, 2, 3),
    'SEARCH_CULL_TYPE': 'genetic',
    'SEARCH_EVAL_TIME': 60,  # time interval between reads of SD.TL,            # TODO: should be dynamic?
    'SEARCH_MAX_DGEN': {'all': 0, 'best': 0, 'dec': 100},
    'SYMMETRY_FIX': '',
    'THETA': 0,   # perpendicular incidence
    'PHI': 0,     # not needed in case of perpendicular incidence
    'ZIP_COMPRESSION_LEVEL': 2,

    # SPECIAL PARAMETERS
    'IV_SHIFT_RANGE': (-3, 3, NO_VALUE),  # NO_VALUE step: from data
    'LAYER_CUTS': 'dz(1.2)',
    'THEO_ENERGIES': (NO_VALUE, NO_VALUE, NO_VALUE),
    'THEO_ENERGIES - no experiments': (20, 800, 3),
    }
