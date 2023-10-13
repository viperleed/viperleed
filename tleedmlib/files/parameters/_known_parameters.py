# -*- coding: utf-8 -*-
"""Module parameters of viperleed.tleedmlib.

Created on Tue Aug 18 16:56:39 2020

@author: Florian Kraushofer
@author: Alexander M. Imre
@author: Michele Riva

Initial version by Florian Kraushofer in 2020, major rewrite by
Alexander Imre & Michele Riva in June 2023.

Functions for reading from and writing to the PARAMETERS file.
"""

import ast
from collections.abc import Sequence
from dataclasses import dataclass, field
from functools import partialmethod
import logging
from pathlib import Path
import re
import shutil

import numpy as np

from viperleed.tleedmlib.base import (strip_comments, splitSublists,
                                      readVector, readIntRange,
                                      recombineListElements)
from viperleed.tleedmlib.classes import rparams
from viperleed.tleedmlib.files.woods_notation import readWoodsNotation
from viperleed.tleedmlib import periodic_table
from viperleed.tleedmlib.sections._sections import TLEEDMSection as Section

from .errors import (
    ParameterError, ParameterValueError, ParameterParseError,
    ParameterIntConversionError, ParameterFloatConversionError,
    ParameterBooleanConversionError, ParameterNotRecognizedError,
    ParameterNumberOfInputsError, ParameterRangeError,
    ParameterUnknownFlagError, ParameterNeedsFlagError
    )


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
    'R_FACTOR_LEGACY', 'R_FACTOR_SMOOTH', 'R_FACTOR_TYPE',
    'SCREEN_APERTURE', 'SEARCH_BEAMS', 'SEARCH_CONVERGENCE', 'SEARCH_CULL',
    'SEARCH_MAX_GEN', 'SEARCH_POPULATION', 'SEARCH_START', 'SITE_DEF',
    'SUPERLATTICE', 'SUPPRESS_EXECUTION', 'SYMMETRIZE_INPUT', 'SYMMETRY_BULK',
    'SYMMETRY_CELL_TRANSFORM', 'SYMMETRY_EPS', 'SYMMETRY_FIND_ORI',
    'SYMMETRY_FIX', 'TENSOR_INDEX', 'TENSOR_OUTPUT', 'THEO_ENERGIES',
    'TL_VERSION', 'TL_IGNORE_CHECKSUM',
    'T_DEBYE', 'T_EXPERIMENT', 'V0_IMAG', 'V0_REAL',
    'V0_Z_ONSET', 'VIBR_AMP_SCALE', 'ZIP_COMPRESSION_LEVEL',
    )



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
    'compression_level': 'ZIP_COMPRESSION_LEVEL',
    'compression': 'ZIP_COMPRESSION_LEVEL',
    }


for known_param in _KNOWN_PARAMS:
    _PARAM_ALIAS[known_param.lower().replace('_', '')] = known_param
