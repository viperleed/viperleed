# -*- coding: utf-8 -*-
"""Moule _limits of viperleed.calc.classes.rparams.

Created on 2023-10-23, originally Jun 13 2019

@author: Florian Kraushofer (@fkraushofer)

Defines the limiting values of 'simple' user PARAMETERS. Not-so-simple
parameters, which are defined as their own classes in package special,
also take care of their own default values. This module was originally
part of the rparams.py module, refactored by Michele Riva in Oct 2023.
"""

import logging

_LOGGER = logging.getLogger('tleedm.rparams')

                                                                                # TODO: fill dict of parameter limits here (e.g. LMAX etc.)
# parameter limits
# either tuple of (min, max) or list of allowed values                          # TODO: allowed would be cleaner as set. It's not great that things are mixed though. Would be better to have a separate global
PARAM_LIMITS = {
    'LMAX': (1, 18),
    'INTPOL_DEG': ['3', '5'],
    }
