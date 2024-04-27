"""Moule limits of viperleed.calc.classes.rparams.

Defines the limiting values of 'simple' user PARAMETERS. Not-so-simple
parameters, which are defined as their own classes in package special,
also take care of their own default values. This module was originally
(2019-06-13) part of the rparams.py module, refactored by Michele Riva
in Oct 2023.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-23'
__license__ = 'GPLv3+'

                                                                                # TODO: fill dict of parameter limits here (e.g. LMAX etc.)
# parameter limits
# either tuple of (min, max) or list of allowed values                          # TODO: allowed would be cleaner as set. It's not great that things are mixed though. Would be better to have a separate global
PARAM_LIMITS = {
    'LMAX': (1, 18),
    'INTPOL_DEG': ['3', '5'],
    }
