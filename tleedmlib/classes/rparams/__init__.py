"""Package rparams of viperleed.tleedmlib.classes.

Created on 2023-10-23

@author: Michele Riva (@michele-riva)
@author: Florian Kraushofer (@fkraushofer)
@author: Alexander M. Imre (@amimre)

This package is a result of the refactor of the rparams.py
module, originally by @fkraushofer. Contains definition of
Rparams, the class containing most of the information during
a viperleed.tleedm run. It also defines specific classes for
not-so-simple user parameters (in package special), used as
attributes of Rparams. Defaults for the simple parameters are
defined in module _defaults. Their variability limits in module
_limits. Special parameters define their own defaults and limits.
"""

# Important note: import first stuff from .special, as it is used
# in _rparams and would otherwise lead to cyclic import issues
from .special.layer_cuts import LayerCuts
from ._domain_params import DomainParameters
from ._rparams import Rparams
