"""Package rparams of viperleed.calc.classes.

This package is a result of the refactor of the rparams.py
module, originally by @fkraushofer. Contains definition of
Rparams, the class containing most of the information during
a viperleed.calc run. It also defines specific classes for
not-so-simple user parameters (in package special), used
as attributes of Rparams. Defaults for the parameters are
defined in module defaults. Their variability limits in
module limits. Some of the special parameters define their
own defaults and limits.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-23'
__license__ = 'GPLv3+'

from .domain_params import DomainParameters
from .rparams import Rparams
from .special.energy_range import EnergyRange, IVShiftRange, TheoEnergies
from .special.l_max import LMax
from .special.layer_cuts import LayerCuts
from .special.search_cull import SearchCull
from .special.symmetry_eps import SymmetryEps
