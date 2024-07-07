"""Package slab of viperleed.calc.classes.

This package is a result of the refactor of what used to be the slab
module. It is now split into multiple submodule to ease maintenance.
The API of the package remains unchanged.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-02-21'
__license__ = 'GPLv3+'

from .bulk_slab import BulkSlab
from .errors import AlreadyMinimalError
from .errors import AtomsTooCloseError
from .errors import EmptySlabError
from .errors import InvalidUnitCellError
from .errors import MissingBulkSlabError
from .errors import MissingElementsError
from .errors import MissingLayersError
from .errors import MissingSublayersError
from .errors import NoBulkRepeatError
from .errors import NoVacuumError
from .errors import NotEnoughVacuumError
from .errors import SlabError
from .errors import TooFewLayersError
from .errors import VacuumError
from .errors import WrongVacuumPositionError
from .surface_slab import SurfaceSlab as Slab
