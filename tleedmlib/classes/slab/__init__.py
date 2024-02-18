"""Package slab of viperleed.tleedmlib.classes.

Created on 2023-02-21

@author: Michele Riva
@author: Florian Kraushofer

This package is a result of the refactor of what used to be the slab
module. It is now split into multiple submodule to ease maintenance.
The API of the package remains unchanged.
"""

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
from .errors import NotEnoughVacuumError
from .errors import NoVacuumError
from .errors import SlabError
from .errors import TooFewLayersError
from .errors import VacuumError
from .errors import WrongVacuumPositionError
from .surface_slab import SurfaceSlab as Slab
