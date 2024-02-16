"""Package slab of viperleed.tleedmlib.classes.

Created on 2023-02-21

@author: Michele Riva
@author: Florian Kraushofer

This package is a result of the refactor of what used to be the slab
module. It is now split into multiple submodule to ease maintenance.
The API of the package remains unchanged.
"""

from .bulk_slab import BulkSlab
from .slab_errors import AlreadyMinimalError
from .slab_errors import AtomsTooCloseError
from .slab_errors import EmptySlabError
from .slab_errors import InvalidUnitCellError
from .slab_errors import MissingBulkSlabError
from .slab_errors import MissingElementsError
from .slab_errors import MissingLayersError
from .slab_errors import MissingSublayersError
from .slab_errors import NoBulkRepeatError
from .slab_errors import NotEnoughVacuumError
from .slab_errors import NoVacuumError
from .slab_errors import SlabError
from .slab_errors import TooFewLayersError
from .slab_errors import VacuumError
from .slab_errors import WrongVacuumPositionError
from .surface_slab import SurfaceSlab as Slab
