"""Package slab of viperleed.tleedmlib.classes.

Created on 2023-02-21

@author: Michele Riva
@author: Florian Kraushofer

This package is a result of the refactor of what used to be the slab
module. It is now split into multiple submodule to ease maintenance.
The API of the package remains unchanged.
"""

from .bulk_slab import BulkSlab
from .slab_errors import SlabError, AlreadyMinimalError
from .surface_slab import SurfaceSlab as Slab
