"""Module search_backends of viperleed.calc.
"""

__authors__ = (
    "Alexander M. Imre (@amimre)",
)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-04-15"
__license__ = "GPLv3+"

from enum import Enum

class SearchBackend(Enum):
    """Enum for search backends."""

    TENSERLEED = "TensErLEED"
    VLJ = "ViPErLEED-JAX"

    def __str__(self):
        return self.value
