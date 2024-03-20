"""Helper module tags of viperleed.tests.calc.

Defines the CaseTag class, an IntEnum useful as tags for @case
functions. This used to be part of viperleed.tests.helpers.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__created__ = '2024-03-19'

from enum import IntEnum, auto


class CaseTag(IntEnum):
    """Enumeration of tags to use for cases."""
    BULK = auto()
    BULK_PROPERTIES = auto()
    LAYER_INFO = auto()
    NEAREST_NEIGHBOURS = auto()
    NEED_ROTATION = auto()
    NO_INFO = auto()
    NON_MINIMAL_CELL = auto()
    RAISES = auto()
    SURFACE_ATOMS = auto()
    THICK_BULK = auto()
    VACUUM_GAP_MIDDLE = auto()  # Not at top
    VACUUM_GAP_SMALL = auto()   # < 5A
    VACUUM_GAP_ZERO = auto()
