"""Test helper module tag of tests/calc/bookkeeper.

Defines the BookkeeperTag class, useful for marking test cases.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-19'
__license__ = 'GPLv3+'

from enum import IntEnum, auto


class BookkeeperTag(IntEnum):
    """Tags for cases for the bookkeeper utility."""
    AUTO_FIX = auto()        # Entry fields can be auto-fixed
    AUTO_FIX_ENTRY = auto()  # The whole entry can be auto-fixed
    CANT_FIX = auto()        # Entry cannot be auto-fixed
    EMPTY = auto()           # history.info file is empty
    ENTRY = auto()           # Use only for HistoryInfoEntry
    BOOKKEEPER = auto()      # Use only for bookkeeper.bookkeeper tests
    HISTORY = auto()         # Use only for history.info file tests
    MISS_MANDATORY = auto()  # Some important entry fields are missing
    MULTI_ENTRY = auto()     # Contains more than one entry
    NEEDS_NO_FIX = auto()    # Entry fields do not require any fixing
    NO_ISSUES = auto()       # Case has no problems at all
    NOT_PRESERVED = auto()   # Contents are not maintained
    OLD = auto()             # Entry has some legacy format
    RAISES = auto()          # Entry causes exceptions when parsed
