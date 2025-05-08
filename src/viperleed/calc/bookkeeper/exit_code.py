"""Module exit_code of viperleed.calc.bookkeeper.

Defines the BookkeeperExitCode enumeration for the result of
calling Bookkeeper.run.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-05'  # Used to be in bookkeeper.py
__license__ = 'GPLv3+'

from enum import IntEnum


class BookkeeperExitCode(IntEnum):
    """Exit code of the bookkeeper."""
    SUCCESS = 0
    NOTHING_TO_DO = -1
    FAIL = 1

    @classmethod
    def from_codes(cls, exit_codes):
        """Return an overall exit code from multiple ones."""
        if not exit_codes:
            raise ValueError('At least one exit code needed.')
        if any(c is cls.FAIL for c in exit_codes):
            return cls.FAIL
        if all(c is cls.NOTHING_TO_DO for c in exit_codes):
            return cls.NOTHING_TO_DO
        return cls.SUCCESS
