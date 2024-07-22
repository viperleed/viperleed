"""Tests for module viperleed.calc.bookkeeper.mode."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-22'
__license__ = 'GPLv3+'


from viperleed.calc.bookkeeper.mode import BookkeeperMode


def test_bookkeeper_mode_enum():
    """Check values of bookkeeper mode Enum."""
    assert BookkeeperMode.ARCHIVE is BookkeeperMode('archive')
    assert BookkeeperMode.CLEAR is BookkeeperMode('clear')
    assert BookkeeperMode.DISCARD is BookkeeperMode('discard')
    assert BookkeeperMode.DISCARD_FULL is BookkeeperMode('discard_full')
