"""Tests for module viperleed.calc.bookkeeper.mode."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-22'
__license__ = 'GPLv3+'

from pytest_cases import parametrize

from viperleed.calc.bookkeeper.mode import BookkeeperMode as Mode


def test_bookkeeper_mode_enum():
    """Check values of bookkeeper mode Enum."""
    assert Mode.ARCHIVE is Mode('archive')
    assert Mode.CLEAR is Mode('clear')
    assert Mode.DISCARD is Mode('discard')
    assert Mode.DISCARD_FULL is Mode('discard_full')
    assert Mode.FIX is Mode('fix')


long_flags = {
    Mode.ARCHIVE: '--archive',
    Mode.CLEAR: '--clear',
    Mode.DISCARD: '--discard',
    Mode.DISCARD_FULL: '--discard-full',
    Mode.FIX: '--fix',
    }
short_flags = {
    Mode.ARCHIVE: '-a',
    Mode.CLEAR: '-c',
    Mode.DISCARD: '-d',
    Mode.DISCARD_FULL: '-df',
    Mode.FIX: None,
    }

@parametrize('mode,expect', long_flags.items())
def test_long_flag(mode, expect):
    """Check the .long_flag property."""
    assert getattr(mode, 'long_flag') == expect


@parametrize('mode,long_flag', long_flags.items())
def test_flags(mode, long_flag):
    """Check the .flags property."""
    short_flag = short_flags[mode]
    expect = (
        (short_flag, long_flag) if short_flag is not None
        else (long_flag,)
        )
    assert getattr(mode, 'flags') == expect
