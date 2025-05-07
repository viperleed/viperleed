"""Tests for BookkeeperExitCode of viperleed.calc.bookkeeper.bookkeeper."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-04-08'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.bookkeeper import BookkeeperExitCode as ExitCode


def test_value():
    """Check the .value of all exit codes."""
    assert not ExitCode.SUCCESS.value
    assert ExitCode.NOTHING_TO_DO.value < 0
    assert ExitCode.FAIL.value > 0


class TestFromCodes:
    """Tests for the from_codes class method."""

    _valid = {
        'all nothing': ((ExitCode.NOTHING_TO_DO,)*3, ExitCode.NOTHING_TO_DO),
        'all success': ((ExitCode.SUCCESS,)*5, ExitCode.SUCCESS),
        'one fail': ((ExitCode.SUCCESS, ExitCode.FAIL), ExitCode.FAIL),
        'one nothing': ((ExitCode.SUCCESS, ExitCode.NOTHING_TO_DO),
                        ExitCode.SUCCESS),
        }
    @parametrize('codes,expect', _valid.values(), ids=_valid)
    def test_valid(self, codes, expect):
        """Test the successful combination of multiple exit codes."""
        combined = ExitCode.from_codes(codes)
        assert combined is expect

    _invalid = {
        'no codes': ((), ValueError),
        'not iterable': (object(), TypeError),
        }

    @parametrize('codes,exc', _invalid.values(), ids=_invalid)
    def test_raises(self, codes, exc):
        """Test complaints with invalid codes to be combined."""
        with pytest.raises(exc):
            ExitCode.from_codes(codes)
