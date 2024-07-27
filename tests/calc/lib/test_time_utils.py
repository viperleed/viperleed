"""Tests for module time_utils.py of viperleed/calc/lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-26'
__license__ = 'GPLv3+'

import time

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.lib.time_utils import ExecutionTimer
from viperleed.calc.lib.time_utils import ExpiringOnCountTimer
from viperleed.calc.lib.time_utils import ExpiringTimer
from viperleed.calc.lib.time_utils import ExpiringTimerWithDeadline
from viperleed.calc.lib.time_utils import _elapsed_time_as_str


class MockTimer:
    """Helper for mocking the passage of time."""

    def __init__(self, start_time):
        """Initialize fake timer at `start_time`."""
        self._current_time = start_time

    def __call__(self):
        """Return the current time."""
        return self._current_time

    def sleep(self, seconds):
        """Sleep for a few seconds."""
        self.advance(seconds)

    def advance(self, seconds):
        """Add a few seconds to the current time."""
        self._current_time += seconds


def patch_timer(patch, timer_cls, *args, **kwargs):
    """Return a fake time module and a patched instance of timer_cls."""
    mock = MockTimer(start_time=kwargs.get('started_at', 0))
    patch.setattr(timer_cls, 'now', mock)
    timer = timer_cls(*args, **kwargs)
    return mock, timer


class TestExecutionTimer:
    """Collection of tests for the ExecutionTimer class."""

    def test_init(self):
        """Check correct initialization of the timer."""
        timer = ExecutionTimer()
        assert timer.started_at >= 0

    _sleep = {
        1.234: 1.234,
        0.1: '0.10 seconds',
        15.33: '15.33 seconds',
        65.4: '1:05 minutes',
        193.4: '3:13 minutes',
        6214.8: '1:43 hours',
        }

    @parametrize('sleep,expect', _sleep.items(), ids=_sleep.values())
    def test_how_long_str(self, sleep, expect, monkeypatch):
        """Check the expected result of checking the elapsed time interval."""
        as_string = isinstance(expect, str)
        with monkeypatch.context() as patch:
            mock, timer = patch_timer(patch, ExecutionTimer)
            mock.sleep(sleep)
            elapsed_time = timer.how_long(as_string=as_string)
        assert isinstance(elapsed_time, type(expect))
        assert elapsed_time == expect

    def test_restart(self, monkeypatch):
        """Check correct outcome of calling restart."""
        with monkeypatch.context() as patch:
            mock, timer = patch_timer(patch, ExecutionTimer, started_at=12.3)
            mock.sleep(15)
            assert timer.how_long() == 15
            mock.advance(27)
            timer.restart()
            assert timer.started_at == 54.3

    def test_synchronize_with_same(self, monkeypatch):
        """Check synchronization of two timers."""
        with monkeypatch.context() as patch:
            mock, _ = patch_timer(patch, ExecutionTimer)
            timer1 = ExecutionTimer()
            timer2 = ExecutionTimer(started_at=-5)
            mock.sleep(15)
            assert timer1.started_at != timer2.started_at
            assert timer1.how_long() == 15
            assert timer2.how_long() == 20
            mock.sleep(15)
            timer2.synchronize_with(timer1)
            assert timer1.started_at == timer2.started_at
            assert timer1.how_long() == timer2.how_long()

    def test_synchronize_typeerror(self):
        """Check complaints when synchronizing with a non-timer."""
        timer = ExecutionTimer()
        with pytest.raises(TypeError):
            timer.synchronize_with('not a timer')

# def test_expiring_on_count_timer_init():
    # timer = ExpiringOnCountTimer(interval=5, count_start=0)
    # assert timer.previous_count == 0

# def test_expiring_on_count_timer_count_expired():
    # timer = ExpiringOnCountTimer(interval=5, count_start=0)
    # assert not timer.count_expired(4)
    # assert timer.count_expired(6)

# def test_expiring_timer_init():
    # timer = ExpiringTimer(interval=5)
    # assert timer.interval == 5

# def test_expiring_timer_has_expired():
    # mock_time = MockTimer(start_time=0)
    # ExecutionTimer.now = mock_time
    # timer = ExpiringTimer(interval=5)
    # assert not timer.has_expired()
    # mock_time.advance(6)
    # assert timer.has_expired()

# def test_expiring_timer_restart():
    # timer = ExpiringTimer(interval=5)
    # time.sleep(0.1)
    # timer.restart()
    # assert timer.started_at >= 0

# def test_expiring_timer_with_deadline_init():
    # timer = ExpiringTimerWithDeadline(interval=5, deadline=10)
    # assert timer.interval == 5

# def test_expiring_timer_with_deadline_has_reached_deadline():
    # mock_time = MockTimer(start_time=0)
    # ExecutionTimer.now = mock_time
    # timer = ExpiringTimerWithDeadline(interval=5, deadline=10)
    # mock_time.advance(9)
    # assert not timer.has_reached_deadline()
    # mock_time.advance(2)
    # assert timer.has_reached_deadline()

# def test_elapsed_time_as_str():
    # assert _elapsed_time_as_str(3661) == "1:01 hours"
    # assert _elapsed_time_as_str(61) == "1:01 minutes"
    # assert _elapsed_time_as_str(0.01) == "0.01 seconds"
