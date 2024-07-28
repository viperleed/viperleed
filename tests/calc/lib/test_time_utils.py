"""Tests for module time_utils.py of viperleed/calc/lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-26'
__license__ = 'GPLv3+'

from collections import namedtuple
from dataclasses import dataclass
import sys
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


@dataclass(order=True, frozen=True)
class MockCountable:
    """Come object that can be counted."""

    value: float = 0.0

    def __add__(self, other):
        """Add this countable with another one."""
        try:
            return MockCountable(self.value + other.value)
        except (TypeError, ValueError, AttributeError):
            return NotImplemented



def get_tester(tester_name):
    """Return one of the test classes below from its name."""
    this_module = sys.modules[__name__]
    return getattr(this_module, tester_name)


class TestExecutionTimer:
    """Collection of tests for the ExecutionTimer class."""

    timer_cls = ExecutionTimer
    timer_args = tuple()

    def patch_timer(self, patch, *args, **kwargs):
        """Return a fake time module and a patched instance of timer_cls."""
        args = args or self.timer_args
        mock = MockTimer(start_time=kwargs.get('started_at', 0))
        patch.setattr(self.timer_cls, 'now', mock)
        timer = self.make_timer(*args, **kwargs)
        return mock, timer

    def make_timer(self, *args, **kwargs):
        """Return a time of the type tested by this class."""
        args = args or self.timer_args
        return self.timer_cls(*args, **kwargs)

    def test_init(self):
        """Check correct initialization of the timer."""
        timer = self.make_timer()
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
            mock, timer = self.patch_timer(patch)
            mock.sleep(sleep)
            elapsed_time = timer.how_long(as_string=as_string)
        assert isinstance(elapsed_time, type(expect))
        assert elapsed_time == expect

    def test_restart(self, monkeypatch):
        """Check correct outcome of calling restart."""
        with monkeypatch.context() as patch:
            mock, timer = self.patch_timer(patch, started_at=12.3)
            mock.sleep(15)
            assert timer.how_long() == 15
            mock.advance(27)
            timer.restart()
            assert timer.started_at == 54.3

    _other_timer = {
        ExecutionTimer: 'TestExecutionTimer',
        ExpiringOnCountTimer: 'TestExpiringOnCountTimer',
        }

    @parametrize('other,test_cls', _other_timer.items(), ids=_other_timer)
    def test_synchronize_with(self, other, test_cls, monkeypatch):
        """Check synchronization of two timers."""
        other_args = getattr(get_tester(test_cls), 'timer_args')
        with monkeypatch.context() as patch:
            mock, _ = self.patch_timer(patch)
            # Patch also the second timer with the same mock
            patch.setattr(other, 'now', mock)
            timer1 = self.make_timer()
            timer2 = other(*other_args, started_at=-5)
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
        timer = self.make_timer()
        with pytest.raises(TypeError):
            timer.synchronize_with('not a timer')


class TestExpiringOnCountTimer(TestExecutionTimer):
    """Collection of tests for a timer that also counts 'object intervals'."""

    timer_cls = ExpiringOnCountTimer
    timer_args = (5, 0)  # interval, count_start

    def test_init(self):
        """Check correct initialization of the timer."""
        super().test_init()
        timer = self.make_timer()
        assert timer.previous_count == 0  # From self.timer_args

    _InitArgs = namedtuple('_InitArgs', ('interval,count_start'))
    _CountInfo = namedtuple('_CountInfo', ('not_expires,expires,previous'))
    _counts = {
        int: (
            _InitArgs(5, 0),
            _CountInfo(4, 6, 0),
            _CountInfo(7, 12, 6),
            _CountInfo(16, 30, 12),
            ),
        MockCountable: (
            _InitArgs(5.4, 1.2),
            _CountInfo(2.9, 7.1, 1.2),
            _CountInfo(10.4, 19.2, 7.1),
            _CountInfo(23.6, 30, 19.2),
            ),
        str: (
            _InitArgs('aaaa', ''),
            _CountInfo('aa', 'abcd', ''),
            _CountInfo('aaaaaaa', 'abcplussomemore', 'abcd'),
            _CountInfo('aaaab', 'z', 'abcplussomemore'),
            ),
        }

    @parametrize(type_=_counts)
    def test_count_expired(self, type_):
        """Check appropriate detection of count expiration."""
        init_args, *_check = self._counts[type_]
        timer = self.timer_cls(interval=type_(init_args.interval),
                               count_start=type_(init_args.count_start))
        for not_expires, expires, previous in _check:
            assert not timer.count_expired(type_(not_expires))
            assert timer.count_expired(type_(expires))
            assert timer.previous_count == type_(previous)

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
