"""Tests for module time_utils.py of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-26'
__license__ = 'GPLv3+'

from collections import namedtuple
import datetime
import sys
import time

import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.dataclass_utils import frozen
from viperleed.calc.lib.time_utils import DateTimeFormat
from viperleed.calc.lib.time_utils import ExecutionTimer
from viperleed.calc.lib.time_utils import ExpiringOnCountTimer
from viperleed.calc.lib.time_utils import ExpiringTimer
from viperleed.calc.lib.time_utils import ExpiringTimerWithDeadline


def get_tester(tester_name):
    """Return one of the test classes below from its name."""
    this_module = sys.modules[__name__]
    return getattr(this_module, tester_name)


def patch_time_module(patch, current_time=0):
    """Temporarily replace the time module of viperleed.calc.lib.time_utils."""
    mock_module = MockTime(start_time=current_time)
    patch.setattr('viperleed.calc.lib.time_utils.time', mock_module)
    return mock_module


class MockTime:
    """A fake version of the time module."""

    def __init__(self, start_time):
        """Initialize fake timer at `start_time`."""
        self._current_time = start_time

    def __call__(self):
        """Return the current time."""
        return self._current_time

    def gmtime(self):
        """Return the current time in the GMT zone."""
        now = datetime.datetime.fromtimestamp(self.now())
        return now.timetuple()

    def localtime(self):
        """Return the current, local time."""
        return self.gmtime()

    def now(self):
        """Return the current time."""
        return self()

    def sleep(self, seconds):
        """Sleep for a few seconds."""
        self._current_time += seconds

    strftime = time.strftime


@frozen(order=True)
class MockCountable:
    """Come object that can be counted."""

    value: float = 0.0

    def __add__(self, other):
        """Add this countable with another one."""
        try:
            return MockCountable(self.value + other.value)
        except (TypeError, ValueError, AttributeError):
            return NotImplemented


class TestExecutionTimer:
    """Collection of tests for the ExecutionTimer class."""

    timer_args = tuple()
    timer_cls = ExecutionTimer

    def patch_timer(self, patch, *args, **kwargs):
        """Return a fake time module and a patched instance of timer_cls."""
        args = args or self.timer_args
        mock = MockTime(start_time=kwargs.get('started_at', 0))
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
        start, sleep_first, sleep_second = 12.3, 15, 27
        with monkeypatch.context() as patch:
            mock, timer = self.patch_timer(patch, started_at=start)
            mock.sleep(sleep_first)
            assert timer.how_long() == sleep_first
            mock.sleep(sleep_second)
            timer.restart()
            assert timer.started_at == start + sleep_first + sleep_second

    _other_timer = {
        ExecutionTimer: 'TestExecutionTimer',
        ExpiringOnCountTimer: 'TestExpiringOnCountTimer',
        ExpiringTimer: 'TestExpiringTimer',
        }

    @parametrize('other,test_cls', _other_timer.items(), ids=_other_timer)
    def test_synchronize_with(self, other, test_cls, monkeypatch):
        """Check synchronization of two timers."""
        other_args = getattr(get_tester(test_cls), 'timer_args')
        timer2_delay, sleep_for = -5, 15
        with monkeypatch.context() as patch:
            mock, _ = self.patch_timer(patch)
            # Patch also the second timer with the same mock
            patch.setattr(other, 'now', mock)
            timer1 = self.make_timer()
            timer2 = other(*other_args, started_at=timer2_delay)
            mock.sleep(sleep_for)
            assert timer1.started_at != timer2.started_at
            assert timer1.how_long() == sleep_for
            assert timer2.how_long() == sleep_for - timer2_delay
            mock.sleep(sleep_for)
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

    _init_args = namedtuple('_init_args', ('interval,count_start'))
    timer_args = _init_args(interval=5, count_start=0)
    timer_cls = ExpiringOnCountTimer

    def test_init(self):
        """Check correct initialization of the timer."""
        super().test_init()
        timer = self.make_timer()
        assert timer.previous_count == self.timer_args.count_start

    _count_info = namedtuple('_count_info', ('not_expires,expires,previous'))
    _counts = {
        int: (
            _init_args(5, 0),
            _count_info(4, 6, 0),
            _count_info(7, 12, 6),
            _count_info(16, 30, 12),
            ),
        MockCountable: (
            _init_args(5.4, 1.2),
            _count_info(2.9, 7.1, 1.2),
            _count_info(10.4, 19.2, 7.1),
            _count_info(23.6, 30, 19.2),
            ),
        str: (
            _init_args('aaaa', ''),
            _count_info('aa', 'abcd', ''),
            _count_info('aaaaaaa', 'abcplussomemore', 'abcd'),
            _count_info('aaaab', 'z', 'abcplussomemore'),
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


class TestExpiringTimer(TestExecutionTimer):
    """Tests for a timer that periodically expires in time."""

    _init_args = namedtuple('_init_args', ('interval'))
    timer_args = _init_args(interval=5)
    timer_cls = ExpiringTimer

    def test_init(self):
        """Check correct initialization."""
        super().test_init()
        timer = self.make_timer()
        assert timer.interval == self.timer_args.interval

    def test_has_expired(self, monkeypatch):
        """Check correct expiration of the timer."""
        interval = self.timer_args.interval
        with monkeypatch.context() as patch:
            mock, timer = self.patch_timer(patch, expire_once=True)
            # Pylint seems to have problems inferring timer comes
            # from self.timer_cls, which is not an ExecutionTimer.
            # It thus emits various no-member (E1101).
            # pylint: disable-next=no-member
            assert timer.has_expired()  # because of expire_once
            mock.sleep(0.25*interval)
            assert not timer.has_expired()  # pylint: disable=no-member
            mock.sleep(1.05*interval)
            assert timer.has_expired()      # pylint: disable=no-member


class TestExpiringTimerWithDeadline(TestExecutionTimer):
    """Tests for a periodically expiring timer that also 'stops' for good."""

    _init_args = namedtuple('_init_args', ('interval,deadline'))
    timer_args = _init_args(interval=5, deadline=32)
    timer_cls = ExpiringTimerWithDeadline

    def test_has_reached_deadline(self, monkeypatch):
        """Check correct detection of deadline-reached condition."""
        with monkeypatch.context() as patch:
            mock, timer = self.patch_timer(patch, started_at=19)
            # Pylint seems to have problems inferring timer comes
            # from self.timer_cls, which is not an ExecutionTimer.
            # It thus emits various no-member (E1101).
            # pylint: disable=no-member
            assert not timer.has_expired()
            assert not timer.has_reached_deadline()
            mock.sleep(12)
            assert timer.has_expired()               # interval is 5
            assert not timer.has_reached_deadline()  # deadline is 32
            mock.sleep(12)
            assert timer.has_expired()
            assert not timer.has_reached_deadline()
            mock.sleep(12)
            assert timer.has_expired()
            assert timer.has_reached_deadline()


class TestDateTimeFormat:  # pylint: disable=too-few-public-methods
    """Tests for the DateTimeFormat enumeration class."""

    mock_time_stamp = time.mktime((2023, 1, 1, 12, 1, 38, 6, 1, 0))
    _formatted = {
        'file_suffix': '230101-120138',
        'iso': '2023-01-01 12:01:38',
        'log_contents': '2023-01-01 12:01:38',
        'time': '12:01:38',
        }

    @parametrize('fmt,expect', _formatted.items(), ids=_formatted)
    @parametrize(use_gmt=(True, False))
    def test_datetime_format(self, fmt, expect, use_gmt, monkeypatch):
        """Check correct date-time formatting."""
        with monkeypatch.context() as patch:
            patch_time_module(patch, current_time=self.mock_time_stamp)
            fmt_enum = DateTimeFormat[fmt.upper()]
            assert fmt_enum.now(use_gmt=use_gmt) == expect
