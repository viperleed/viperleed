"""Module time_utils.py of viperleed.calc.lib.

Contains functions and definitions related to our handling of time
formats, and extensions to the time and datetime standard libraries.
It also defined the ExecutionTimer and its subclasses, to keep track
of the time spent doing stuff.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-07-26'
__license__ = 'GPLv3+'

from datetime import timedelta
from timeit import default_timer  # This is the most accurate available


_ONE_MINUTE = 60
_ONE_HOUR = 60 * _ONE_MINUTE


class ExecutionTimer:
    """A class that keeps track of how long passes."""

    now = default_timer

    def __init__(self, started_at=None):
        """Initialize instance."""
        self._start = 0
        self.restart(started_at=started_at)

    @property
    def started_at(self):
        """Return the reference time for this timer."""
        return self._start

    def how_long(self, as_string=False):
        """Return how long has passed, in seconds.

        Parameters
        ----------
        as_string : bool, optional
            Return `interval` as formatted string, instead of a float.
            Default is False.

        Returns
        -------
        interval : str or float
            The time interval passed since this timer was started.
            When a string (`as_string == True`), it is formatted as:
                - "H:MM hours"     if interval is longer than 1 h
                - "M:SS minutes"   if interval is less more than 1 min
                - "S.mm seconds"   if interval is shorter than 1 min
        """
        interval = self.now() - self.started_at
        return _elapsed_time_as_str(interval) if as_string else float(interval)

    def restart(self, started_at=None):
        """Start this timer again."""
        self._start = self.now() if started_at is None else started_at

    def synchronize_with(self, other):
        """Set this timer's start moment identical to the one of other."""
        try:
            self.restart(started_at=other.started_at)
        except AttributeError as exc:
            raise TypeError('Can only synchronize with ExecutionTimer, '
                            f'not {type(other).__name__!r}') from exc


class ExpiringOnCountTimer(ExecutionTimer):
    """A timer that has a non-time-related interval for expiration."""

    def __init__(self, interval, count_start, started_at=None):
        """Initialize instance with a time interval.

        Parameters
        ----------
        interval : object
            How often should this timer expire. Should support
            addition with `count_start`.
        count_start : object
            The starting value for the stuff to be "counted". Must
            support subtraction with other objects to be counted.
        started_at : float, optional
            When this timer was started. If not given, the timer
            is started right now. Default is None.

        Returns
        -------
        None.
        """
        super().__init__(started_at=started_at)
        self._interval = interval
        self._count_started = count_start
        self._previous_count = count_start

    @property
    def previous_count(self):
        """Return the starting counts before the last expiration."""
        return self._previous_count

    def count_expired(self, count_new):
        """Return whether count_new causes this timer to expire.

        Parameters
        ----------
        count_new : object
            The quantity to use for checking whether this interval
            has expired. Should support comparing to self.interval
            via ordering operators. If this is the case, `count_new`
            is saved internally as the starting point for the next
            expiration.

        Returns
        -------
        expired : bool
            Whether count_new made the timer expire.
        """
        prev_count = self._count_started
        expired = count_new > prev_count + self._interval
        if expired:
            self._previous_count = prev_count
            self._count_started = count_new
        return expired


class ExpiringTimer(ExecutionTimer):
    """A timer that expires after a certain amount of time."""

    def __init__(self, interval, started_at=None, expire_once=False):
        """Initialize instance with a time interval.

        Parameters
        ----------
        interval : float
            How much time in seconds should this timer take to expire.
        started_at : float, optional
            When this timer was started. If not given, the timer
            is started right now. Default is None.
        expire_once : bool, optional
            Whether this timer should certainly expire the first time
            `has_expired` is called. After that it will expire after
            `interval` seconds since the last call to `has_expired`.

        Returns
        -------
        None.
        """
        self._interval = interval
        self._interval_started = 0  # super() sets it via self.restart
        super().__init__(started_at=started_at)
        if expire_once:
            # Subtracting a bit more ensures we expire right away.
            # One interval is not enough, as we check ">", not ">="
            self._interval_started -= 1.1*self.interval

    @property
    def interval(self):
        """Return the interval of this timer."""
        return self._interval

    def has_expired(self):
        """Return whether this interval has gone by and how long in total.

        If the timer has expired its internal clock for the
        interval is also reset to right now.

        Returns
        -------
        has_expired : bool
            Whether enough time has passed for this timer to expire.
        """
        now = self.now()
        how_long_interval = now - self._interval_started
        expired = how_long_interval > self.interval
        if expired:
            self._interval_started = now
        return expired

    def restart(self, started_at=None):
        """Start this timer again."""
        super().restart(started_at=started_at)
        self._interval_started = self.started_at
def _elapsed_time_as_str(interval):
    """Return an elapsed-time string from an interval in seconds.

    Parameters
    ----------
    interval : float
        Number of seconds passed.

    Returns
    -------
    interval_str : str
        A formatted version of interval. The current formats are:
            - "H:MM hours"     if interval is longer than 1 h
            - "M:SS minutes"   if interval is less more than 1 min
            - "S.mm seconds"   if interval is shorter than 1 min
    """
    delta = timedelta(seconds=interval)
    hours = delta.seconds // _ONE_HOUR
    minutes = (delta.seconds % _ONE_HOUR) // _ONE_MINUTE
    seconds = delta.seconds % _ONE_MINUTE
    milliseconds = delta.microseconds // 1000
    if hours:
        return f'{hours}:{minutes:02d} hours'
    if minutes:
        return f'{minutes}:{seconds:02d} minutes'
    return f'{seconds}.{milliseconds//10:02d} seconds'
