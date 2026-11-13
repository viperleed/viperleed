"""Pytest plugin for the --check-slow option of the test suite."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-27'
__license__ = 'GPLv3+'

import pytest

from . import tiers

_CALL_PHASE = 'call'
_SLOW_MARKER = 'slow'
CHECK_SLOW_OPTION = '--check-slow'
_PYTEST_COV_PLUGIN = '_cov'  # Registered by pytest-cov only if --cov


def _upper(threshold):
    """Return the duration from which an unmarked test is listed."""
    return (1 + tiers.SLOW_BAND) * threshold


def _lower(threshold):
    """Return the duration below which a marked test is a candidate."""
    return (1 - tiers.SLOW_BAND) * threshold


def pytest_addoption(parser):
    """Add the --check-slow option."""
    parser.addoption(
        CHECK_SLOW_OPTION,
        action='store_true',
        help=('At the end, list the tests not marked "slow" that took '
              f'{_upper(tiers.SLOW_THRESHOLD):.2f} s or more, and the '
              f'"slow" ones that took less than '
              f'{_lower(tiers.SLOW_THRESHOLD):.2f} s (call phase), as '
              'candidates to be marked or unmarked in tests/tiers.py. '
              'Only reports; also with pytest-xdist.'),
        )


def _coverage_is_on(config):
    """Return whether pytest-cov measures coverage in this run."""
    return (config.pluginmanager.hasplugin(_PYTEST_COV_PLUGIN)
            and not config.getoption('no_cov', default=False))


def pytest_configure(config):
    """Start recording durations if --check-slow was given.

    Raises
    ------
    pytest.UsageError
        If --check-slow is given while coverage is on: it makes tests
        2-3 times slower, so that the durations are not meaningful.
    """
    # With pytest-xdist, only the controller records
    if (getattr(config, 'workerinput', None) is not None
            or not config.option.check_slow):
        return
    if _coverage_is_on(config):
        message = (
            f'{CHECK_SLOW_OPTION} needs --no-cov: coverage makes the '
            'tests 2-3 times slower, and the durations meaningless'
            )
        raise pytest.UsageError(message)
    config.pluginmanager.register(
        _DurationRecorder(),
        _DurationRecorder.name,
        )


class _DurationRecorder:
    """Plugin that records the call durations of the tests that run."""

    name = 'check-slow-recorder'

    def __init__(self):
        """Initialize an instance with no durations."""
        self.durations = {}  # nodeid: seconds
        self.slow_ids = set()

    def pytest_runtest_logreport(self, report):
        """Record the call duration of a report."""
        # With pytest-xdist, the controller receives all the reports.
        # Subtests report under the id of their test: keep the longest
        if report.when != _CALL_PHASE or report.skipped:
            return
        self.durations[report.nodeid] = max(
            self.durations.get(report.nodeid, 0),
            report.duration,
            )
        if _SLOW_MARKER in report.keywords:
            self.slow_ids.add(report.nodeid)

    def candidates(self, threshold):
        """Return the tests whose marking may need to change.

        Parameters
        ----------
        threshold : float
            Duration (seconds) from which a test counts as slow.
            Tests within tiers.SLOW_BAND of it keep their marking.

        Returns
        -------
        to_mark : list
            (duration, node id) of tests that are not marked slow,
            but took (1 + SLOW_BAND) * `threshold` or more, longest
            first.
        to_unmark : list
            (duration, node id) of tests marked slow that took less
            than (1 - SLOW_BAND) * `threshold`, shortest first.
        """
        upper, lower = _upper(threshold), _lower(threshold)
        to_mark = sorted(
            ((time, nodeid) for nodeid, time in self.durations.items()
             if nodeid not in self.slow_ids and time >= upper),
            reverse=True,
            )
        to_unmark = sorted(
            (time, nodeid) for nodeid, time in self.durations.items()
            if nodeid in self.slow_ids and time < lower
            )
        return to_mark, to_unmark


def _write_slow_candidates(terminalreporter, config):
    """Write the lists of tests whose marking may need to change."""
    recorder = config.pluginmanager.get_plugin(_DurationRecorder.name)
    threshold = tiers.SLOW_THRESHOLD
    to_mark, to_unmark = recorder.candidates(threshold)
    terminalreporter.section('check-slow (report only)')
    terminalreporter.write_line(
        f'Call phase, threshold {threshold} s, band '
        f'{tiers.SLOW_BAND:.0%}: tests in between keep their marking. '
        'Edit tests/tiers.py (SLOW_CASES) or the marks; timings are '
        'noisy: re-run a candidate before changing it.'
        )
    to_mark_title = (
        'Not marked "slow", but took the upper limit or more '
        f'({_upper(threshold):.2f} s):'
        )
    to_unmark_title = (
        'Marked "slow", but took less than the lower limit '
        f'({_lower(threshold):.2f} s):'
        )
    for title, tests in (
            (to_mark_title, to_mark),
            (to_unmark_title, to_unmark),
            ):
        terminalreporter.write_line(f'{title} {len(tests)}')
        for time, nodeid in tests:
            terminalreporter.write_line(f'  {time:7.2f} s  {nodeid}')


def pytest_terminal_summary(terminalreporter, config):
    """Report the tests whose marking may need to change."""
    if config.option.check_slow:
        _write_slow_candidates(terminalreporter, config)
