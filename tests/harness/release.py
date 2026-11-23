"""Pytest plugin for the --release option of the test suite."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-28'
__license__ = 'GPLv3+'

import pytest

from . import markexpr

_DESELECTED = pytest.StashKey()
RELEASE_OPTION = '--release'


def pytest_addoption(parser):
    """Add the --release option."""
    parser.addoption(
        RELEASE_OPTION,
        action='store_true',
        help=('Release tier: also run the tests marked "release", which '
              'are deselected otherwise, even with -m tensorleed.'),
        )


def pytest_configure(config):
    """Add 'not release' to the marker expression unless --release.

    This must be done in each pytest-xdist worker: the options that
    xdist hands to a worker are those parsed from the command line,
    not the ones modified by the controller's pytest_configure.
    """
    if config.option.release:
        return
    markexpr.exclude_marker(config, 'release')


def pytest_deselected(items):
    """Count the deselected release cases.

    Nothing is counted with -k, as it also deselects items for
    other reasons than the marker expression.
    """
    if not items:
        return
    config = items[0].config
    if config.option.release or config.option.keyword:
        return
    release = sum(1 for item in items if item.get_closest_marker('release'))
    stash = config.stash
    stash[_DESELECTED] = stash.get(_DESELECTED, 0) + release


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Say why no test ran if release cases were deselected."""
    count = config.stash.get(_DESELECTED, 0)
    if not count or exitstatus != pytest.ExitCode.NO_TESTS_COLLECTED:
        return
    terminalreporter.write_line(
        f'{count} release case(s) deselected: they only run with '
        f'{RELEASE_OPTION} (and with a -m that allows tensorleed tests)'
        )
