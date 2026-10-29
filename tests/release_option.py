"""Pytest plugin for the --release option of the test suite."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-28'
__license__ = 'GPLv3+'

from . import markexpr

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
