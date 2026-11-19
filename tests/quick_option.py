"""Pytest plugin for the --quick option of the test suite."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-27'
__license__ = 'GPLv3+'

from . import markexpr

QUICK_OPTION = '--quick'


def pytest_addoption(parser):
    """Add the --quick option."""
    parser.addoption(
        QUICK_OPTION,
        action='store_true',
        help=('Quick tier: skip the tests marked "slow", in addition '
              'to what the -m expression (if any) selects.'),
        )


def pytest_configure(config):
    """Add 'not slow' to the marker expression if --quick was given.

    This must be done in each pytest-xdist worker: the options that
    xdist hands to a worker are those parsed from the command line,
    not the ones modified by the controller's pytest_configure.
    """
    if not config.option.quick:
        return
    markexpr.exclude_marker(config, 'slow')
