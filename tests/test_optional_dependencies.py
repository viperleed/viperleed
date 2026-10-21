"""Tests for module optional_dependencies of viperleed.tests."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-21'
__license__ = 'GPLv3+'

from . import optional_dependencies

_TEST_FILE = """
import pytest

def test_fixture(missing_dependency):
    missing_dependency('fake-fixture-dep', 'some detail')

def test_unaffected():
    pass
"""


def _run(pytester, *args):
    """Run pytest on the test file with the plugin."""
    pytester.makepyfile(_TEST_FILE)
    return pytester.runpytest_inprocess(
        *args,
        '-p', 'no:cacheprovider',
        plugins=[optional_dependencies],
        )


def test_option_in_help(pytester):
    """Check that --require-all is documented in the help."""
    result = _run(pytester, '--help')
    result.stdout.fnmatch_lines(
        [f'*{optional_dependencies.REQUIRE_ALL_OPTION}*']
        )


def test_skip_without_option(pytester):
    """Check that missing dependencies skip, naming the dependency."""
    result = _run(pytester, '-rs')
    result.assert_outcomes(passed=1, skipped=1)
    result.stdout.fnmatch_lines([
        '*optional dependency fake-fixture-dep not available: some*',
        ])


def test_fail_with_option(pytester):
    """Check that --require-all turns skips into failures."""
    result = _run(pytester, optional_dependencies.REQUIRE_ALL_OPTION)
    result.assert_outcomes(passed=1, failed=1)
    result.stdout.fnmatch_lines([
        '*optional dependency fake-fixture-dep not available: some*',
        ])
