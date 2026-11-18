"""Tests for module optional_dependencies of viperleed.tests."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-21'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest
from pytest_cases import parametrize

from . import optional_dependencies

_ALLOW = optional_dependencies.ALLOW_MISSING_OPTION
_REQUIRE = optional_dependencies.REQUIRE_ALL_OPTION

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


def _run_xdist(pytester, monkeypatch, *args):
    """Run pytest with pytest-xdist on the test file with the plugin."""
    root = Path(__file__).resolve().parent.parent
    monkeypatch.setenv('PYTHONPATH', str(root))
    pytester.makepyfile(_TEST_FILE)
    return pytester.runpytest_subprocess(
        *args,
        '-n', '2',
        '-p', 'no:cacheprovider',
        '-p', optional_dependencies.__name__,
        )


@parametrize(option=(_REQUIRE, _ALLOW))
def test_option_in_help(pytester, option):
    """Check that the options are documented in the help."""
    result = _run(pytester, '--help')
    result.stdout.fnmatch_lines([f'*{option}*'])


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


_ALLOWED_DEP = (
    f'{_ALLOW}=fake-fixture-dep',
    f'{_ALLOW}=other-dep,fake-fixture-dep',
    f'{_ALLOW}= other-dep , fake-fixture-dep ',
    (f'{_ALLOW}=other-dep', f'{_ALLOW}=fake-fixture-dep'),
    )


@parametrize(allow=_ALLOWED_DEP)
def test_allow_missing(pytester, allow):
    """Check that --allow-missing exempts the named dependencies."""
    allow = (allow,) if isinstance(allow, str) else allow
    result = _run(pytester, _REQUIRE, *allow, '-rs')
    result.assert_outcomes(passed=1, skipped=1)
    result.stdout.fnmatch_lines([
        f'*fake-fixture-dep not available: some*(allowed by {_ALLOW})*',
        ])


@parametrize(allow=(f'{_ALLOW}=other-dep', f'{_ALLOW}=fake-fixture'))
def test_allow_missing_other(pytester, allow):
    """Check that other missing dependencies still fail."""
    result = _run(pytester, _REQUIRE, allow)
    result.assert_outcomes(passed=1, failed=1)
    result.stdout.fnmatch_lines([f'*(required by {_REQUIRE})*'])


def test_allow_missing_alone(pytester):
    """Check that --allow-missing alone does not change the skips."""
    result = _run(pytester, f'{_ALLOW}=fake-fixture-dep', '-rs')
    result.assert_outcomes(passed=1, skipped=1)
    result.stdout.no_fnmatch_line('*allowed by*')


@pytest.mark.slow  # About 1 s (pytester, subprocess)
@parametrize(allow=((), (f'{_ALLOW}=fake-fixture-dep',)))
def test_allow_missing_xdist(pytester, monkeypatch, allow):
    """Check --require-all and --allow-missing with pytest-xdist."""
    result = _run_xdist(pytester, monkeypatch, _REQUIRE, *allow)
    if allow:
        result.assert_outcomes(passed=1, skipped=1)
    else:
        result.assert_outcomes(passed=1, failed=1)
