"""Tests for module optional_dependencies of viperleed.tests."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-21'
__license__ = 'GPLv3+'

import os
from pathlib import Path

import pytest
from pytest_cases import parametrize

from . import optional_dependencies

_ALLOW = optional_dependencies.ALLOW_MISSING_OPTION
_REQUIRE = optional_dependencies.REQUIRE_ALL_OPTION

_TEST_FILE = """
import pytest

def test_fixture(missing_dependency):
    missing_dependency('gfortran', 'some detail')

def test_unaffected():
    pass"""


def _run(pytester, *args, source=_TEST_FILE):
    """Run pytest on a test file with the plugin."""
    pytester.makepyfile(source)
    return pytester.runpytest_inprocess(
        *args,
        '-p', 'no:cacheprovider',
        plugins=[optional_dependencies],
        )


def _run_xdist(pytester, mocker, *args):
    """Run pytest with pytest-xdist on the test file with the plugin."""
    root = Path(__file__).resolve().parent.parent
    mocker.patch.dict(os.environ, {'PYTHONPATH': str(root)})
    pytester.makepyfile(_TEST_FILE)
    return pytester.runpytest_subprocess(
        *args,
        '-n', '2',
        '-p', 'no:cacheprovider',
        '-p', optional_dependencies.__name__,
        )


_UNREGISTERED = 'not-a-dependency'
_UNREGISTERED_FILE = f"""
def test_unregistered(missing_dependency):
    missing_dependency({_UNREGISTERED!r})
"""


@parametrize(option=(_REQUIRE, _ALLOW))
def test_option_in_help(pytester, option):
    """Check that the options are documented in the help."""
    result = _run(pytester, '--help')
    result.stdout.fnmatch_lines([f'*{option}*'])


def test_help_lists_known(pytester):
    """Check that --require-all lists the known dependencies."""
    result = _run(pytester, '--help')
    known = '/'.join(optional_dependencies.KNOWN_DEPENDENCIES)
    result.stdout.fnmatch_lines([f'*({known})*'])


@parametrize(args=((), (_REQUIRE,)))
def test_unregistered_name(pytester, args):
    """Check that a test fails on an unregistered dependency name."""
    result = _run(pytester, *args, source=_UNREGISTERED_FILE)
    result.assert_outcomes(failed=1)
    result.stdout.fnmatch_lines([
        f'*{_UNREGISTERED!r} is not one of the KNOWN_DEPENDENCIES*',
        ])


def test_skip_without_option(pytester):
    """Check that missing dependencies skip, naming the dependency."""
    result = _run(pytester, '-rs')
    result.assert_outcomes(passed=1, skipped=1)
    result.stdout.fnmatch_lines([
        '*optional dependency gfortran not available: some*',
        ])


def test_fail_with_option(pytester):
    """Check that --require-all turns skips into failures."""
    result = _run(pytester, optional_dependencies.REQUIRE_ALL_OPTION)
    result.assert_outcomes(passed=1, failed=1)
    result.stdout.fnmatch_lines([
        '*optional dependency gfortran not available: some*',
        ])


_ALLOWED_DEP = (
    (f'{_ALLOW}=gfortran',),
    (f'{_ALLOW}=mpifort,gfortran',),
    (f'{_ALLOW}= mpifort , gfortran ',),
    (f'{_ALLOW}=mpifort', f'{_ALLOW}=gfortran'),
    )


@parametrize(allow=_ALLOWED_DEP)
def test_allow_missing(pytester, allow):
    """Check that --allow-missing exempts the named dependencies."""
    result = _run(pytester, _REQUIRE, *allow, '-rs')
    result.assert_outcomes(passed=1, skipped=1)
    result.stdout.fnmatch_lines([
        f'*gfortran not available: some*(allowed by {_ALLOW})*',
        ])


@parametrize(allow=(f'{_ALLOW}=mpifort', f'{_ALLOW}=mpifort,'))
def test_allow_missing_other(pytester, allow):
    """Check that other missing dependencies still fail."""
    result = _run(pytester, _REQUIRE, allow)
    result.assert_outcomes(passed=1, failed=1)
    result.stdout.fnmatch_lines([f'*(required by {_REQUIRE})*'])


def test_allow_missing_alone(pytester):
    """Check that --allow-missing alone does not change the skips."""
    result = _run(pytester, f'{_ALLOW}=gfortran', '-rs')
    result.assert_outcomes(passed=1, skipped=1)
    result.stdout.no_fnmatch_line('*allowed by*')


@pytest.mark.slow  # About 1 s (pytester, subprocess)
@parametrize('allow,outcome', (
    ((), {'failed': 1}),
    ((f'{_ALLOW}=gfortran',), {'skipped': 1}),
    ))
def test_allow_missing_xdist(pytester, mocker, allow, outcome):
    """Check --require-all and --allow-missing with pytest-xdist."""
    result = _run_xdist(pytester, mocker, _REQUIRE, *allow)
    result.assert_outcomes(passed=1, **outcome)


_UNKNOWN = (
    (f'{_ALLOW}=gfortan',),
    (f'{_ALLOW}=gfortran,fake-dep',),
    (f'{_ALLOW}=mpifort', f'{_ALLOW}=gfortan'),
    )


@parametrize(allow=_UNKNOWN)
def test_allow_missing_unknown(pytester, allow):
    """Check that unknown names for --allow-missing are an error."""
    result = _run(pytester, _REQUIRE, *allow)
    assert result.ret == pytest.ExitCode.USAGE_ERROR
    result.stderr.fnmatch_lines([
        (f'*{_ALLOW}: unknown optional dependencies *. '
         'Valid names are gfortran, mpifort.'),
        ])


@pytest.mark.slow  # About 1 s (pytester, subprocess)
def test_allow_missing_unknown_xdist(pytester, mocker):
    """Check that unknown names are an error with pytest-xdist."""
    result = _run_xdist(pytester, mocker, _REQUIRE, f'{_ALLOW}=gfortan')
    assert result.ret == pytest.ExitCode.USAGE_ERROR
    result.stderr.fnmatch_lines(['*unknown optional dependencies gfortan.*'])
