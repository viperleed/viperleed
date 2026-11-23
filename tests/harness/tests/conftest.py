"""Test configuration for tests/harness/tests.

Defines the fixtures and constants shared by the tests of the tier
plugins (tier_marks, check_slow and release), which run sample tests
with pytester.

Fixtures
--------
release_cases
    An empty dict in place of tiers.RELEASE_CASES.
run (factory)
    Run sample tests with the tier options.
slow_cases
    An empty dict in place of tiers.SLOW_CASES.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

import os
from pathlib import Path

import pytest

from ... import conftest as root_conftest
from .. import tiers

MODULE = 'test_sample.py::'  # Node-id prefix of the sample tests
XDIST_ARGS = ('-n', '2')

_INI = """
[pytest]
markers =
    slow
    tensorleed
    release
"""
_OTHER_FILE = """
def test_unrelated():
    pass
"""
_REPO_ROOT = Path(__file__).resolve().parents[3]
_SLOW_PLUGIN = 'slow_cases_plugin'
_TEST_FILE = """
import pytest

@pytest.mark.parametrize('case', ['a', 'b-c'])
def test_cases(case):
    pass

def test_other():
    pass
"""
_SAMPLE_FILES = {'test_sample': _TEST_FILE, 'test_unrelated': _OTHER_FILE}


@pytest.fixture(name='release_cases')
def fixture_release_cases(monkeypatch):
    """Replace tiers.RELEASE_CASES with a dict that tests can fill."""
    cases = {}
    monkeypatch.setattr(tiers, 'RELEASE_CASES', cases)
    return cases


@pytest.fixture(name='run')
def fixture_run(pytester, monkeypatch, slow_cases, release_cases):
    """Return a function that runs sample tests with the tier options.

    The function takes the command-line arguments, the `slow` cases
    to put in tiers.SLOW_CASES, the `release` cases to put in
    tiers.RELEASE_CASES, the `files` to write (default: the
    sample tests), and whether to keep the `cache` (default: no). It
    runs in the current process, unless `own_process` is True (default:
    no) or the arguments request pytest-xdist. Use `own_process` for
    runs with pytest-cov, which must not interfere with the coverage
    of the outer run. A separate process gets `slow` handed over with
    a plugin.
    """
    def _run(*args, slow=None, release=None, files=None, cache=False,
             own_process=False):
        slow_cases.update(slow or {})
        release_cases.update(release or {})
        pytester.makeini(_INI)
        pytester.makepyfile(**(files or _SAMPLE_FILES))
        options = list(args)
        if '-m' not in args:
            options.extend(('-m', ''))
        if not cache:
            options.extend(('-p', 'no:cacheprovider'))
        if XDIST_ARGS[0] not in args and not own_process:
            return pytester.runpytest_inprocess(
                *options,
                plugins=[root_conftest],
                )
        pytester.makepyfile(**{
            _SLOW_PLUGIN: f'from tests.harness import tiers\n'
                          f'tiers.SLOW_CASES.update({slow_cases!r})\n'
                          f'tiers.RELEASE_CASES.update({release_cases!r})\n',
            })
        monkeypatch.setenv(
            'PYTHONPATH',
            os.pathsep.join((str(_REPO_ROOT), str(pytester.path))),
            )
        return pytester.runpytest_subprocess(
            *options,
            '-p', _SLOW_PLUGIN,
            '-p', root_conftest.__name__,
            )
    return _run


@pytest.fixture(name='slow_cases')
def fixture_slow_cases(monkeypatch):
    """Replace tiers.SLOW_CASES with a dict that tests can fill."""
    cases = {}
    monkeypatch.setattr(tiers, 'SLOW_CASES', cases)
    return cases
