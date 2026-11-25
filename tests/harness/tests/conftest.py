"""Test configuration for tests/harness/tests.

Most tests of the harness run sample tests with pytester, in the
current process or in a subprocess. This module defines what these
runs share, and the fixtures and constants of the tests of the tier
plugins (tier_marks, check_slow and release).

Fixtures
--------
release_cases
    An empty dict in place of tiers.RELEASE_CASES.
run (factory)
    Run sample tests with the tier options.
slow_cases
    An empty dict in place of tiers.SLOW_CASES.
subprocess_path (autouse)
    Make package tests importable in pytest subprocesses.

Functions
---------
run_pytest
    Run pytest on the files of pytester, with a plugin of the tests.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

import os

import pytest
from pytest_cases import fixture

from ... import conftest as root_conftest
from .. import qt
from .. import tiers

# The markers of the test suite, for the ini file of sample tests
INI = """
[pytest]
markers =
    slow
    tensorleed
    release
"""
MODULE = 'test_sample.py::'  # Node-id prefix of the sample tests
# Options of the nested pytest runs (see run_pytest). The cache is
# disabled, unless a test needs it (e.g., for --last-failed). The qt
# plugin of the harness is always loaded: otherwise, pytest-qt (if
# installed) is loaded automatically, and aborts the run if PyQt5 is
# not installed.
PYTESTER_NO_CACHE = ('-p', 'no:cacheprovider')
PYTESTER_PLUGINS = ('-p', qt.__name__)
XDIST_ARGS = ('-n', '2')
# Extra arguments of runs without and with pytest-xdist, keyed by id,
# for parametrize. The runs with pytest-xdist are slow.
SERIAL_AND_XDIST = {
    'serial': (),
    'xdist': pytest.param(XDIST_ARGS, marks=pytest.mark.slow),
    }

_OTHER_FILE = """
def test_unrelated():
    pass
"""
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


@fixture(name='release_cases')
def fixture_release_cases(monkeypatch):
    """Replace tiers.RELEASE_CASES with a dict that tests can fill."""
    cases = {}
    monkeypatch.setattr(tiers, 'RELEASE_CASES', cases)
    return cases


@fixture(name='run')
def fixture_run(pytester, slow_cases, release_cases):
    """Return a function that runs sample tests with the tier options.

    The function takes the command-line arguments, the `slow` cases
    to put in tiers.SLOW_CASES, the `release` cases to put in
    tiers.RELEASE_CASES, the `files` to write (default: the
    sample tests), and the `cache` and `own_process` arguments of
    run_pytest. A separate process gets `slow` and `release` handed
    over with a plugin.
    """
    def _run(*args, slow=None, release=None, files=None, cache=False,
             own_process=False):
        slow_cases.update(slow or {})
        release_cases.update(release or {})
        pytester.makeini(INI)
        pytester.makepyfile(**(files or _SAMPLE_FILES))
        options = list(args)
        if '-m' not in args:
            options.extend(('-m', ''))
        own_process = _in_own_process(args, own_process)
        if own_process:
            pytester.makepyfile(**{
                _SLOW_PLUGIN:
                    'from tests.harness import tiers\n'
                    f'tiers.SLOW_CASES.update({slow_cases!r})\n'
                    f'tiers.RELEASE_CASES.update({release_cases!r})\n',
                })
            options.extend(('-p', _SLOW_PLUGIN))
        return run_pytest(pytester, *options, cache=cache,
                          own_process=own_process)
    return _run


@fixture(name='slow_cases')
def fixture_slow_cases(monkeypatch):
    """Replace tiers.SLOW_CASES with a dict that tests can fill."""
    cases = {}
    monkeypatch.setattr(tiers, 'SLOW_CASES', cases)
    return cases


@fixture(name='subprocess_path', autouse=True)
def fixture_subprocess_path(monkeypatch, repo_root):
    """Make package tests importable in pytest subprocesses.

    The root of the repository goes first in PYTHONPATH, so that the
    subprocesses import the plugins of this copy of the tests. The
    folder of the pytester files is added by pytester itself.
    """
    python_path = [str(repo_root)]
    if os.environ.get('PYTHONPATH'):
        python_path.append(os.environ['PYTHONPATH'])
    monkeypatch.setenv('PYTHONPATH', os.pathsep.join(python_path))


def run_pytest(pytester, *args, plugin=root_conftest, cache=False,
               own_process=False):
    """Run pytest on the files of `pytester` with a plugin of the tests.

    Parameters
    ----------
    pytester : pytest.Pytester
        The pytester whose files are tested.
    *args : str
        Command-line arguments.
    plugin : module, optional
        The plugin to load. Default is tests/conftest.py, which loads
        all the plugins of the test suite. The qt plugin of the harness
        is loaded in any case (see PYTESTER_PLUGINS).
    cache : bool, optional
        Whether to keep the cache (e.g., for --lf). Default is False.
    own_process : bool, optional
        Whether to run pytest in a separate process. Default is False.
        A separate process is used anyway if `args` request pytest-
        xdist. Use it for runs with pytest-cov, which must not
        interfere with the coverage of the outer run.

    Returns
    -------
    result : pytest.RunResult
        The outcome of the run.
    """
    options = [*args, *PYTESTER_PLUGINS]
    if not cache:
        options.extend(PYTESTER_NO_CACHE)
    if _in_own_process(args, own_process):
        return pytester.runpytest_subprocess(*options, '-p', plugin.__name__)
    return pytester.runpytest_inprocess(*options, plugins=[plugin])


def _in_own_process(args, own_process):
    """Return whether pytest with `args` runs in a separate process."""
    return own_process or XDIST_ARGS[0] in args
