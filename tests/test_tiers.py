"""Tests for the tier data and options of the test suite."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

from importlib import import_module
import os
from pathlib import Path

import pytest

from . import conftest
from . import tiers

_TEST_FILE = """
import pytest

@pytest.mark.parametrize('case', ['a', 'b-c'])
def test_cases(case):
    pass

def test_other():
    pass
"""
_FAILING_FILE = """
import pytest

@pytest.mark.parametrize('case', ['a', 'b-c'])
def test_cases(case):
    assert case != 'a'

def test_other():
    pass
"""
_OTHER_FILE = """
def test_unrelated():
    pass
"""
_SAMPLE_FILES = {'test_sample': _TEST_FILE, 'test_unrelated': _OTHER_FILE}
_INI = """
[pytest]
markers =
    slow
    tensorleed
"""
_MODULE = 'test_sample.py::'
_REPO_ROOT = Path(__file__).resolve().parent.parent
_SLOW_PLUGIN = 'slow_cases_plugin'
_XDIST_ARGS = ('-n', '2')


@pytest.fixture(name='slow_cases')
def fixture_slow_cases(monkeypatch):
    """Replace tiers.SLOW_CASES with a dict that tests can fill."""
    cases = {}
    monkeypatch.setattr(tiers, 'SLOW_CASES', cases)
    return cases


@pytest.fixture(name='run')
def fixture_run(pytester, monkeypatch, slow_cases):
    """Return a function that runs sample tests with the tier options.

    The function takes the command-line arguments, the `slow` cases
    to put in tiers.SLOW_CASES, the `files` to write (default: the
    sample tests), and whether to keep the `cache` (default: no). It
    runs in a subprocess if the arguments request pytest-xdist, and
    fails (not skips) if it is missing. There, `slow` is handed over
    with a plugin.
    """
    def _run(*args, slow=None, files=None, cache=False):
        slow_cases.update(slow or {})
        pytester.makeini(_INI)
        pytester.makepyfile(**(files or _SAMPLE_FILES))
        options = [
            *args,
            '-m', '',
            ]
        if not cache:
            options.extend(('-p', 'no:cacheprovider'))
        if _XDIST_ARGS[0] not in args:
            return pytester.runpytest_inprocess(
                *options,
                plugins=[conftest],
                )
        import_module('xdist')  # Fail if missing
        pytester.makepyfile(**{
            _SLOW_PLUGIN: f'from tests import tiers\n'
                          f'tiers.SLOW_CASES.update({slow_cases!r})\n',
            })
        monkeypatch.setenv(
            'PYTHONPATH',
            os.pathsep.join((str(_REPO_ROOT), str(pytester.path))),
            )
        return pytester.runpytest_subprocess(
            *options,
            '-p', _SLOW_PLUGIN,
            '-p', conftest.__name__,
            )
    return _run


_STALE = [
    pytest.param(
        {_MODULE + 'test_renamed': ('a',)},
        f'*{_MODULE}test_renamed (no such test)',
        '(no such case)',
        id='test',
        ),
    pytest.param(
        {_MODULE + 'test_cases': ('a', 'renamed')},
        f'*{_MODULE}test_cases[[]renamed] (no such case)',
        '[a]',
        id='case',
        ),
    ]
_NO_ERROR = [
    pytest.param(
        {_MODULE + 'test_cases': ('a', 'b-c')},
        (),
        4,
        id='valid',
        ),
    pytest.param(
        {_MODULE + 'test_renamed': ('a',)},
        ('test_unrelated.py',),
        1,
        id='module-not-collected',
        ),
    pytest.param(
        {_MODULE + 'test_cases': ('a',)},
        (_MODULE + 'test_other',),
        1,
        id='selected-by-node-id',
        ),
    pytest.param(  # Still counted as collected
        {_MODULE + 'test_cases': ('a',)},
        ('-k', 'test_other'),
        1,
        id='deselected-by-keyword',
        ),
    ]


class TestStaleSlowCases:
    """Tests for the error on SLOW_CASES entries that match nothing."""

    @pytest.mark.parametrize(('entries', 'line', 'absent'), _STALE)
    @pytest.mark.parametrize(
        'args',
        [
            pytest.param((), id='serial'),
            pytest.param(_XDIST_ARGS, marks=pytest.mark.slow, id='xdist'),
            ],
        )
    def test_reported(self, run, entries, line, absent, args):
        """Check that a stale entry is a usage error naming it."""
        result = run(*args, slow=entries)
        assert result.ret is pytest.ExitCode.USAGE_ERROR
        result.stderr.fnmatch_lines(
            [
                '*Stale entries in SLOW_CASES*Update the module.',
                line,
                ],
            )
        assert absent not in result.stderr.str()

    @pytest.mark.parametrize(('entries', 'args', 'passed'), _NO_ERROR)
    def test_no_error(self, run, entries, args, passed):
        """Check that valid entries and subset runs are no error."""
        result = run(*args, slow=entries)
        result.assert_outcomes(passed=passed)

    @pytest.mark.parametrize('args', [('--lf',), ('--lf', '.')])
    def test_last_failed(self, run, args):
        """Check that --lf, which drops passed items, is no error."""
        slow = {_MODULE + 'test_cases': ('a', 'b-c')}
        files = {'test_sample': _FAILING_FILE}
        first = run(slow=slow, files=files, cache=True)
        first.assert_outcomes(passed=2, failed=1)
        result = run(*args, slow=slow, files=files, cache=True)
        result.assert_outcomes(failed=1)
        assert result.ret is pytest.ExitCode.TESTS_FAILED
