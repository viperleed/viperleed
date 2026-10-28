"""Tests for the tier data and options of the test suite."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

import os
from pathlib import Path

import pytest

from . import check_slow
from . import conftest
from . import release_option
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
    release
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
        if _XDIST_ARGS[0] not in args and not own_process:
            return pytester.runpytest_inprocess(
                *options,
                plugins=[conftest],
                )
        pytester.makepyfile(**{
            _SLOW_PLUGIN: f'from tests import tiers\n'
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


_SLOW_FILE = """
import time

import pytest

def test_fast():
    pass

def test_unmarked_slow():
    time.sleep(0.7)

@pytest.mark.slow
def test_marked_fast():
    pass

@pytest.mark.slow
def test_marked_slow():
    time.sleep(0.7)
"""
_SLOW_FILES = {'test_durations': _SLOW_FILE}
_DURATIONS = 'test_durations.py::'
_BELOW = 'took less than the threshold'


@pytest.mark.slow  # About 1 s each (pytester)
class TestCheckSlow:
    """Tests for the --check-slow option."""

    @staticmethod
    def _check_report(result):
        """Check the lists of --check-slow, and that nothing failed."""
        result.assert_outcomes(passed=4)  # Report only
        assert result.ret is pytest.ExitCode.OK
        out = result.stdout.str()
        head, _, rest = out.partition(_BELOW)
        assert _DURATIONS + 'test_unmarked_slow' in head
        assert _DURATIONS + 'test_marked_fast' in rest
        for name in ('test_fast', 'test_marked_slow'):
            assert f'{_DURATIONS}{name}\n' not in out

    def test_candidates(self, run):
        """Check that tests to mark and to unmark are listed."""
        result = run(check_slow.CHECK_SLOW_OPTION, files=_SLOW_FILES)
        self._check_report(result)

    def test_needs_no_cov(self, run):
        """Check that --check-slow with coverage on is a usage error."""
        result = run(
            check_slow.CHECK_SLOW_OPTION,
            '--cov=.',
            '--cov-report=',
            files=_SLOW_FILES,
            own_process=True,
            )
        assert result.ret is pytest.ExitCode.USAGE_ERROR
        result.stderr.fnmatch_lines(['*--check-slow needs --no-cov*'])

    def test_no_cov(self, run):
        """Check that --no-cov overrides --cov for --check-slow."""
        result = run(
            check_slow.CHECK_SLOW_OPTION,
            '--cov=.',
            '--cov-report=',
            '--no-cov',
            files=_SLOW_FILES,
            own_process=True,
            )
        self._check_report(result)

    def test_without_option(self, run):
        """Check that nothing is reported without the option."""
        result = run(files=_SLOW_FILES)
        result.assert_outcomes(passed=4)
        assert result.ret is pytest.ExitCode.OK
        assert check_slow.CHECK_SLOW_OPTION.lstrip('-') not in (
            result.stdout.str()
            )

    def test_xdist(self, run):
        """Check that --check-slow also works with pytest-xdist."""
        result = run(
            check_slow.CHECK_SLOW_OPTION,
            *_XDIST_ARGS,
            files=_SLOW_FILES,
            )
        self._check_report(result)


_RELEASE_FILE = """
import pytest

@pytest.mark.tensorleed
@pytest.mark.parametrize('case', ['a', 'b-c'])
def test_tl(case):
    pass

def test_plain():
    pass
"""
_RELEASE_FILES = {'test_sample': _RELEASE_FILE}
_RELEASE = {_MODULE + 'test_tl': ('b-c',)}
_RELEASE_OPTION = release_option.RELEASE_OPTION
# Numbers of tests that run, out of test_tl[a], test_tl[b-c], test_plain
_RELEASE_SELECTION = [
    pytest.param((), 2, id='default'),
    pytest.param(('-m', 'tensorleed'), 1, id='tensorleed'),
    pytest.param(('-m', 'not tensorleed'), 1, id='not-tensorleed'),
    pytest.param((_RELEASE_OPTION,), 3, id='release'),
    pytest.param(
        (_RELEASE_OPTION, '-m', 'tensorleed'), 2, id='release-tensorleed',
        ),
    pytest.param((_RELEASE_OPTION, '--quick'), 3, id='release-quick'),
    ]


class TestReleaseOption:
    """Tests for the release marker and the --release option."""

    @pytest.mark.parametrize(('args', 'passed'), _RELEASE_SELECTION)
    @pytest.mark.parametrize(
        'xdist',
        [
            pytest.param((), id='serial'),
            pytest.param(_XDIST_ARGS, marks=pytest.mark.slow, id='xdist'),
            ],
        )
    def test_selection(self, run, args, passed, xdist):
        """Check that release cases only run with --release."""
        result = run(*args, *xdist, release=_RELEASE, files=_RELEASE_FILES)
        result.assert_outcomes(passed=passed)

    @pytest.mark.parametrize(
        ('args', 'present'),
        [
            pytest.param((), False, id='default'),
            pytest.param(('-m', 'tensorleed'), False, id='tensorleed'),
            pytest.param((_RELEASE_OPTION,), True, id='release'),
            pytest.param(
                ('-m', 'tensorleed', _RELEASE_OPTION),
                True,
                id='release-tensorleed',
                ),
            ],
        )
    def test_release_case_collected(self, run, args, present):
        """Check that only --release collects a release case."""
        result = run('--collect-only', '-q', *args, release=_RELEASE,
                     files=_RELEASE_FILES)
        assert ('test_tl[b-c]' in result.stdout.str()) is present
    @pytest.mark.parametrize(('entries', 'line', 'absent'), _STALE)
    @pytest.mark.parametrize(
        'args',
        [
            pytest.param((), id='serial'),
            pytest.param(_XDIST_ARGS, marks=pytest.mark.slow, id='xdist'),
            ],
        )
    def test_stale_entry(self, run, entries, line, absent, args):
        """Check that a stale RELEASE_CASES entry is a usage error."""
        result = run(*args, release=entries)
        assert result.ret is pytest.ExitCode.USAGE_ERROR
        result.stderr.fnmatch_lines(
            [
                '*Stale entries in RELEASE_CASES*Update the module.',
                line,
                ],
            )
        assert absent not in result.stderr.str()
