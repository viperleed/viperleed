"""Tests for module tier_marks of tests.harness."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

import pytest

from .conftest import MODULE
from .conftest import XDIST_ARGS

_FAILING_FILE = """
import pytest

@pytest.mark.parametrize('case', ['a', 'b-c'])
def test_cases(case):
    assert case != 'a'

def test_other():
    pass
"""
# The tables are the `slow` and `release` arguments of the run fixture
_STALE = [
    pytest.param(
        {'slow': {MODULE + 'test_renamed': ('a',)}},
        f'*{MODULE}test_renamed (no such test)',
        '(no such case)',
        id='slow-test',
        ),
    pytest.param(
        {'slow': {MODULE + 'test_cases': ('a', 'renamed')}},
        f'*{MODULE}test_cases[[]renamed] (no such case)',
        '[a]',
        id='slow-case',
        ),
    pytest.param(  # 'case' is the argument of test_cases
        {'release': {'case': ('a', 'renamed')}},
        '*case[[]renamed] (no such case)',
        '[a]',
        id='release-case',
        ),
    ]
_NO_ERROR = [
    pytest.param(
        {MODULE + 'test_cases': ('a', 'b-c')},
        (),
        4,
        id='valid',
        ),
    pytest.param(
        {MODULE + 'test_renamed': ('a',)},
        ('test_unrelated.py',),
        1,
        id='module-not-collected',
        ),
    pytest.param(
        {MODULE + 'test_cases': ('a',)},
        (MODULE + 'test_other',),
        1,
        id='selected-by-node-id',
        ),
    pytest.param(  # Still counted as collected
        {MODULE + 'test_cases': ('a',)},
        ('-k', 'test_other'),
        1,
        id='deselected-by-keyword',
        ),
    ]


class TestStaleSlowCases:
    """Tests for the error on tier-table entries that match nothing."""

    @pytest.mark.parametrize('args', [('--lf',), ('--lf', '.')])
    def test_last_failed(self, run, args):
        """Check that --lf, which drops passed items, is no error."""
        slow = {MODULE + 'test_cases': ('a', 'b-c')}
        files = {'test_sample': _FAILING_FILE}
        first = run(slow=slow, files=files, cache=True)
        first.assert_outcomes(passed=2, failed=1)
        result = run(*args, slow=slow, files=files, cache=True)
        result.assert_outcomes(failed=1)
        assert result.ret is pytest.ExitCode.TESTS_FAILED

    @pytest.mark.parametrize(('entries', 'args', 'passed'), _NO_ERROR)
    def test_no_error(self, run, entries, args, passed):
        """Check that valid entries and subset runs are no error."""
        result = run(*args, slow=entries)
        result.assert_outcomes(passed=passed)

    @pytest.mark.parametrize(('tables', 'line', 'absent'), _STALE)
    @pytest.mark.parametrize(
        'args',
        [
            pytest.param((), id='serial'),
            pytest.param(XDIST_ARGS, marks=pytest.mark.slow, id='xdist'),
            ],
        )
    def test_reported(self, run, tables, line, absent, args):
        """Check that a stale entry is a usage error naming it."""
        (table,) = tables
        result = run(*args, **tables)
        assert result.ret is pytest.ExitCode.USAGE_ERROR
        result.stderr.fnmatch_lines(
            [
                f'*Stale entries in {table.upper()}_CASES*Update the module.',
                line,
                ],
            )
        assert absent not in result.stderr.str()

    def test_unused_fixture_is_no_error(self, run):
        """Check that a fixture that no test uses is not stale."""
        result = run(release={'renamed': ('a',)})
        result.assert_outcomes(passed=4)
