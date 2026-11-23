"""Tests for module tier_marks of tests.harness."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from .conftest import MODULE
from .conftest import SERIAL_AND_XDIST

_FAILING_FILE = """
import pytest

@pytest.mark.parametrize('case', ['a', 'b-c'])
def test_cases(case):
    assert case != 'a'

def test_other():
    pass
"""


class TestStaleSlowCases:
    """Tests for the error on tier-table entries that match nothing."""

    _no_error = {
        'valid': (
            {MODULE + 'test_cases': ('a', 'b-c')},
            (),
            4,
            ),
        'module-not-collected': (
            {MODULE + 'test_renamed': ('a',)},
            ('test_unrelated.py',),
            1,
            ),
        'selected-by-node-id': (
            {MODULE + 'test_cases': ('a',)},
            (MODULE + 'test_other',),
            1,
            ),
        'deselected-by-keyword': (  # Still counted as collected
            {MODULE + 'test_cases': ('a',)},
            ('-k', 'test_other'),
            1,
            ),
        }
    # The tables are the `slow` and `release` arguments of the run fixture
    _stale = {
        'slow-test': (
            {'slow': {MODULE + 'test_renamed': ('a',)}},
            f'*{MODULE}test_renamed (no such test)',
            '(no such case)',
            ),
        'slow-case': (
            {'slow': {MODULE + 'test_cases': ('a', 'renamed')}},
            f'*{MODULE}test_cases[[]renamed] (no such case)',
            '[a]',
            ),
        'release-case': (  # 'case' is the argument of test_cases
            {'release': {'case': ('a', 'renamed')}},
            '*case[[]renamed] (no such case)',
            '[a]',
            ),
        }

    # Not as a keyword argument, which would change the ids (args0...)
    @parametrize('args', [('--lf',), ('--lf', '.')])
    def test_last_failed(self, run, args):
        """Check that --lf, which drops passed items, is no error."""
        slow = {MODULE + 'test_cases': ('a', 'b-c')}
        files = {'test_sample': _FAILING_FILE}
        first = run(slow=slow, files=files, cache=True)
        first.assert_outcomes(passed=2, failed=1)
        result = run(*args, slow=slow, files=files, cache=True)
        result.assert_outcomes(failed=1)
        assert result.ret is pytest.ExitCode.TESTS_FAILED

    @parametrize('entries,args,passed', _no_error.values(), ids=_no_error)
    def test_no_error(self, run, entries, args, passed):
        """Check that valid entries and subset runs are no error."""
        result = run(*args, slow=entries)
        result.assert_outcomes(passed=passed)

    @parametrize('tables,line,absent', _stale.values(), ids=_stale)
    @parametrize('args', SERIAL_AND_XDIST.values(), ids=SERIAL_AND_XDIST)
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
