"""Tests for the --quick option of the test-suite configuration."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-22'
__license__ = 'GPLv3+'

import re

import pytest

from .. import quick
from .conftest import INI
from .conftest import XDIST_ARGS
from .conftest import run_pytest

_TEST_FILE = """
import pytest

def test_plain():
    pass

@pytest.mark.slow
def test_slow():
    pass

@pytest.mark.tensorleed
def test_tensorleed():
    pass

@pytest.mark.slow
@pytest.mark.tensorleed
def test_slow_tensorleed():
    pass
"""
_ALL = frozenset(('plain', 'slow', 'tensorleed', 'slow_tensorleed'))
_FAST = {'plain', 'tensorleed'}
_PASSED = re.compile(r'::test_(\w+) PASSED|PASSED \S*::test_(\w+)')
_QUICK = quick.QUICK_OPTION


def _selected(pytester, *args):
    """Return the names of the tests that run with `args`."""
    pytester.makeini(INI)
    pytester.makepyfile(_TEST_FILE)
    result = run_pytest(pytester, *args, '-v')
    passed = (_PASSED.search(line) for line in result.outlines)
    return {match.group(1) or match.group(2) for match in passed if match}


@pytest.mark.parametrize(
    ('args', 'expected'),
    [
        ((), _ALL),
        (('-m', ''), _ALL),
        ((_QUICK,), _FAST),
        ((_QUICK, '-m', 'tensorleed'), {'tensorleed'}),
        ((_QUICK, '-m', ''), _FAST),
        (('-m', 'tensorleed'), {'tensorleed', 'slow_tensorleed'}),
        ((_QUICK, '-m', 'not tensorleed'), {'plain'}),
        ((_QUICK, '-m', 'tensorleed or not tensorleed'), _FAST),
        ],
    ids=(
        'default',
        'empty expression',
        'quick',
        'quick and tensorleed',
        'quick and empty expression',
        'tensorleed',
        'quick and not tensorleed',
        'quick and or expression',
        ),
    )
def test_quick_selection(pytester, args, expected):
    """Check that --quick skips slow tests, on top of any -m."""
    assert _selected(pytester, *args) == expected


@pytest.mark.slow  # About 0.7 s (pytester, subprocess)
def test_quick_selection_xdist(pytester):
    """Check that --quick also skips slow tests with pytest-xdist."""
    selected = _selected(pytester, _QUICK, *XDIST_ARGS)
    assert selected == {'plain', 'tensorleed'}
