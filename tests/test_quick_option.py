"""Tests for the --quick option of the test-suite configuration."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-22'
__license__ = 'GPLv3+'

from pathlib import Path
import re

import pytest

from . import conftest

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
_INI = """
[pytest]
markers =
    slow
    tensorleed
"""
_ALL = frozenset(('plain', 'slow', 'tensorleed', 'slow_tensorleed'))
_PASSED = re.compile(r'::test_(\w+) PASSED|PASSED \S*::test_(\w+)')
_QUICK = conftest.QUICK_OPTION


def _selected(pytester, *args, run=None):
    """Return the names of the tests that run with `args`."""
    pytester.makeini(_INI)
    pytester.makepyfile(_TEST_FILE)
    common = ('-v', '-p', 'no:cacheprovider')
    result = (
        run(*args, *common, '-p', conftest.__name__)
        if run else
        pytester.runpytest_inprocess(
            *args,
            *common,
            plugins=[conftest],
            )
        )
    passed = (_PASSED.search(line) for line in result.outlines)
    return {match.group(1) or match.group(2) for match in passed if match}


@pytest.mark.parametrize(
    ('args', 'expected'),
    [
        ((), _ALL),
        (('-m', ''), _ALL),
        ((_QUICK,), {'plain', 'tensorleed'}),
        ((_QUICK, '-m', 'tensorleed'), {'tensorleed'}),
        ((_QUICK, '-m', ''), {'plain', 'tensorleed'}),
        (('-m', 'tensorleed'), {'tensorleed', 'slow_tensorleed'}),
        ((_QUICK, '-m', 'not tensorleed'), {'plain'}),
        ],
    ids=(
        'default',
        'empty expression',
        'quick',
        'quick and tensorleed',
        'quick and empty expression',
        'tensorleed',
        'quick and not tensorleed',
        ),
    )
def test_quick_selection(pytester, args, expected):
    """Check that --quick skips slow tests, on top of any -m."""
    assert _selected(pytester, *args) == expected


def test_quick_selection_xdist(pytester, monkeypatch):
    """Check that --quick also skips slow tests with pytest-xdist."""
    pytest.importorskip('xdist')
    root = Path(__file__).resolve().parent.parent
    monkeypatch.setenv('PYTHONPATH', str(root))
    selected = _selected(
        pytester,
        _QUICK,
        '-n',
        '2',
        run=pytester.runpytest_subprocess,
        )
    assert selected == {'plain', 'tensorleed'}
