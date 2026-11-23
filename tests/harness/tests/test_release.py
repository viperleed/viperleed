"""Tests for module release of tests.harness."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

from pytest_cases import parametrize

from .. import release
from .conftest import MODULE
from .conftest import SERIAL_AND_XDIST

_RELEASE_FILE = """
import pytest

@pytest.fixture(params=['a', 'b-c'])
def tl_case(request):
    return request.param

@pytest.mark.tensorleed
def test_tl(tl_case):
    pass

@pytest.mark.tensorleed
def test_tl_other(tl_case):
    pass

def test_plain():
    pass
"""
_RELEASE_FILES = {'test_sample': _RELEASE_FILE}
_RELEASE = {'tl_case': ('b-c',)}
_RELEASE_OPTION = release.RELEASE_OPTION
_RELEASE_CASES = ('test_tl[b-c]', 'test_tl_other[b-c]')


class TestReleaseOption:
    """Tests for the release marker and the --release option."""

    # Whether the release cases are collected with the arguments
    _collected = {
        'default': ((), False),
        'tensorleed': (('-m', 'tensorleed'), False),
        'release': ((_RELEASE_OPTION,), True),
        'release-tensorleed': (('-m', 'tensorleed', _RELEASE_OPTION), True),
        }
    # Whether a run with the arguments hints at deselected release cases
    _hint = {
        'node-id': ((MODULE + 'test_tl[b-c]',), True),
        'marker': (('-m', 'release'), True),
        'no-release': ((MODULE + 'test_tl[a]',), False),
        'some-ran': ((), False),
        'with-option': ((MODULE + 'test_tl[b-c]', _RELEASE_OPTION), False),
        }
    # Numbers of tests that run, out of test_plain and the cases 'a' and
    # 'b-c' of test_tl and test_tl_other
    _release_selection = {
        'default': ((), 3),
        'tensorleed': (('-m', 'tensorleed'), 2),
        'not-tensorleed': (('-m', 'not tensorleed'), 1),
        'tensorleed-or-not': (('-m', 'tensorleed or not tensorleed'), 3),
        'release': ((_RELEASE_OPTION,), 5),
        'release-tensorleed': ((_RELEASE_OPTION, '-m', 'tensorleed'), 4),
        'release-quick': ((_RELEASE_OPTION, '--quick'), 5),
        }

    @parametrize('args,hint', _hint.values(), ids=_hint)
    def test_hint_if_nothing_ran(self, run, args, hint):
        """Check the hint if only release cases were deselected."""
        result = run(*args, release=_RELEASE, files=_RELEASE_FILES)
        assert ('release case(s) deselected' in result.stdout.str()) is hint

    @parametrize('args,present', _collected.values(), ids=_collected)
    def test_release_case_collected(self, run, args, present):
        """Check that only --release collects the release cases."""
        result = run(
            '--collect-only',
            '-q',
            *args,
            release=_RELEASE,
            files=_RELEASE_FILES,
            )
        out = result.stdout.str()
        assert all((case in out) is present for case in _RELEASE_CASES)

    @parametrize('args,passed', _release_selection.values(),
                 ids=_release_selection)
    @parametrize('xdist', SERIAL_AND_XDIST.values(), ids=SERIAL_AND_XDIST)
    def test_selection(self, run, args, passed, xdist):
        """Check that release cases only run with --release."""
        result = run(*args, *xdist, release=_RELEASE, files=_RELEASE_FILES)
        result.assert_outcomes(passed=passed)
