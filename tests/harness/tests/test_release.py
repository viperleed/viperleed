"""Tests for module release of tests.harness."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

import pytest

from .. import release
from .conftest import MODULE
from .conftest import XDIST_ARGS

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
# Numbers of tests that run, out of test_plain and the cases 'a' and
# 'b-c' of test_tl and test_tl_other
_RELEASE_SELECTION = [
    pytest.param((), 3, id='default'),
    pytest.param(('-m', 'tensorleed'), 2, id='tensorleed'),
    pytest.param(('-m', 'not tensorleed'), 1, id='not-tensorleed'),
    pytest.param(
        ('-m', 'tensorleed or not tensorleed'), 3, id='tensorleed-or-not',
        ),
    pytest.param((_RELEASE_OPTION,), 5, id='release'),
    pytest.param(
        (_RELEASE_OPTION, '-m', 'tensorleed'), 4, id='release-tensorleed',
        ),
    pytest.param((_RELEASE_OPTION, '--quick'), 5, id='release-quick'),
    ]


class TestReleaseOption:
    """Tests for the release marker and the --release option."""

    @pytest.mark.parametrize(
        ('args', 'hint'),
        [
            pytest.param((MODULE + 'test_tl[b-c]',), True, id='node-id'),
            pytest.param(('-m', 'release'), True, id='marker'),
            pytest.param((MODULE + 'test_tl[a]',), False, id='no-release'),
            pytest.param((), False, id='some-ran'),
            pytest.param(
                (MODULE + 'test_tl[b-c]', _RELEASE_OPTION), False,
                id='with-option',
                ),
            ],
        )
    def test_hint_if_nothing_ran(self, run, args, hint):
        """Check the hint if only release cases were deselected."""
        result = run(*args, release=_RELEASE, files=_RELEASE_FILES)
        assert ('release case(s) deselected' in result.stdout.str()) is hint

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

    @pytest.mark.parametrize(('args', 'passed'), _RELEASE_SELECTION)
    @pytest.mark.parametrize(
        'xdist',
        [
            pytest.param((), id='serial'),
            pytest.param(XDIST_ARGS, marks=pytest.mark.slow, id='xdist'),
            ],
        )
    def test_selection(self, run, args, passed, xdist):
        """Check that release cases only run with --release."""
        result = run(*args, *xdist, release=_RELEASE, files=_RELEASE_FILES)
        result.assert_outcomes(passed=passed)
