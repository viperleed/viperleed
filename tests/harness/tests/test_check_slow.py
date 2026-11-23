"""Tests for module check_slow of tests.harness."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-23'
__license__ = 'GPLv3+'

import pytest

from .. import check_slow
from .. import tiers
from .conftest import XDIST_ARGS

_BELOW = 'took less than the lower limit'
_DURATIONS = 'test_durations.py::'
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

_THRESHOLD = tiers.SLOW_THRESHOLD
_LOWER = _THRESHOLD * (1 - tiers.SLOW_BAND)
_MARGIN = 0.01 * _THRESHOLD
_UPPER = _THRESHOLD * (1 + tiers.SLOW_BAND)
_BAND_DURATIONS = {  # Names starting with 'marked' are marked slow
    'unmarked_above': _UPPER + _MARGIN,
    'unmarked_below': _UPPER - _MARGIN,
    'unmarked_at_threshold': _THRESHOLD,
    'marked_above': _LOWER + _MARGIN,
    'marked_below': _LOWER - _MARGIN,
    'marked_at_threshold': _THRESHOLD,
    }


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

    @pytest.mark.slow  # About 1.5 s (pytester)
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

    @pytest.mark.slow  # About 1.5 s (pytester)
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

    @pytest.mark.slow  # About 1.5 s (pytester)
    def test_without_option(self, run):
        """Check that nothing is reported without the option."""
        result = run(files=_SLOW_FILES)
        result.assert_outcomes(passed=4)
        assert result.ret is pytest.ExitCode.OK
        assert check_slow.CHECK_SLOW_OPTION.lstrip('-') not in (
            result.stdout.str()
            )

    @pytest.mark.slow  # About 1.5 s (pytester)
    def test_xdist(self, run):
        """Check that --check-slow also works with pytest-xdist."""
        result = run(
            check_slow.CHECK_SLOW_OPTION,
            *XDIST_ARGS,
            files=_SLOW_FILES,
            )
        self._check_report(result)


class TestCheckSlowBand:
    """Tests for the band of --check-slow around the threshold."""

    @staticmethod
    def _candidates():
        """Return the candidates for the durations of the band."""
        # pylint: disable-next=protected-access
        recorder = check_slow._DurationRecorder()
        recorder.durations.update(_BAND_DURATIONS)
        recorder.slow_ids.update(
            n for n in _BAND_DURATIONS if n.startswith('marked')
            )
        return recorder.candidates(tiers.SLOW_THRESHOLD)

    def test_to_mark(self):
        """Check that unmarked tests are listed from the upper limit."""
        to_mark, _ = self._candidates()
        assert [nodeid for _, nodeid in to_mark] == ['unmarked_above']

    def test_to_unmark(self):
        """Check that marked tests are listed below the lower limit."""
        _, to_unmark = self._candidates()
        assert [nodeid for _, nodeid in to_unmark] == ['marked_below']
