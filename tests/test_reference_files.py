"""Tests for the reference-file comparison of the test suite."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest

from . import conftest
from . import reference_files

_REWRITE = reference_files.REWRITE_OPTION
_TEST_FILE = """
def test_compare(check_reference, tmp_path):
    written = tmp_path / 'written'
    written.write_text(WRITTEN, encoding='utf-8')
    check_reference(written, REFERENCE)
"""
_XDIST_OPTION = '-n'
_CHANGED = 'other\n'
_BASE = 'line 1\nline 2\nline 3\n'
_TIMESTAMP_A = 'Ag(100) 261023-142446\n'
_TIMESTAMP_B = 'Ag(100) 261024-090000\n'
_PATH_A = "'../../seSernelius' |exchange-correlation file\n"
_PATH_B = "'/long/absolute/path/seSernelius' |exchange-correlation file\n"


def _run(pytester, written, reference, *args):
    """Run a test comparing `written` text to the file `reference`."""
    pytester.makepyfile(
        f'WRITTEN = {written!r}\nREFERENCE = {str(reference)!r}\n'
        + _TEST_FILE
        )
    run = (
        pytester.runpytest_subprocess
        if _XDIST_OPTION in args else
        pytester.runpytest_inprocess
        )
    return run(*args, '-p', 'no:cacheprovider', '-p', conftest.__name__)


@pytest.fixture(name='subprocess_path')
def fixture_subprocess_path(monkeypatch):
    """Make the test-suite plugins importable in a subprocess."""
    root = Path(__file__).resolve().parent.parent
    monkeypatch.setenv('PYTHONPATH', str(root))


def _reference(pytester, text):
    """Store `text` as a reference file and return its path."""
    path = pytester.path / 'refs' / 'reference.txt'
    path.parent.mkdir(exist_ok=True)
    path.write_text(text, encoding='utf-8')
    return path


class TestCheckReference:
    """Tests for the check_reference fixture."""

    def test_equal(self, pytester):
        """Check that an identical file passes."""
        reference = _reference(pytester, _BASE)
        _run(pytester, _BASE, reference).assert_outcomes(passed=1)

    def test_mismatch_diff(self, pytester):
        """Check that a mismatch fails with a readable diff."""
        reference = _reference(pytester, _BASE)
        result = _run(pytester, 'line 1\nline X\nline 3\n', reference)
        result.assert_outcomes(failed=1)
        result.stdout.fnmatch_lines([
            f'*--- reference: {reference}',
            '*+++ written: *',
            '*-line 2',
            '*+line X',
            ])
        assert reference.read_text(encoding='utf-8') == _BASE  # Not rewritten

    @pytest.mark.parametrize(
        ('written', 'reference'),
        [
            (_TIMESTAMP_A, _TIMESTAMP_B),
            (_PATH_A, _PATH_B),
            ],
        ids=('timestamp', 'path'),
        )
    def test_timestamp_or_path_difference_fails(self, pytester,
                                                written, reference):
        """Check that timestamps and paths are compared verbatim."""
        ref = _reference(pytester, reference)
        _run(pytester, written, ref).assert_outcomes(failed=1)
        assert ref.read_text(encoding='utf-8') == reference  # Not rewritten

    def test_missing_reference(self, pytester):
        """Check that a missing reference fails, naming the option."""
        reference = pytester.path / 'refs' / 'missing.txt'
        result = _run(pytester, _BASE, reference)
        result.assert_outcomes(failed=1)
        result.stdout.fnmatch_lines([f'*Missing reference file*{_REWRITE}*'])
        assert not reference.exists()

    def test_rewrite(self, pytester):
        """Check that the option rewrites and reports, not passes."""
        reference = _reference(pytester, _BASE)
        written = 'line 1\nline X 261023-142446\n'
        result = _run(pytester, written, reference, _REWRITE, '-rs')
        result.assert_outcomes(passed=0, skipped=1)
        result.stdout.fnmatch_lines([
            '*rewritten reference files*',
            f'*reference rewritten: {reference}*',
            ])
        assert reference.read_text(encoding='utf-8') == written

    def test_rewrite_creates_missing(self, pytester):
        """Check that the option creates a missing reference."""
        reference = pytester.path / 'new' / 'folder' / 'reference.txt'
        result = _run(pytester, _BASE, reference, _REWRITE)
        result.assert_outcomes(skipped=1)
        assert reference.read_text(encoding='utf-8') == _BASE

    def test_rewrite_equal_is_a_pass(self, pytester):
        """Check that an up-to-date reference is not touched."""
        reference = _reference(pytester, _BASE)
        result = _run(pytester, _BASE, reference, _REWRITE)
        result.assert_outcomes(passed=1)
        assert result.ret is pytest.ExitCode.OK
        assert reference.read_text(encoding='utf-8') == _BASE

    @pytest.mark.slow  # Subprocess, about 1 s
    @pytest.mark.usefixtures('subprocess_path')
    def test_rewrite_xdist(self, pytester):
        """Check rewriting is reported with pytest-xdist."""
        reference = _reference(pytester, _BASE)
        result = _run(pytester, _CHANGED, reference, _REWRITE,
                      _XDIST_OPTION, '2')
        result.assert_outcomes(passed=0, skipped=1)
        result.stdout.fnmatch_lines([
            '*rewritten reference files*',
            f'{reference}',
            ])
        assert reference.read_text(encoding='utf-8') == _CHANGED

    @pytest.mark.slow  # Subprocess, about 1 s
    @pytest.mark.usefixtures('subprocess_path')
    def test_mismatch_xdist(self, pytester):
        """Check a mismatch fails with pytest-xdist."""
        reference = _reference(pytester, _BASE)
        result = _run(pytester, _CHANGED, reference, _XDIST_OPTION, '2')
        result.assert_outcomes(failed=1)

    @pytest.mark.usefixtures('subprocess_path')
    @pytest.mark.parametrize(
        'args',
        [
            pytest.param((), id='plain'),
            pytest.param(('--no-summary',), id='no-summary'),
            pytest.param(('-p', 'no:terminal'), id='no-terminal'),
            pytest.param(
                (_XDIST_OPTION, '2'),
                id='xdist',
                marks=pytest.mark.slow,  # Subprocess, about 1 s
                ),
            ],
        )
    def test_rewrite_exit_status(self, pytester, args):
        """Check that a run that rewrote a reference exits non-zero."""
        reference = _reference(pytester, _BASE)
        result = _run(pytester, _CHANGED, reference, _REWRITE, *args)
        assert result.ret is pytest.ExitCode.TESTS_FAILED
        assert reference.read_text(encoding='utf-8') == _CHANGED
