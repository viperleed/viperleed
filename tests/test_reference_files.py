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
    written.write_text(WRITTEN)
    check_reference(written, REFERENCE)
"""
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
        if '-n' in args else
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
    path.write_text(text)
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
        assert reference.read_text() == _BASE  # Not rewritten

    @pytest.mark.parametrize(
        ('written', 'reference'),
        [
            (_TIMESTAMP_A, _TIMESTAMP_B),
            (_PATH_A, _PATH_B),
            (_TIMESTAMP_A + _PATH_A, _TIMESTAMP_B + _PATH_B),
            ],
        ids=('timestamp', 'seSernelius path', 'both'),
        )
    def test_masked_difference_passes(self, pytester, written, reference):
        """Check that differences in the masked parts are ignored."""
        ref = _reference(pytester, reference)
        _run(pytester, written, ref).assert_outcomes(passed=1)
        assert ref.read_text() == reference  # Not rewritten

    def test_mask_is_narrow(self, pytester):
        """Check that other timestamp-like text is not masked."""
        reference = _reference(pytester, 'x 261023-142446 y\n')
        result = _run(pytester, 'x 261024-142446 y\n', reference)
        result.assert_outcomes(failed=1)

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
        assert reference.read_text() == written  # Unmasked

    def test_rewrite_creates_missing(self, pytester):
        """Check that the option creates a missing reference."""
        reference = pytester.path / 'new' / 'folder' / 'reference.txt'
        result = _run(pytester, _BASE, reference, _REWRITE)
        result.assert_outcomes(skipped=1)
        assert reference.read_text() == _BASE

    def test_rewrite_equal_is_a_pass(self, pytester):
        """Check that an up-to-date reference is not touched."""
        reference = _reference(pytester, _TIMESTAMP_B)
        result = _run(pytester, _TIMESTAMP_A, reference, _REWRITE)
        result.assert_outcomes(passed=1)
        assert reference.read_text() == _TIMESTAMP_B

    @pytest.mark.slow  # Subprocess, about 1 s
    @pytest.mark.usefixtures('subprocess_path')
    def test_rewrite_xdist(self, pytester):
        """Check rewriting is reported with pytest-xdist."""
        reference = _reference(pytester, _BASE)
        result = _run(pytester, 'other\n', reference, _REWRITE, '-n', '2')
        result.assert_outcomes(passed=0, skipped=1)
        result.stdout.fnmatch_lines([
            '*rewritten reference files*',
            f'{reference}',
            ])
        assert reference.read_text() == 'other\n'

    @pytest.mark.slow  # Subprocess, about 1 s
    @pytest.mark.usefixtures('subprocess_path')
    def test_mismatch_xdist(self, pytester):
        """Check a mismatch fails with pytest-xdist."""
        reference = _reference(pytester, _BASE)
        result = _run(pytester, 'other\n', reference, '-n', '2')
        result.assert_outcomes(failed=1)
