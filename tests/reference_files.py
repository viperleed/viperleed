"""Pytest plugin to compare written files to stored references.

Provides the fixture `check_reference`, which compares a file written
by a test to its stored reference after masking the pieces that are not
reproducible (see MASKS), and the option --rewrite-references.

With the option, a reference that differs from the written file (or
does not exist yet) is overwritten with the written file instead of
failing. The test then does not pass silently: if it would otherwise
pass, its outcome is reported as skipped ("reference rewritten:
<path>"), and the terminal summary lists all the rewritten files. The
paths travel as user_properties of the test report, so that the list
is complete also with pytest-xdist. Review the diff of the rewritten
references (e.g., with git) before committing them.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

import difflib
from pathlib import Path
import re

import pytest

REWRITE_OPTION = '--rewrite-references'
REWRITTEN_KEY = 'rewritten_reference'  # In user_properties of a report
_ENCODING = 'utf-8'

# Masks are (compiled pattern, replacement) pairs. They cover only the
# pieces of the written input files that depend on the time or on the
# location of the run:
# - The timestamp (rp.timestamp, format %y%m%d-%H%M%S) that ends the
#   header line "<systemName> <timestamp>" of delta-input, of
#   AUXLATGEO (refcalc-FIN), and of line 2 of eeasisss-input.
# - The path of the seSernelius file in eeasisss-input, which is
#   relative or absolute depending on the depth of the working folder.
TIMESTAMP_MASK = (
    re.compile(r' \d{6}-\d{6}$', re.MULTILINE),
    ' <timestamp>',
    )
SESERNELIUS_MASK = (
    re.compile(r"^'[^\n]*seSernelius' \|exchange-correlation file$",
               re.MULTILINE),
    "'<seSernelius>' |exchange-correlation file",
    )
MASKS = (TIMESTAMP_MASK, SESERNELIUS_MASK)


def pytest_addoption(parser):
    """Add the --rewrite-references option."""
    parser.addoption(
        REWRITE_OPTION,
        action='store_true',
        help=('Overwrite the stored reference files that differ from '
              'what the tests write, instead of failing. Tests that '
              'rewrite a reference are reported as skipped, and the '
              'rewritten files are listed at the end. Review the '
              'changes before committing them.'),
        )


def apply_masks(text, masks):
    """Return `text` with the (pattern, replacement) `masks` applied."""
    for pattern, replacement in masks:
        text = pattern.sub(replacement, text)
    return text


@pytest.fixture(name='check_reference')
def fixture_check_reference(request):
    """Return a function to compare a written file to its reference."""
    rewrite = request.config.getoption('rewrite_references')

    def _check(written, reference, masks=MASKS):
        """Check that `written` equals `reference` after masking.

        Parameters
        ----------
        written : str or Path
            The file written by the test.
        reference : str or Path
            The stored reference file.
        masks : Sequence, optional
            (compiled pattern, replacement) pairs applied to both
            texts before comparing. Default is MASKS.

        Raises
        ------
        pytest.fail.Exception
            If the masked texts differ, or the reference is missing
            and --rewrite-references was not given.
        """
        written, reference = Path(written), Path(reference)
        written_text = written.read_text(encoding=_ENCODING)
        if reference.is_file():
            reference_text = reference.read_text(encoding=_ENCODING)
        else:
            reference_text = None
        if reference_text is not None:
            masked_ref = apply_masks(reference_text, masks)
            masked_written = apply_masks(written_text, masks)
            if masked_ref == masked_written:
                return
        if rewrite:
            reference.parent.mkdir(parents=True, exist_ok=True)
            reference.write_text(written_text, encoding=_ENCODING)
            request.node.user_properties.append(
                (REWRITTEN_KEY, str(reference))
                )
            return
        if reference_text is None:
            pytest.fail(f'Missing reference file {reference}. Run with '
                        f'{REWRITE_OPTION} to create it from {written}.',
                        pytrace=False)
        diff = difflib.unified_diff(
            masked_ref.splitlines(),
            masked_written.splitlines(),
            fromfile=f'reference: {reference}',
            tofile=f'written: {written}',
            lineterm='',
            )
        pytest.fail('Written file differs from its reference (after '
                    'masking). Run with '
                    f'{REWRITE_OPTION} to update the reference.\n'
                    + '\n'.join(diff),
                    pytrace=False)

    return _check


def _rewritten_paths(report):
    """Return the references that the test of `report` rewrote."""
    return [value for key, value in report.user_properties
            if key == REWRITTEN_KEY]


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Report a test that rewrote a reference as skipped, not passed."""
    outcome = yield
    report = outcome.get_result()
    if call.when != 'call' or not report.passed:
        return
    rewritten = _rewritten_paths(report)
    if not rewritten:
        return
    _, line, _ = item.location
    report.outcome = 'skipped'
    report.longrepr = (
        str(item.path), line + 1,
        'Skipped: reference rewritten: ' + ', '.join(rewritten),
        )


def pytest_terminal_summary(terminalreporter):
    """List the references rewritten in this run."""
    rewritten = sorted({
        path
        for reports in terminalreporter.stats.values()
        for report in reports
        if getattr(report, 'when', None) == 'call'
        for path in _rewritten_paths(report)
        })
    if not rewritten:
        return
    terminalreporter.section('rewritten reference files')
    for path in rewritten:
        terminalreporter.write_line(path)
    terminalreporter.write_line(
        'Review these changes before committing them.'
        )
