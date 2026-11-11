"""Pytest plugin to compare written files to stored references.

Provides the fixture `check_reference`, which compares a file written
by a test to its stored reference, text for text, and the option
--rewrite-references. Tests must therefore pin whatever depends on the
time or on the location of the run (e.g., timestamps and paths).

With the option, a reference that differs from the written file (or
does not exist yet) is overwritten with the written file instead of
failing. The test then does not pass silently: if it would otherwise
pass, its outcome is reported as skipped ("reference rewritten:
<path>"), and the terminal summary lists all the rewritten files. The
paths travel as user_properties of the test report, so that the list
is complete also with pytest-xdist. A run that rewrote any reference
(and would otherwise succeed) ends with exit status 1, so that a
rewrite is never mistaken for a pass, also with --no-summary and
-p no:terminal. Review the diff of the rewritten
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

import pytest

REWRITE_OPTION = '--rewrite-references'
REWRITTEN_KEY = 'rewritten_reference'  # In user_properties of a report
REWRITE_EXIT_STATUS = pytest.ExitCode.TESTS_FAILED
_ENCODING = 'utf-8'


def pytest_addoption(parser):
    """Add the --rewrite-references option."""
    parser.addoption(
        REWRITE_OPTION,
        action='store_true',
        help=('Overwrite the stored reference files that differ from '
              'what the tests write, instead of failing. Tests that '
              'rewrite a reference are reported as skipped, the '
              'rewritten files are listed at the end, and the run '
              'exits with a non-zero status. Review the changes '
              'before committing them.'),
        )


@pytest.fixture(name='check_reference')
def fixture_check_reference(request):
    """Return a function to compare a written file to its reference."""
    rewrite = request.config.getoption('rewrite_references')

    def _check(written, reference):
        """Check that `written` has the same text as `reference`.

        Parameters
        ----------
        written : str or Path
            The file written by the test.
        reference : str or Path
            The stored reference file.

        Raises
        ------
        pytest.fail.Exception
            If the texts differ, or the reference is missing
            and --rewrite-references was not given.
        """
        written, reference = Path(written), Path(reference)
        written_text = written.read_text(encoding=_ENCODING)
        reference_text = (
            reference.read_text(encoding=_ENCODING)
            if reference.is_file() else
            None
            )
        if reference_text == written_text:
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
            reference_text.splitlines(),
            written_text.splitlines(),
            fromfile=f'reference: {reference}',
            tofile=f'written: {written}',
            lineterm='',
            )
        pytest.fail('Written file differs from its reference. '
                    f'Run with {REWRITE_OPTION} to update the reference.\n'
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


class _RewriteCollector:
    """Collect the rewritten references and set the exit status."""

    def __init__(self):
        """Initialize an empty collection."""
        self.rewritten = set()

    def pytest_runtest_logreport(self, report):
        """Collect the references rewritten by a test."""
        # Also called for the reports that xdist workers send to the
        # controller, and independent of the terminal plugin.
        if report.when == 'call':
            self.rewritten.update(_rewritten_paths(report))

    @pytest.hookimpl(trylast=True)
    def pytest_sessionfinish(self, session, exitstatus):
        """Make a run that rewrote references exit non-zero."""
        if hasattr(session.config, 'workerinput'):  # xdist worker
            return
        if self.rewritten and exitstatus == pytest.ExitCode.OK:
            session.exitstatus = REWRITE_EXIT_STATUS

    def pytest_terminal_summary(self, terminalreporter):
        """List the references rewritten in this run."""
        if not self.rewritten:
            return
        terminalreporter.section('rewritten reference files')
        for path in sorted(self.rewritten):
            terminalreporter.write_line(path)
        terminalreporter.write_line(
            'Review these changes before committing them.'
            )


def pytest_configure(config):
    """Register the collector of rewritten references."""
    config.pluginmanager.register(_RewriteCollector(), 'rewrite-collector')
