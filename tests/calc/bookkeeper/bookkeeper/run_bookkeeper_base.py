"""Module run_bookkeeper_base of tests/calc/bokkeeper/bookkeeper.

Defines the _TestBookkeeperRunBase class that collects basic
functionality for testing execution of the bookkeeper in
various modes.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-11-15'
__license__ = 'GPLv3+'

import logging

from viperleed.calc.bookkeeper.bookkeeper import BookkeeperExitCode
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history.meta import _METADATA_NAME
from viperleed.calc.bookkeeper.log import BOOKIE_LOGFILE
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP

from ..conftest import MOCK_INPUT_CONTENT
from ..conftest import MOCK_OUT_CONTENT
from ..conftest import MOCK_ORIG_CONTENT
from ..conftest import MOCK_STATE_FILES
from ..conftest import MOCK_WORKHISTORY


class _TestBookkeeperRunBase:
    """Base class for checking correct execution of bookkeeper."""

    mode = None

    def check_metadata_exists(self, history_folder):
        """Test that the metadata file is present in `history_folder`."""
        assert (history_folder / _METADATA_NAME).is_file()

    def check_input_files_in_history(self, *run, check_input_contents=True):
        """Make sure that input files were stored in history."""
        *_, history_run_path = run
        expected_contents = (MOCK_INPUT_CONTENT,)
        if self.mode.uses_ori_files_as_fallback:
            expected_contents += (MOCK_ORIG_CONTENT,)
        for file in MOCK_STATE_FILES:
            archived_input = history_run_path / file
            assert archived_input.is_file()
            if check_input_contents:
                contents = archived_input.read_text()
                assert any(e in contents for e in expected_contents)

    def check_history_exists(self, bookkeeper, history_run_path):
        """Test that history_path and directory/history.info exist."""
        assert history_run_path.is_dir()
        assert (bookkeeper.cwd / HISTORY_INFO_NAME).exists()
        self.check_metadata_exists(history_run_path)

    def check_history_folder_empty(self, bookkeeper, *_):
        """Test that an empty history folder exists, except bookkeeper.log."""
        cwd = bookkeeper.cwd
        assert (cwd / DEFAULT_HISTORY).exists()
        history_contents = (f for f in (cwd / DEFAULT_HISTORY).iterdir()
                            if f.name != BOOKIE_LOGFILE)
        assert not any(history_contents)

    def check_no_warnings(self, caplog, contain=None, exclude_msgs=()):
        """Check that there are no warnings or errors."""
        if contain is None:
            def _record_faulty(record):
                return record.levelno >= logging.WARNING
        else:
            def _record_faulty(record):
                return (record.levelno >= logging.WARNING
                        and contain in str(record))
        faulty = [
            rec
            for rec in caplog.records
            if not any(exclude in str(rec) for exclude in exclude_msgs)
            and _record_faulty(rec)
            ]
        assert not any(faulty), f'Found: {faulty[0].getMessage()!r}'

    def check_out_files_in_history(self, *run):
        """Check that the expected state files are stored in 'OUT'."""
        *_, history_run_path = run
        for file in MOCK_STATE_FILES:
            hist_file = history_run_path / DEFAULT_OUT / f'{file}_OUT'
            assert hist_file.is_file()
            assert MOCK_OUT_CONTENT in hist_file.read_text()

    def check_out_files_untouched(self, bookkeeper, *_):
        """Ensure all expected files are found in OUT."""
        out = bookkeeper.cwd / DEFAULT_OUT
        for file in MOCK_STATE_FILES:
            out_file = out / f'{file}_OUT'
            assert out_file.is_file
            assert MOCK_OUT_CONTENT in out_file.read_text()

    def check_root_after_archive(self, *after_archive,
                                 check_input_contents=True):
        """Make sure the root is structured as expected after archiving."""
        self.check_root_inputs_renamed_to_ori(
            *after_archive,
            check_input_contents=check_input_contents
            )
        self.check_out_files_untouched(*after_archive)
        if check_input_contents:
            self.check_root_inputs_replaced_by_out(*after_archive)

    def check_root_inputs_renamed_to_ori(self, bookkeeper, *_,
                                         check_input_contents=True):
        """Check that the input files have now a _ori suffix."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            ori_contents = (cwd / f'{file}_ori').read_text()
            assert MOCK_INPUT_CONTENT in ori_contents
            if not check_input_contents:
                continue
            input_content = (cwd / file).read_text()
            assert MOCK_OUT_CONTENT in input_content

    def check_root_inputs_replaced_by_out(self, bookkeeper, *_):
        """Check that the input files in root come from OUT."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            out_content = (cwd / file).read_text()
            assert MOCK_OUT_CONTENT in out_content

    def check_root_inputs_untouched(self, bookkeeper, *_):
        """Check the the original state files have not been moved."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            assert (cwd / file).is_file()
            input_content = (cwd / file).read_text()
            assert MOCK_INPUT_CONTENT in input_content

    def check_root_is_clean(self, bookkeeper, *_):
        """Check that no calc output is present in the main directory."""
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            assert not (cwd / f'{file}_ori').is_file()
        assert not (cwd / DEFAULT_SUPP).exists()
        assert not (cwd / DEFAULT_OUT).exists()
        assert not any(cwd.glob('*.log'))

    def check_workhistory_archived(self, bookkeeper, *_):
        """Ensure the workhistory folders have been stored correctly."""
        history = bookkeeper.history.path
        # pylint: disable-next=protected-access           # OK in tests
        workhistory = bookkeeper._workhistory.path
        assert not workhistory.is_dir()
        for ori_name, hist_name in MOCK_WORKHISTORY.items():
            if hist_name is None:  # File should be deleted
                continue
            moved_dir = history/hist_name
            moved_file = moved_dir/'file'
            assert moved_dir.is_dir()
            assert moved_file.is_file()
            assert ori_name in moved_file.read_text()

    def run_after_archive_and_check(self, after_archive,
                                    check_input_contents=True,
                                    **kwargs):
        """Check that running bookkeeper after ARCHIVE does basic stuff."""
        bookkeeper, *_ = after_archive
        # bookkeeper should not think that it needs archiving
        assert not bookkeeper.archiving_required
        kwargs.setdefault('mode', self.mode)
        bookkeeper.run(**kwargs)
        bookkeeper.update_from_cwd(silent=True)
        self.check_history_exists(*after_archive)
        self.check_out_files_in_history(*after_archive)
        self.check_input_files_in_history(
            *after_archive,
            check_input_contents=check_input_contents,
            )
        if self.mode is BookkeeperMode.ARCHIVE:
            self.check_root_after_archive(
                *after_archive,
                check_input_contents=check_input_contents,
                )
        else:
            self.check_root_is_clean(*after_archive)
        # Check that the workhistory directories are
        # also where they should be
        self.check_workhistory_archived(bookkeeper)

    def run_after_calc_exec_and_check(self, after_calc_execution,
                                      check_archiving_required=True,
                                      **kwargs):
        """Check that running bookkeeper after calc does some basic stuff."""
        bookkeeper, *_ = after_calc_execution
        if check_archiving_required:
            # bookkeeper should think that it needs archiving
            assert bookkeeper.archiving_required
        kwargs.setdefault('mode', self.mode)
        bookkeeper.run(**kwargs)
        bookkeeper.update_from_cwd(silent=True)
        self.check_history_exists(*after_calc_execution)
        self.check_out_files_in_history(*after_calc_execution)
        self.check_workhistory_archived(*after_calc_execution)

    def run_and_check_prerun_archiving(self, after_calc_execution, caplog):
        """Execute bookkeeper, and verify it has archived a calc run."""
        self.run_after_calc_exec_and_check(after_calc_execution)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))
        self.check_root_is_clean(*after_calc_execution)

        # Original SHOULD NOT be replaced by output:
        # ARCHIVE does not run only if the run crashed,
        # in which case we don't want to overwrite
        self.check_root_inputs_untouched(*after_calc_execution)

    def run_archive_after_calc_and_check(self, after_calc_execution, caplog,
                                         check_archiving_required=True):
        """Check correct storage of history files in ARCHIVE mode."""
        kwargs = {
            'check_archiving_required': check_archiving_required,
            'mode': 'archive',
            }
        self.run_after_calc_exec_and_check(after_calc_execution, **kwargs)
        self.check_root_after_archive(*after_calc_execution)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))

    def run_before_calc_exec_and_check(self, before_calc_execution,
                                       check_archiving_required=True,
                                       **kwargs):
        """Check that running bookkeeper before calc does almost nothing."""
        bookkeeper = before_calc_execution
        if check_archiving_required:
            # bookkeeper should not think that it needs archiving
            assert not bookkeeper.archiving_required
        kwargs.setdefault('mode', self.mode)
        exit_code = bookkeeper.run(**kwargs)
        # Bookkeeper should not do anything (except for logging)
        self.check_history_folder_empty(before_calc_execution)
        self.check_root_inputs_untouched(before_calc_execution)
        assert exit_code is not BookkeeperExitCode.FAIL
