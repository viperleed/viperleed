"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for running bookkeeper in ARCHIVE mode.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.meta import _METADATA_NAME
from viperleed.calc.bookkeeper.errors import BookkeeperUnexpectedError
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.constants import DEFAULT_HISTORY
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP
from viperleed.calc.constants import DEFAULT_WORK
from viperleed.calc.constants import LOG_PREFIX

from ....helpers import filesystem_from_dict
from ..conftest import MOCK_ORIG_CONTENT
from .run_bookkeeper_base import _TestBookkeeperRunBase


class TestBookkeeperArchive(_TestBookkeeperRunBase):
    """Tests for correct behavior of ARCHIVE bookkeeper runs."""

    mode = BookkeeperMode.ARCHIVE

    def test_archive_after_calc_exec(self,
                                     after_calc_execution,
                                     caplog,
                                     **kwargs):
        """Check correct storage of history files in ARCHIVE mode."""
        self.run_archive_after_calc_and_check(after_calc_execution,
                                              caplog,
                                              **kwargs)

    def test_archive_domains(self,
                             domains_after_calc_execution,
                             caplog,
                             mocker):
        """Check the ARCHIVE behavior without passing explicit domain paths."""
        # This is, e.g., a manual call if the ARCHIVE after calc fails.
        *main_run, domains = domains_after_calc_execution
        self.patch_for_domains(mocker)
        self.run_archive_after_calc_and_check(main_run, caplog)
        self.check_domains_archived(main_run, domains)

    def test_archive_twice(self, after_archive, caplog):
        """Bookkeeper ARCHIVE after ARCHIVE should not do anything."""
        warnings = ('metadata',)
        self.run_again_and_check_nothing_changed(after_archive, caplog,
                                                 acceptable_warnings=warnings)

    def test_archive_twice_domains(self, archived_domains, caplog, mocker):
        """Bookkeeper ARCHIVE after ARCHIVE should not do anything."""
        self.patch_for_domains(mocker)
        self.test_archive_twice(archived_domains, caplog)

    def test_archive_with_edited_file(self,
                                      after_calc_with_edited_file,
                                      caplog):
        """Check expected directory tree after running ARCHIVE."""
        bookkeeper, _, expect = after_calc_with_edited_file
        self._run_bookkeeper(bookkeeper, {}, caplog)
        assert self.collect_root_contents(bookkeeper) == expect

    def test_domains_explicit(self,
                              domains_after_calc_execution,
                              caplog,
                              mocker):
        """Check the result of archiving with explicit domain paths given."""
        *main_run, domains = domains_after_calc_execution
        self.patch_for_domains(mocker)
        self.run_archive_after_calc_and_check(main_run,
                                              caplog,
                                              domains=domains)
        self.check_domains_archived(main_run, domains)

    def test_domains_explicit_no_previous_metadata(
        self,
        domains_after_calc_execution,
        caplog,
        mocker,
        ):
        """Check the result of archiving with explicit domain paths given."""
        main_bookie, *_, domains = domains_after_calc_execution
        # Remove all metadata files from the various history folders
        for root_path in (main_bookie.cwd, *domains):
            history = root_path/DEFAULT_HISTORY
            for folder in history.glob('**'):
                try:
                    (folder/_METADATA_NAME).unlink()
                except FileNotFoundError:
                    pass
        self.test_domains_explicit(domains_after_calc_execution,
                                   caplog,
                                   mocker)

    def test_run_before_calc_exec(self, before_calc_execution, caplog):
        """Check no archiving happens before calc runs."""
        self.run_before_calc_exec_and_check(before_calc_execution, caplog)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))

    @parametrize(workname=(DEFAULT_WORK, 'some_other_folder'))
    def test_issue_344(self, workname, after_calc_execution, caplog):
        """Ensure Issue #344 does not present itself anymore."""
        bookkeeper, *_ = after_calc_execution
        cwd = bookkeeper.cwd

        # The issue came up if (i) a STATE_FILE was missing from OUT
        vibrocc_out = 'VIBROCC'
        try:
            (cwd/DEFAULT_OUT/vibrocc_out).unlink()
        except FileNotFoundError:  # Probably old-style with _OUT
            vibrocc_out_p = next((cwd/DEFAULT_OUT).glob('VIBROCC_OUT*'))
            vibrocc_out_p.unlink()
            vibrocc_out = vibrocc_out_p.name
        # (ii) there was a work directory
        work_tree = {
            f'{LOG_PREFIX}_timestamp.log': '',
            DEFAULT_OUT: {},
            DEFAULT_SUPP: {},
            }
        work_dir = cwd/workname
        work_dir.mkdir()
        filesystem_from_dict(work_tree, work_dir)
        self.test_archive_after_calc_exec(after_calc_execution,
                                          caplog,
                                          missing_out_files={'VIBROCC'})
        self._check_file_contents(cwd/'VIBROCC', MOCK_ORIG_CONTENT)

    def test_raises_on_missing_state_files(self, after_calc_execution, mocker):
        """Check that missing state files after a run cause exceptions."""
        has_out_suffixed = self.has_out_suffixed(*after_calc_execution)
        bookkeeper, *_ = after_calc_execution
        # pylint: disable-next=protected-access           # OK in tests
        mocker.patch.object(bookkeeper._root,
                            'ensure_has_unlabled_inputs',
                            return_value=True)
        with pytest.raises(BookkeeperUnexpectedError):
            bookkeeper.run(self.mode)
        # While the above raises, it should happen only at the end:
        # the rest of the archiving and handling of the root folder
        # should proceed as for a "good" run.
        # Perform the same checks as run_archive_after_calc_and_check
        self.check_has_archived(after_calc_execution,
                                has_out_suffixed=has_out_suffixed,
                                mode=self.mode)
        self.check_root_after_archive(*after_calc_execution,
                                      out_suffixed=has_out_suffixed)

    def test_raises_on_missing_state_files_domains(self,
                                                   domains_after_calc_execution,
                                                   mocker):
        """Check that missing state files after a run cause exceptions."""
        *main_run, domains = domains_after_calc_execution
        self.patch_for_domains(mocker)
        self.test_raises_on_missing_state_files(main_run, mocker)
        self.check_domains_archived(main_run, domains)
