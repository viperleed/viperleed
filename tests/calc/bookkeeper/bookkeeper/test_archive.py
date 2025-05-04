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

from viperleed.calc.bookkeeper.history.meta import _METADATA_NAME
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.constants import DEFAULT_HISTORY

from .run_bookkeeper_base import _TestBookkeeperRunBase


class TestBookkeeperArchive(_TestBookkeeperRunBase):
    """Tests for correct behavior of ARCHIVE bookkeeper runs."""

    mode = BookkeeperMode.ARCHIVE

    def test_archive_after_calc_exec(self, after_calc_execution, caplog):
        """Check correct storage of history files in ARCHIVE mode."""
        self.run_archive_after_calc_and_check(after_calc_execution, caplog)

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
                except OSError:
                    pass
        self.test_domains_explicit(domains_after_calc_execution,
                                   caplog,
                                   mocker)

    def test_run_before_calc_exec(self, before_calc_execution, caplog):
        """Check no archiving happens before calc runs."""
        self.run_before_calc_exec_and_check(before_calc_execution, caplog)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))
