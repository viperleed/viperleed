"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for running bookkeeper in DISCARD_FULL mode.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

from pytest_cases import parametrize

from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.bookkeeper.bookkeeper import BookkeeperExitCode
from viperleed.calc.bookkeeper.history.errors import MetadataMismatchError
from viperleed.calc.bookkeeper.history.errors import CantRemoveEntryError
from viperleed.calc.bookkeeper.log import BOOKIE_LOGFILE
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.constants import DEFAULT_HISTORY

from ....helpers import filesystem_to_dict
from ..conftest import MOCK_WORKHISTORY
from .run_bookkeeper_base import _TestBookkeeperRunBase
from .test_discard import TestBookkeeperDiscard as _TestDiscard


class TestBookkeeperDiscardFull(_TestBookkeeperRunBase):
    """Tests for correct behavior of DISCARD_FULL bookkeeper runs."""

    mode = BookkeeperMode.DISCARD_FULL

    @staticmethod
    def check_workhistory_folders_deleted(history_path):
        """Ensure "new" workhistory folders are not in `history_path`."""
        for formerly_archived in MOCK_WORKHISTORY.values():
            if not formerly_archived:
                continue
            assert not (history_path/formerly_archived).is_dir()

    _consistency_check_fails = (
        FileNotFoundError,
        MetadataMismatchError,
        CantRemoveEntryError,
        )

    @parametrize(error=_consistency_check_fails)
    def test_consistency_check_fails(self, error, tmp_path, mocker):
        """Check failure if last-entry consistency checks fail."""
        bookkeeper = Bookkeeper(tmp_path)
        # Patch a few methods to be able to test the relevant piece
        # of code. We update_from_cwd here as this creates the handler
        # of history.info, patch its methods, then patch away them
        # .prepare_info_file method to avoid replacing the patched one.
        bookkeeper.update_from_cwd()
        mocker.patch.object(bookkeeper.history.info, 'may_remove_last_entry')
        mocker.patch.object(bookkeeper.history, 'prepare_info_file')
        check = mocker.patch.object(bookkeeper.history,
                                    'check_last_folder_consistent',
                                    side_effect=error)
        exit_code = bookkeeper.run(self.mode)
        check.assert_called_once()
        assert exit_code is not BookkeeperExitCode.SUCCESS

    def test_discard_full_after_archive(self, after_archive, caplog):
        """Check correct behavior of DISCARD_FULL on an ARCHIVEd calc run.

        Should revert state to before the run and purge the history.

        Parameters
        ----------
        after_archive: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        bookkeeper, *_, history_run_path, _ = after_archive
        last_entry = bookkeeper.history.info.last_entry
        bookkeeper.run(mode=self.mode, requires_user_confirmation=False)
        if not last_entry.can_be_removed:
            # Should prevent the removal of the history
            assert history_run_path.is_dir()
            return
        assert not history_run_path.is_dir()
        self.check_root_reverted_to_previous_calc_run(*after_archive)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))

        # "Sibling" folders that were archived from
        # workhistory should also have been removed.
        self.check_workhistory_folders_deleted(bookkeeper.history.path)

    def test_discard_full_after_calc_exec(self, after_calc_execution):
        """Check that a non-ARCHIVEd calc run cannot be DISCARD_FULL-ed."""
        bookkeeper, *_ = after_calc_execution

        # Collect the contents of the root directory before running.
        # Skip the bookkeeper.log, i.e., the only one that changes
        skip = {BOOKIE_LOGFILE}
        before_run = filesystem_to_dict(bookkeeper.cwd, skip=skip)
        exit_code = bookkeeper.run(mode=self.mode,
                                   requires_user_confirmation=False)
        assert exit_code is BookkeeperExitCode.FAIL

        # Now make sure that the contents are identical
        after_run = filesystem_to_dict(bookkeeper.cwd, skip=skip)
        assert after_run == before_run

    @staticmethod
    def get_tree_from_archived(archived):
        """Return a dict of the root after DISCARD from the ARCHIVE one."""
        discarded = _TestDiscard.get_tree_from_archived(archived)
        discarded[DEFAULT_HISTORY] = {}
        return discarded

    def test_discard_full_archived_with_edited(self,
                                               after_calc_with_edited_file,
                                               after_bookkeper_run,
                                               caplog):
        """Check root after DISCARDing an ARCHIVEd run with user edits."""
        # Run in archive mode first to mimic the root structure
        bookkeeper, _, archived = after_bookkeper_run(
            after_calc_with_edited_file,
            'archive',
            )
        caplog.clear()  # The ARCHIVE logs
        bookkeeper.run(mode=self.mode, requires_user_confirmation=False)

        # See viperleed/pull/198#issuecomment-2506549827
        expect = self.get_tree_from_archived(archived)
        assert self.collect_root_contents(bookkeeper) == expect

    def test_discard_full_domains(self, archived_domains, caplog):
        """Check correct behavior of DISCARD_FULL on an ARCHIVEd domain run."""
        *main_run, domains = archived_domains
        *_, mocker = main_run
        self.patch_for_domains(mocker)
        self.test_discard_full_after_archive(main_run, caplog)

        # And repeat checks also for the domain subfolders
        for domain_path, domain_info in domains.items():
            history_run_path = domain_info['history_run']
            assert not history_run_path.is_dir()
            mock_bookie = mocker.MagicMock(cwd=domain_path)
            self.check_root_reverted_to_previous_calc_run(mock_bookie)
            self.check_workhistory_folders_deleted(history_run_path.parent)

    def test_run_before_calc_exec(self, before_calc_execution, caplog):
        """Check correct overwriting of input files in DISCARD_FULL mode.

        This should do the same as a normal DISCARD.

        Parameters
        ----------
        before_calc_execution: fixture
        caplog: fixture

        Returns
        -------
        None.
        """
        self.run_before_calc_exec_and_check(before_calc_execution,
                                            caplog,
                                            requires_user_confirmation=False)
        self.check_no_warnings(
            caplog,
            exclude_msgs=('remove',)   # Catches two expected warnings
            )

    @parametrize(confirmed=(True, False))
    def test_user_confirmation(self, confirmed, after_archive, mocker):
        """Check successful execution when user does/does not confirm."""
        bookkeeper, *_ = after_archive
        mocker.patch.object(bookkeeper,
                            '_user_confirmed',
                            return_value=confirmed)
        # Patch away stuff that may prevent removal,
        # as we only want to check the exit code
        mocker.patch.object(bookkeeper, '_check_may_discard_full')
        mocker.patch.object(type(bookkeeper.history.info), 'remove_last_entry')

        # The default is asking for user confirmation. Do not purposely
        # provide an explicit keyword argument to test the default.
        code = bookkeeper.run(mode=self.mode)
        expect = (BookkeeperExitCode.SUCCESS if confirmed
                  else BookkeeperExitCode.NOTHING_TO_DO)
        assert code is expect
