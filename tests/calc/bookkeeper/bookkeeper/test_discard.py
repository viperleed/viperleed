"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for running bookkeeper in DISCARD mode.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

from copy import deepcopy

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.constants import EDITED_SUFFIX
from viperleed.calc.bookkeeper.constants import ORI_SUFFIX
from viperleed.calc.bookkeeper.constants import STATE_FILES
from viperleed.calc.bookkeeper.history.entry.notes_field import _DISCARDED
from viperleed.calc.bookkeeper.history.info import HistoryInfoFile
from viperleed.calc.bookkeeper.mode import BookkeeperMode as Mode
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP

from .run_bookkeeper_base import _TestBookkeeperRunBase
from .test_clear import TestBookkeeperClear as _TestClear


@fixture(name='after_discard')
def fixture_after_discard(after_archive, after_bookkeper_run):
    """Prepare a directory like the one after DISCARD was executed."""
    return after_bookkeper_run(after_archive, Mode.DISCARD)


@fixture
def fixture_dummy_edited_file(tmp_path):
    """Add a dummy _edited-suffixed file to tmp_path."""
    edited_file = tmp_path / f'some_dummy_file_that_was{EDITED_SUFFIX}'
    edited_file.touch()


with_dummy_edited_file = pytest.mark.usefixtures('fixture_dummy_edited_file')


class TestBookkeeperDiscard(_TestBookkeeperRunBase):
    """Tests for correct behavior of DISCARD bookkeeper runs."""

    mode = Mode.DISCARD

    def check_domains_discarded(self, domains_run):
        """Ensure the root and all domains have been discarded."""
        *main_run, domains = domains_run
        *_, mocker = main_run
        self.check_root_discarded(*main_run)
        # Repeat the check for all subdomains
        for domain_path in domains:
            mock_bookie = mocker.MagicMock(cwd=domain_path)
            mock_bookie.history.info = info = HistoryInfoFile(domain_path)
            info.read()
            self.check_root_discarded(mock_bookie)

    def check_last_entry_discarded(self, bookkeeper, *_):
        """Ensure the last history.info entry has a DISCARDED tag."""
        info = bookkeeper.history.info
        assert info.last_entry_was_discarded
        assert _DISCARDED in info.path.read_text()

    def check_root_discarded(self, bookkeeper, *_):
        """Ensure the most recent run in bookkeeper.cwd is discarded."""
        self.check_root_reverted_to_previous_calc_run(bookkeeper)
        # A 'DISCARDED' note should be in history.info
        self.check_last_entry_discarded(bookkeeper)

    @with_dummy_edited_file
    def test_run_before_calc_exec(self, before_calc_execution, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        self.run_before_calc_exec_and_check(before_calc_execution, caplog)
        self.check_no_warnings(
            caplog,
            exclude_msgs=(
                'Failed to mark as discarded the last entry',
                'user-edited',
                ),
            )
        self.check_complained_about_edited(caplog)

    @with_dummy_edited_file
    def test_discard_after_archive(self, after_archive, caplog):
        """Check reverting of state when DISCARDing an ARCHIVEd calc run."""
        self.run_after_archive_and_check(after_archive, caplog)
        self.check_root_discarded(*after_archive)

        # Some history.info fields are knowingly faulty,
        # but we can still DISCARD the entry.
        faulty_entry_logs = (
            'Found entry with',
            'Could not understand',
            'Faulty entry is',
            'metadata',
            'user-edited',
            )
        self.check_no_warnings(caplog, exclude_msgs=faulty_entry_logs)
        self.check_complained_about_edited(caplog)

    def test_discard_after_archive_domains(self, archived_domains, caplog):
        """Check behavior of DISCARD after ARCHIVE in a multi-domain tree."""
        *main_run, _ = archived_domains
        *_, mocker = main_run
        self.patch_for_domains(mocker)
        self.run_after_archive_and_check(main_run, caplog)
        self.check_domains_discarded(archived_domains)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))

    @with_dummy_edited_file
    def test_discard_after_calc_exec(self, after_calc_execution,
                                     caplog, mocker):
        """Check behavior of DISCARD after a non-ARCHIVEd calc run.

        This may happen, for example, if the previous (calc or
        bookkeeper) execution crashed.

        Parameters
        ----------
        after_calc_execution : fixture
            A bookkeeper and information on a root directory right
            after viperleed.calc has run, and before any bookkeeper
            execution (even the default --archive).
        caplog : fixture
            The pytest.caplog fixture.
        mocker : fixture
            The pytest-mock mocker fixture.

        Returns
        -------
        None.
        """
        bookkeeper, *_ = after_calc_execution
        mock_discard = mocker.patch.object(bookkeeper.history.info,
                                           'discard_last_entry')
        self.run_and_check_prerun_archiving(
            after_calc_execution,
            caplog,
            exclude_warnings=('user-edited',),
            )

        # A 'DISCARDED' note should be in history.info...
        self.check_last_entry_discarded(*after_calc_execution)
        # ...but should have been added already when archiving, not
        # as a result of a call to history.info.discard_last_entry.
        mock_discard.assert_not_called()
        self.check_complained_about_edited(caplog)

    @staticmethod
    def get_tree_from_archived(archived):
        """Return a dict of the root after DISCARD from the ARCHIVE one."""
        discarded = deepcopy(archived)
        _deleted = (
            DEFAULT_OUT,
            DEFAULT_SUPP,
            *(file for file in discarded if file.endswith('.log')),
            )
        for file in _deleted:
            try:
                del discarded[file]
            except KeyError:
                pass
        for file in STATE_FILES:
            try:
                discarded[file] = discarded.pop(f'{file}{ORI_SUFFIX}')
            except KeyError:
                pass
        return discarded

    def test_discard_after_calc_with_edited(self,
                                            after_calc_with_edited_file,
                                            caplog):
        """Check root after DISCARDing a non-ARCHIVEd run with user edits."""
        bookkeeper, _, archived = after_calc_with_edited_file
        self._run_bookkeeper(bookkeeper, {}, caplog)

        # See viperleed/pull/198#issuecomment-2508005204.
        # DISCARD reverts everything, except for adding a history entry
        expect = self.get_tree_from_archived(archived)
        # However, we don't pull back original_inputs of _edited files.
        # This is the only difference from an ARCHIVE-then-DISCARD case
        del expect['POSCAR']
        assert self.collect_root_contents(bookkeeper) == expect

    test_discard_archived_with_edited = (
        _TestClear.test_clear_archived_with_edited
        )

    def test_discard_domains(self, domains_after_calc_execution, caplog):
        """Check behavior of DISCARD after a non-ARCHIVEd domain run."""
        bookkeeper, *_, mocker, _ = domains_after_calc_execution
        mock_discard = mocker.patch.object(bookkeeper.history.info,
                                           'discard_last_entry')
        self.run_and_check_prerun_archiving_domains(
            domains_after_calc_execution,
            caplog,
            )
        self.check_domains_discarded(domains_after_calc_execution)
        mock_discard.assert_not_called()

    @parametrize(mode=(Mode.ARCHIVE, Mode.CLEAR, Mode.DISCARD))
    def test_one_domain_already_processed(self,
                                          mode,
                                          manual_run_one_domain,
                                          caplog):
        """Check DISCARD when one domain subfolder was manually processed."""
        domain_run, already_processed = manual_run_one_domain(mode)
        pre_discarded = mode is Mode.DISCARD
        *_, mocker, _ = domain_run
        caplog.clear()
        self.run_and_check_prerun_archiving_domains(
            domain_run,
            caplog,
            exclude_warnings=('already discarded',) if pre_discarded else (),
            already_processed=(already_processed,),
            )
        self.check_domains_discarded(domain_run)
        # Check also the one that was manually processed.
        mock_bookie = mocker.MagicMock(cwd=already_processed)
        # The root folder should be free of calc results, irrespective
        # of whether the domain was processed or not (since we run in
        # discard mode in the root).
        self.check_root_is_clean(mock_bookie)
        self.check_root_inputs_untouched(mock_bookie)
        self.check_no_edited_files(mock_bookie)

    @with_dummy_edited_file
    def test_twice(self, after_discard, caplog):
        """Check that double discarding does nothing."""
        warnings = (
            'metadata',
            'already discarded',
            'user-edited',
            )
        self.run_again_and_check_nothing_changed(after_discard, caplog,
                                                 acceptable_warnings=warnings)
        self.check_complained_about_edited(caplog)

    test_twice_domains = with_dummy_edited_file(_TestClear.test_twice_domains)
