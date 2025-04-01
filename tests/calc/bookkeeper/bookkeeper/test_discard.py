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

from viperleed.calc.bookkeeper.constants import EDITED_SUFFIX
from viperleed.calc.bookkeeper.constants import ORI_SUFFIX
from viperleed.calc.bookkeeper.constants import STATE_FILES
from viperleed.calc.bookkeeper.history.entry.notes_field import _DISCARDED
from viperleed.calc.bookkeeper.mode import BookkeeperMode
from viperleed.calc.constants import DEFAULT_OUT
from viperleed.calc.constants import DEFAULT_SUPP

from .run_bookkeeper_base import _TestBookkeeperRunBase


@fixture(name='after_discard')
def fixture_after_discard(after_archive, after_bookkeper_run):
    """Prepare a directory like the one after DISCARD was executed."""
    return after_bookkeper_run(after_archive, BookkeeperMode.DISCARD)


@fixture
def with_dummy_edited_file(tmp_path):
    """Add a dummy _edited-suffixed file to tmp_path."""
    edited_file = tmp_path / f'some_dummy_file_that_was{EDITED_SUFFIX}'
    edited_file.touch()


class TestBookkeeperDiscard(_TestBookkeeperRunBase):
    """Tests for correct behavior of DISCARD bookkeeper runs."""

    mode = BookkeeperMode.DISCARD

    def check_last_entry_discarded(self, bookkeeper, *_):
        """Ensure the last history.info entry has a DISCARDED tag."""
        info = bookkeeper.history.info
        assert info.last_entry_was_discarded
        assert _DISCARDED in info.path.read_text()

    @pytest.mark.usefixtures('with_dummy_edited_file')
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

    @pytest.mark.usefixtures('with_dummy_edited_file')
    def test_discard_after_archive(self, after_archive, caplog):
        """Check reverting of state when DISCARDing an ARCHIVEd calc run."""
        self.run_after_archive_and_check(after_archive, caplog)
        self.check_root_reverted_to_previous_calc_run(*after_archive)

        # A 'DISCARDED' note should be in history.info
        self.check_last_entry_discarded(*after_archive)

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

    @pytest.mark.usefixtures('with_dummy_edited_file')
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

    def _get_discarded_tree_from_archived(self, archived):
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
        expect = self._get_discarded_tree_from_archived(archived)
        # However, we don't pull back original_inputs of _edited files.
        # This is the only difference from an ARCHIVE-then-DISCARD case
        del expect['POSCAR']
        assert self.collect_root_contents(bookkeeper) == expect

    def test_discard_archived_with_edited(self,
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
        self._run_bookkeeper(bookkeeper, {}, caplog)

        # See viperleed/pull/198#issuecomment-2506549827
        expect = self._get_discarded_tree_from_archived(archived)
        assert self.collect_root_contents(bookkeeper) == expect

    @pytest.mark.usefixtures('with_dummy_edited_file')
    def test_discard_twice(self, after_discard, caplog):
        """Check that double discarding does nothing."""
        warnings = (
            'metadata',
            'already discarded',
            'user-edited',
            )
        self.run_again_and_check_nothing_changed(after_discard, caplog,
                                                 acceptable_warnings=warnings)
        self.check_complained_about_edited(caplog)
