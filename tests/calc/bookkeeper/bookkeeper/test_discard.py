"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for running bookkeeper in DISCARD mode.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

from viperleed.calc.bookkeeper.history.entry.notes_field import _DISCARDED
from viperleed.calc.bookkeeper.mode import BookkeeperMode

from ..conftest import MOCK_INPUT_CONTENT
from ..conftest import MOCK_STATE_FILES
from .run_bookkeeper_base import _TestBookkeeperRunBase


class TestBookkeeperDiscard(_TestBookkeeperRunBase):
    """Tests for correct behavior of DISCARD bookkeeper runs."""

    mode = BookkeeperMode.DISCARD

    def test_run_before_calc_exec(self, before_calc_execution, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        self.run_before_calc_exec_and_check(before_calc_execution)
        self.check_no_warnings(
            caplog,
            exclude_msgs=('Failed to mark last entry as discarded',),
            )

    def test_discard_after_archive(self, after_archive, caplog):
        """Check reverting of state when DISCARDing an ARCHIVEd calc run."""
        self.run_after_archive_and_check(after_archive)

        # Original be replaced by output                                        # TODO: this does something else!
        bookkeeper, *_ = after_archive
        cwd = bookkeeper.cwd
        for file in MOCK_STATE_FILES:
            out_content = (cwd / file).read_text()
            assert MOCK_INPUT_CONTENT in out_content
        # A 'DISCARDED' note should be in history.info
        assert _DISCARDED in bookkeeper.history.info.path.read_text()
        # Some fields are knowingly faulty,
        # but we can still DISCARD them.
        faulty_entry_logs = (
            'Found entry with',
            'Could not understand',
            'Faulty entry is',
            'metadata',
            )
        self.check_no_warnings(caplog, exclude_msgs=faulty_entry_logs)

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
        self.run_and_check_prerun_archiving(after_calc_execution, caplog)

        # A 'DISCARDED' note should be in history.info...
        assert bookkeeper.history.info.last_entry_was_discarded
        assert _DISCARDED in bookkeeper.history.info.path.read_text()
        # ...but should have been added already when archiving, not
        # as a result of a call to history.info.discard_last_entry.
        mock_discard.assert_not_called()
