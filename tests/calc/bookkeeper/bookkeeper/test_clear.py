"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for running bookkeeper in CLEAR mode.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

from viperleed.calc.bookkeeper.mode import BookkeeperMode

from .run_bookkeeper_base import _TestBookkeeperRunBase


class TestBookkeeperClear(_TestBookkeeperRunBase):
    """Tests for correct behavior of CLEAR bookkeeper runs."""

    mode = BookkeeperMode.CLEAR

    def test_run_before_calc_exec(self, before_calc_execution, caplog):
        """Check correct overwriting of input files in CLEAR mode."""
        self.run_before_calc_exec_and_check(before_calc_execution)
        self.check_no_warnings(caplog)

    def test_clear_after_archive(self, after_archive, caplog):
        """Check behavior of CLEAR after ARCHIVE (e.g., manual call)."""
        self.run_after_archive_and_check(after_archive)
        self.check_root_inputs_replaced_by_out(*after_archive)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))

    def test_clear_after_calc_exec(self, after_calc_execution, caplog):
        """Check behavior of CLEAR after a non-ARCHIVEd calc run.

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

        Returns
        -------
        None.
        """
        self.run_and_check_prerun_archiving(after_calc_execution, caplog)
