"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for simulating the default execution of bookkeeper
before and after viperleed.calc.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-11-15'
__license__ = 'GPLv3+'

from viperleed.calc.bookkeeper.bookkeeper import Bookkeeper
from viperleed.calc.constants import DEFAULT_HISTORY

from ..conftest import MOCK_TIMESTAMP
from .run_bookkeeper_base import _TestBookkeeperRunBase


class TestBookkeeperDuringCalc(_TestBookkeeperRunBase):
    """Test the same conditions as when bookkeeper runs around calc."""

    def test_run_around_calc(self,
                             mock_tree_before_calc_execution,
                             mock_tree_after_calc_execution,
                             caplog):
        """Check reuse of bookkeeper in the default calls around calc."""
        tmp_path = mock_tree_before_calc_execution()
        bookkeeper = Bookkeeper(cwd=tmp_path)
        # Before calc, we run in CLEAR mode
        self.run_before_calc_exec_and_check(bookkeeper,
                                            caplog,
                                            check_archiving_required=False,
                                            mode='clear')
        # Then, we simulate a calc run that produces some output
        mock_tree_after_calc_execution()
        # After calc, we run in ARCHIVE mode
        after_calc = (
            bookkeeper,
            bookkeeper.cwd / DEFAULT_HISTORY / f't004.r001_{MOCK_TIMESTAMP}',
            )
        self.run_archive_after_calc_and_check(after_calc,
                                              caplog,
                                              check_archiving_required=False)
