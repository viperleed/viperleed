"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for running bookkeeper in ARCHIVE mode.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

from viperleed.calc.bookkeeper.mode import BookkeeperMode

from ..conftest import MOCK_STATE_FILES
from .run_bookkeeper_base import _TestBookkeeperRunBase


class TestBookkeeperArchive(_TestBookkeeperRunBase):
    """Tests for correct behavior of ARCHIVE bookkeeper runs."""

    mode = BookkeeperMode.ARCHIVE

    def test_archive_after_calc_exec(self, after_calc_execution, caplog):
        """Check correct storage of history files in ARCHIVE mode."""
        self.run_archive_after_calc_and_check(after_calc_execution, caplog)

    def test_archive_again(self, after_archive, caplog):
        """Bookkeeper ARCHIVE after ARCHIVE should not do anything."""
        bookkeeper, *_ = after_archive
        # Write stuff to files to check they are not overwritten
        cwd = bookkeeper.cwd
        sentinel_text = 'something else'
        for file in MOCK_STATE_FILES:
            (cwd / file).write_text(sentinel_text)
        self.run_after_archive_and_check(after_archive,
                                         check_input_contents=False)
        for file in MOCK_STATE_FILES:
            assert (cwd / file).read_text() == sentinel_text
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))

    def test_run_before_calc_exec(self, before_calc_execution, caplog):
        """Check no archiving happens before calc runs."""
        self.run_before_calc_exec_and_check(before_calc_execution)
        self.check_no_warnings(caplog, exclude_msgs=('metadata',))
