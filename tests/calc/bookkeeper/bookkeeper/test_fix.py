"""Tests for module bookkeeper of viperleed.calc.bookkeeper.

Collects tests for running bookkeeper in FIX mode.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-11-15'
__license__ = 'GPLv3+'

from pytest_cases import fixture

from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.mode import BookkeeperMode

from .run_bookkeeper_base import _TestBookkeeperRunBase


@fixture(name='after_fix')
def fixture_after_discard(after_archive, after_bookkeper_run):
    """Prepare a directory like the one after FIX was executed."""
    return after_bookkeper_run(after_archive, BookkeeperMode.FIX)


class TestBookkeeperFix(_TestBookkeeperRunBase):
    """Tests for correct behavior of FIX bookkeeper runs."""

    mode = BookkeeperMode.FIX

    def test_fix_missing_metadata(self, after_archive):
        """Check correct fixing of missing metadata files in history."""
        bookkeeper, *_ = after_archive
        def _metadata_everywhere():
            # pylint: disable-next=protected-access       # OK in tests
            return all(f.has_metadata for f in bookkeeper.history._subfolders)
        assert not _metadata_everywhere()
        bookkeeper.run(self.mode)
        assert _metadata_everywhere()

        # There was nothing to fix in history.info.
        # Make sure we don't clutter with backups.
        assert not any(bookkeeper.cwd.glob(f'{HISTORY_INFO_NAME}.bak*'))

        # Make sure that nothing was touched
        self.check_root_after_archive(*after_archive)

    def test_fix_twice(self, after_fix, caplog):
        """Check that running FIX twice does nothing."""
        warnings = ()  # No missing-metadata-related warnings!
        self.run_again_and_check_nothing_changed(after_fix, caplog,
                                                 acceptable_warnings=warnings)
