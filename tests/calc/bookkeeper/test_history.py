"""Tests for module history of viperleed.calc.bookkeeper."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'


import pytest
from pytest_cases import fixture

from viperleed.calc.bookkeeper.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history import _DISCARDED
from viperleed.calc.bookkeeper.history import NoHistoryEntryError

from .conftest import NOTES_TEST_CONTENT


@fixture(name='history_info_file')
def fixture_history_info(after_run):
    bookkeeper, *_ = after_run
    history_info = bookkeeper.history_info
    return history_info, bookkeeper.cwd / HISTORY_INFO_NAME


class TestHistoryInfoFile:

    def test_history_info_read_contents(self,history_info_file):
        """Check that the history.info file is read correctly."""
        history_info, actual_file = history_info_file
        assert actual_file.exists()
        assert history_info.path == actual_file

    def test_history_info_entry_parsing(self,history_info_file):
        history_info, actual_file = history_info_file
        assert bool(history_info.last_entry) == bool(actual_file.read_text())

    def test_history_info_has_notes(self,history_info_file):
        """Check that the history.info file is read correctly."""
        history_info, actual_file = history_info_file
        if NOTES_TEST_CONTENT in actual_file.read_text():
            assert history_info.last_entry.notes == NOTES_TEST_CONTENT

    def test_history_info_discarded(self,history_info_file):
        """Check that the history.info file is read correctly."""
        history_info, actual_file = history_info_file
        discarded_in_entry = _DISCARDED in actual_file.read_text()
        assert history_info.last_entry_was_discarded == discarded_in_entry

    def test_history_info_regenerate_from_entries(self,history_info_file):
        """Check that the history.info file can be regenerated after parsing."""
        history_info, actual_file = history_info_file
        actual_text = actual_file.read_text()
        assert actual_text == history_info.raw_contents
        last_entry = history_info.last_entry
        if last_entry is None:
            return
        history_info.remove_last_entry()
        history_info.append_entry(last_entry)
        assert actual_text.strip() == history_info.raw_contents.strip()

    def test_history_info_discard_last_entry(self, history_info_file):
        """Check we can discard the last history.info entry."""
        history_info, actual_file = history_info_file
        last_entry = history_info.last_entry
        if last_entry is None or last_entry.discarded:
            return
        history_info.discard_last_entry()
        assert history_info.last_entry_was_discarded

    def test_history_info_remove_last_entry(self, history_info_file):
        """Check we can remove the last history.info entry."""
        history_info, *_ = history_info_file
        # check number of entries before and run checks accordingly
        n_entries = history_info.path.read_text().count('# TENSORS')
        if n_entries == 0:
            assert history_info.last_entry is None
            with pytest.raises(NoHistoryEntryError):
                history_info.remove_last_entry()
        else:
            assert history_info.last_entry is not None
            history_info.remove_last_entry()
            assert history_info.path.read_text().count('# TENSORS') == n_entries - 1
