"""Tests for module history of viperleed.calc.bookkeeper."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize_with_cases
from pytest_cases.filters import has_tags

from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.entry import HistoryInfoEntry
from viperleed.calc.bookkeeper.history.entry.entry import PureCommentEntry
from viperleed.calc.bookkeeper.history.entry.notes_field import _DISCARDED
from viperleed.calc.bookkeeper.history.errors import CantRemoveEntryError
from viperleed.calc.bookkeeper.history.errors import EntrySyntaxError
from viperleed.calc.bookkeeper.history.errors import HistoryInfoError
from viperleed.calc.bookkeeper.history.errors import NoHistoryEntryError
from viperleed.calc.bookkeeper.history.file import HistoryInfoFile

from ....helpers import exclude_tags
from ....helpers import make_obj_raise
from ..conftest import NOTES_TEST_CONTENT
from ..tag import BookkeeperTag as Tag
from . import cases_history
from .entry import cases_entry


# TODO: test that remaking multiple separated entries preserves spacing

all_history_cases = cases_history, cases_entry


@fixture(name='make_history_file', scope='session')
def factory_history_file(tmp_path_factory):
    """Return a HistoryInfoFile with contents."""
    def _make(contents):
        base_path = tmp_path_factory.mktemp(basename='history_info',
                                            numbered=True)
        info = HistoryInfoFile(base_path/HISTORY_INFO_NAME, create_new=True)
        info.path.write_text(contents, encoding='utf-8')
        return info, contents
    return _make


@fixture(name='history_info_file')
@parametrize_with_cases(
    'contents',
    cases=all_history_cases,
    has_tag=Tag.HISTORY,
    filter=exclude_tags(Tag.RAISES),
    )
def fixture_history_info(contents, make_history_file):
    """Return a HistoryInfoFile object."""
    info, *_ = make_history_file(contents)
    info.read()
    return info, contents


class TestHistoryInfoFile:
    """Collection of tests for reading/writing the history.info file."""

    def test_read_contents(self, after_calc_execution):
        """Check that the history.info file is read correctly."""
        bookkeeper, *_ = after_calc_execution
        history_info = bookkeeper.history_info
        actual_file = bookkeeper.cwd / HISTORY_INFO_NAME
        assert actual_file.exists()
        assert history_info.path == actual_file

    def test_entry_parsing(self, history_info_file):
        """Check correct parsing of a one-entry history.info file."""
        history_info, contents = history_info_file
        assert bool(history_info.last_entry) == bool(contents)

    def test_has_notes(self, history_info_file):
        """Check that the history.info file is read correctly."""
        history_info, contents = history_info_file
        has_notes = NOTES_TEST_CONTENT in contents
        last_entry = history_info.last_entry
        if last_entry and not isinstance(last_entry, PureCommentEntry):
            assert bool(last_entry.has_notes) == has_notes
        if has_notes:
            # pylint: disable-next=protected-access       # OK in tests
            notes = history_info.last_entry.notes._get_string_value()
            assert re.fullmatch(rf'{NOTES_TEST_CONTENT}(\n{_DISCARDED})*',
                                notes, re.M)

    def test_discarded(self, history_info_file):
        """Check that the history.info file is read correctly."""
        history_info, contents = history_info_file
        discarded_in_entry = _DISCARDED in contents
        assert history_info.last_entry_was_discarded == discarded_in_entry

    @staticmethod
    def _check_linewise_equal(to_check, expected):
        """Test equality line by line, excluding trailing spaces."""
        to_check = [line.rstrip() for line in to_check.strip().splitlines()]
        expected = [line.rstrip() for line in expected.strip().splitlines()]
        assert to_check == expected

    def test_regenerate_from_entries(self, history_info_file):
        """Check that the history.info file can be regenerated after parsing."""
        history_info, contents = history_info_file
        assert contents == history_info.raw_contents
        last_entry = history_info.last_entry
        try:
            history_info.remove_last_entry()
        except HistoryInfoError:
            assert last_entry is None or not last_entry.can_be_removed
            return
        if history_info.last_entry is last_entry:
            assert isinstance(last_entry, PureCommentEntry)
            return
        history_info.append_entry(last_entry)
        self._check_linewise_equal(contents, history_info.raw_contents)

    def test_discard_last_entry(self, history_info_file, caplog):
        """Check we can discard the last history.info entry."""
        history_info, *_ = history_info_file
        last_entry = history_info.last_entry
        if not last_entry:
            with pytest.raises(NoHistoryEntryError):
                history_info.discard_last_entry()
            return
        history_info.discard_last_entry()
        try:
            was_discarded = last_entry.is_discarded
        except AttributeError:
            return
        assert history_info.last_entry_was_discarded
        if was_discarded:
            # pylint: disable-next=magic-value-comparison
            assert 'already' in caplog.text
            assert history_info.last_entry is last_entry
        else:
            assert history_info.last_entry is not last_entry

    @staticmethod
    def _count_entries(path):
        n_entries = 0
        text = path.read_text()
        for tag in FieldTag:
            if not tag.value:
                continue
            n_entries = text.count(tag.value)
            if n_entries:
                break
        return n_entries

    def test_remove_last_entry(self, history_info_file):
        """Check we can remove the last history.info entry."""
        history_info, *_ = history_info_file
        # check number of entries before and run checks accordingly
        n_entries = self._count_entries(history_info.path)
        last_entry = history_info.last_entry
        if not n_entries and last_entry is None:
            with pytest.raises(NoHistoryEntryError):
                history_info.remove_last_entry()
            with pytest.raises(NoHistoryEntryError):
                history_info.discard_last_entry()
            return
        assert last_entry is not None
        if not last_entry.can_be_removed:
            with pytest.raises(CantRemoveEntryError):
                history_info.remove_last_entry()
            return
        history_info.remove_last_entry()
        n_entries_again = self._count_entries(history_info.path)
        assert n_entries_again == n_entries - 1

    def test_infer_time_format_invalid_entry(self, make_history_file):
        """Check that no time format is detected for invalid file contents."""
        with make_obj_raise(HistoryInfoEntry, EntrySyntaxError, 'from_string'):
            contents = cases_entry.CasesInfoEntryCorrect().case_no_notes()
            info, *_ = make_history_file(contents)
            info.read()
        assert info._time_format is None


class TestHistoryInfoRaises:
    """Tests for HistoryInfoFile conditions that raise exceptions."""

    def test_filenotfound_at_init(self):
        """Check complaints when initialized with a non-existing file."""
        with pytest.raises(FileNotFoundError):
            HistoryInfoFile('a_file_that_does_not_exist.zz_yy',
                            create_new=False)

    @parametrize_with_cases('contents,exc',
                            cases=all_history_cases,
                            filter=has_tags(Tag.RAISES, Tag.HISTORY))
    def test_file_with_invalid_entry(self, contents, exc, make_history_file):
        """Check complaints when reading a info file with invalid entries."""
        info, *_ = make_history_file(contents)
        with pytest.raises(exc):
            info.read()
        # Disable as we really want to check the private member,
        # since the @property returns a default value if unset
        # pylint: disable-next=protected-access
        assert info._time_format is None
