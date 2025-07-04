"""Tests for module info of viperleed.calc.bookkeeper.history."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

from contextlib import nullcontext
import logging
import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_SEPARATOR
from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.entry import HistoryInfoEntry
from viperleed.calc.bookkeeper.history.entry.entry import PureCommentEntry
from viperleed.calc.bookkeeper.history.entry.notes_field import _DISCARDED
from viperleed.calc.bookkeeper.history.errors import CantDiscardEntryError
from viperleed.calc.bookkeeper.history.errors import CantRemoveEntryError
from viperleed.calc.bookkeeper.history.errors import EntrySyntaxError
from viperleed.calc.bookkeeper.history.errors import FixableSyntaxError
from viperleed.calc.bookkeeper.history.errors import FixFailedError
from viperleed.calc.bookkeeper.history.errors import HistoryInfoError
from viperleed.calc.bookkeeper.history.errors import NoHistoryEntryError
from viperleed.calc.bookkeeper.history.info import HistoryInfoFile

from ....helpers import exclude_tags
from ....helpers import has_any_tag
from ....helpers import make_obj_raise
from ....helpers import raises_exception
from ..conftest import NOTES_TEST_CONTENT
from ..tag import BookkeeperTag as Tag
from . import cases_history
from .entry import cases_entry


all_history_cases = cases_history, cases_entry
_MODULE = 'viperleed.calc.bookkeeper.history.info'


@fixture(name='make_history_file', scope='session')
def factory_history_file(tmp_path_factory):
    """Return a HistoryInfoFile with contents."""
    def _make(contents):
        base_path = tmp_path_factory.mktemp(basename='history_info',
                                            numbered=True)
        info = HistoryInfoFile(base_path, create_new=True)
        info.path.write_text(contents, encoding='utf-8')
        return info, contents
    return _make


@fixture(name='simple_info_file')
def fixture_simple_info_file(make_history_file):
    """Return a HistoryInfoFile and file contents with a single entry."""
    sample = cases_entry.CasesInfoEntryCorrect().case_no_notes()
    return make_history_file(sample)


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
        history_info = bookkeeper.history.info
        actual_file = bookkeeper.cwd / HISTORY_INFO_NAME
        assert actual_file.exists()
        assert history_info.path == actual_file

    def test_append_entry(self, history_info_file):
        """Check correct result of appending an entry."""
        entry_str = cases_entry.CasesInfoEntryCorrect().case_no_notes()
        new_entry = HistoryInfoEntry.from_string(entry_str)
        info, contents = history_info_file
        info.append_entry(new_entry, fix_time_format=False)
        *_, last_entry_raw = info.raw_contents.split(HISTORY_INFO_SEPARATOR)
        expect = str(new_entry)
        if not contents:
            expect = expect[1:]  # Skip leading newline
        elif not contents.strip():
            expect = contents + expect
        assert info.last_entry is new_entry
        assert last_entry_raw == expect

    def test_append_entry_fails(self, history_info_file, mocker):
        """Check that failure to append an entry emits log messages."""
        entry_str = cases_entry.CasesInfoEntryCorrect().case_no_notes()
        new_entry = HistoryInfoEntry.from_string(entry_str)
        info, _ = history_info_file
        mocker.patch('pathlib.Path.open', side_effect=OSError)
        mock_log = mocker.patch(f'{_MODULE}.LOGGER.error')
        with pytest.raises(OSError):
            info.append_entry(new_entry, fix_time_format=False)
        mock_log.assert_called_once()
        error_msg, *_ = mock_log.call_args[0]
        assert str(new_entry) in error_msg

    _empty = parametrize_with_cases('contents',
                                    cases=cases_history,
                                    has_tag=Tag.EMPTY)

    @_empty
    def test_empty_file(self, contents, make_history_file):
        """Check that an empty history.info file is handled correctly."""
        if contents is None:
            pytest.skip('Not the contents of a file, but a non-existing one')
        info, *_ = make_history_file(contents)
        # pylint: disable-next=protected-access           # OK in tests
        entries = info._entries
        only_comments = all(isinstance(e, PureCommentEntry)
                            for e in entries)
        assert only_comments if contents else not entries

    def test_entry_parsing(self, history_info_file):
        """Check correct parsing of a one-entry history.info file."""
        history_info, contents = history_info_file
        assert bool(history_info.last_entry) == bool(contents)

    def test_has_notes(self, history_info_file):
        """Check that notes in history.info file entries are read correctly."""
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
        """Check the last_entry_was_discarded property."""
        history_info, contents = history_info_file
        discarded_in_entry = _DISCARDED in contents
        assert history_info.last_entry_was_discarded == discarded_in_entry

    def test_discard_last_entry(self, history_info_file, caplog):
        """Check discarding of the last history.info entry."""
        info, *_ = history_info_file
        was_discarded = info.last_entry_was_discarded
        last_entry = info.last_entry
        _is_comment = isinstance(last_entry, PureCommentEntry)
        exc = (
            NoHistoryEntryError if not last_entry
            else CantDiscardEntryError if _is_comment else None
            )
        context = pytest.raises(exc) if exc else nullcontext()
        with context:
            info.discard_last_entry()
        if exc:
            return
        assert info.last_entry_was_discarded
        if was_discarded:
            # pylint: disable-next=magic-value-comparison
            assert 'already' in caplog.text
            assert info.last_entry is last_entry
        else:
            assert info.last_entry is not last_entry

    def test_init_does_not_read(self, simple_info_file):
        """Check that making a HistoryInfoFile does not read file contents."""
        info, contents = simple_info_file
        assert not info.raw_contents
        assert not info.last_entry
        # pylint: disable-next=protected-access           # OK in tests
        assert not info._time_format
        assert info.raw_contents != contents

        # Now read and check again
        info.read()
        assert info.raw_contents
        assert info.last_entry
        # pylint: disable-next=protected-access           # OK in tests
        assert info._time_format
        assert info.raw_contents == contents

    def test_infer_time_only_once(self, simple_info_file):
        """Check that _infer_time_format only does so once."""
        info, *_ = simple_info_file
        info.read()
        # pylint: disable-next=protected-access           # OK in tests
        format_before = info._time_format
        assert format_before

        # Now erase the file and read again.
        info.path.write_text('', encoding='utf-8')
        info.read()
        # pylint: disable-next=protected-access           # OK in tests
        assert info._time_format is format_before

    @parametrize(exc=(EntrySyntaxError, FixableSyntaxError))
    def test_read_not_raises(self, exc, simple_info_file):
        """Check that .read catches no exceptions."""
        info, *_ = simple_info_file
        with raises_exception(HistoryInfoEntry, exc, 'from_string'):
            # Ensure exception is not caught. HistoryInfoEntry
            # does not normally raise anything: it logs instead.
            # It's a bug to catch it, as it may mean that we do
            # not log an exception. This is for future-proofing
            # of implementation changes.
            info.read()

    def test_regenerate_from_entries(self, history_info_file):
        """Check that a history.info file can be regenerated after parsing."""
        history_info, contents = history_info_file
        assert contents == history_info.raw_contents
        last_entry = history_info.last_entry
        try:
            history_info.remove_last_entry()
        except HistoryInfoError:
            assert last_entry is None or not last_entry.can_be_removed
            return
        history_info.append_entry(last_entry)
        clean_contents = re.sub(rf'{HISTORY_INFO_SEPARATOR}$', '', contents, 1)
        clean_contents = re.sub(rf'{HISTORY_INFO_SEPARATOR.rstrip()}$', '',
                                clean_contents, count=1)
        assert history_info.raw_contents == clean_contents

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
        """Check removal of the last history.info entry, where possible."""
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
        try:
            is_removable = (
                last_entry.can_be_removed
                or (last_entry.is_only_outdated and n_entries == 1)
                or (last_entry.is_only_outdated
                    # pylint: disable-next=protected-access  # OK in test
                    and last_entry.time_format is history_info._time_format)
                )
        except AttributeError:
            assert isinstance(last_entry, PureCommentEntry)
            is_removable = False
        if not is_removable:
            with pytest.raises((CantRemoveEntryError, HistoryInfoError)):
                history_info.remove_last_entry()
            return
        history_info.remove_last_entry()
        n_entries_again = self._count_entries(history_info.path)
        assert n_entries_again == n_entries - 1


class TestHistoryInfoFileFix:
    """Tests for fixing the contents of a history.info file."""

    @staticmethod
    def needs_fix(file):
        """Return whether any entry needs fixing."""
        # pylint: disable-next=protected-access           # OK in tests
        entries = file._entries
        return any(entry.needs_fixing
                   for entry in entries
                   if not isinstance(entry, PureCommentEntry))

    _fixable = parametrize_with_cases(
        'contents',
        cases=all_history_cases,
        filter=has_any_tag(Tag.AUTO_FIX, Tag.AUTO_FIX_ENTRY)
        )

    @_fixable
    def test_fixable(self, contents, make_history_file):
        """Check fixing of a history.info file."""
        info, *_ = make_history_file(contents)
        info.read()
        assert self.needs_fix(info)
        info.fix()
        assert not self.needs_fix(info)

    def test_fix_fails(self, caplog, simple_info_file):
        """Check that failing to fix an entry logs an error message."""
        with caplog.at_level(logging.CRITICAL):
            info, *_ = simple_info_file
            info.read()
        caplog.clear()
        # pylint: disable-next=protected-access           # OK in tests
        entries_before = info._entries.copy()
        with make_obj_raise(HistoryInfoEntry, FixFailedError, 'as_fixed'):
            info.fix()
        assert any(r for r in caplog.records if r.levelno == logging.ERROR)
        assert all(e_after is e_before
                   # pylint: disable-next=protected-access  # OK in tests
                   for e_after, e_before in zip(info._entries, entries_before))

    @parametrize_with_cases('contents', cases=all_history_cases,
                            has_tag=Tag.NEEDS_NO_FIX)
    def test_needs_no_fix(self, contents, make_history_file):
        """Check fixing of stuff that needs no fix."""
        info, *_ = make_history_file(contents)
        info.read()
        assert not self.needs_fix(info)
        # pylint: disable-next=protected-access           # OK in tests
        entries_before = info._entries.copy()
        info.fix()
        assert not self.needs_fix(info)
        assert all(e_after is e_before
                   # pylint: disable-next=protected-access  # OK in tests
                   for e_after, e_before in zip(info._entries, entries_before))
        assert info.raw_contents == contents

    @parametrize_with_cases('contents',
                            cases=cases_entry.CasesInfoEntryPureComment)
    def test_pure_comment(self, contents, make_history_file):
        """Check (no) fixing of comment-only entries."""
        self.test_needs_no_fix(contents, make_history_file)


class TestHistoryInfoRaises:
    """Tests for HistoryInfoFile conditions that raise exceptions."""

    _raise_filenotfound = (
        'fix',
        'read',
        '_do_remove_last_entry',
        )

    @parametrize(method_name=_raise_filenotfound, ids=_raise_filenotfound)
    def test_filenotfound(self, method_name):
        """Check complaints when initialized with a non-existing file."""
        info = HistoryInfoFile('a_file_that_does_not_exist.zz_yy',
                               create_new=False)
        method = getattr(info, method_name)
        with pytest.raises(FileNotFoundError):
            method()

    @parametrize_with_cases('contents,exc',
                            cases=all_history_cases,
                            has_tag=Tag.RAISES)
    def test_file_with_invalid_entry(self, contents, exc, make_history_file):
        """Check complaints when reading an info file with invalid entries."""
        info, *_ = make_history_file(contents)
        with raises_exception(HistoryInfoEntry, exc, 'from_string'):
            info.read()
        # Disable as we really want to check the private member
        # pylint: disable-next=protected-access
        assert info._time_format is None
