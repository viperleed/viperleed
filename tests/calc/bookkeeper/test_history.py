"""Tests for module history of viperleed.calc.bookkeeper."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-02'
__license__ = 'GPLv3+'

from dataclasses import fields
import logging
import re

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases
from pytest_cases.filters import has_tags

from viperleed.calc.bookkeeper.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history import _DISCARDED
from viperleed.calc.bookkeeper.history import _MSG_NOT_UNDERSTOOD_PREFIX
from viperleed.calc.bookkeeper.history import _TAG
from viperleed.calc.bookkeeper.history import CantRemoveEntryError
from viperleed.calc.bookkeeper.history import EntrySyntaxError
from viperleed.calc.bookkeeper.history import HistoryInfoEntry
from viperleed.calc.bookkeeper.history import HistoryInfoError
from viperleed.calc.bookkeeper.history import HistoryInfoFile
from viperleed.calc.bookkeeper.history import NoHistoryEntryError
from viperleed.calc.bookkeeper.history import PureCommentEntry
from viperleed.calc.bookkeeper.history import TimestampFormat

from ...helpers import exclude_tags
from .conftest import NOTES_TEST_CONTENT
from . import cases_bookkeeper
from .cases_bookkeeper import BookkeeperTag as Tag

CorrectEntry = cases_bookkeeper.CasesInfoEntryCorrect


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
    cases=cases_bookkeeper,
    has_tag=Tag.HISTORY,
    filter=exclude_tags(Tag.RAISES),
    )
def fixture_history_info(contents, make_history_file):
    """Return a HistoryInfoFile object."""
    info, *_ = make_history_file(contents)
    info.read()
    return info, contents


class TestHistoryEntry:
    """Collection of tests for history.info entries."""

    auto_fix = parametrize_with_cases('entry_str',
                                      cases=cases_bookkeeper,
                                      has_tag=Tag.AUTO_FIX)
    cant_fix = parametrize_with_cases(
        'entry_str',
        cases=cases_bookkeeper,
        has_tag=Tag.CANT_FIX,
        filter=exclude_tags(Tag.AUTO_FIX,
                            Tag.MULTI_ENTRY,
                            Tag.NEEDS_NO_FIX,
                            Tag.RAISES,
                            Tag.EMPTY)
        )
    need_no_fix = parametrize_with_cases('entry_str',
                                         cases=cases_bookkeeper,
                                         has_tag=Tag.NEEDS_NO_FIX)

    def test_discard_twice(self):
        """Ensure correct behavior of discarding an entry twice."""
        contents = CorrectEntry().case_no_notes()
        entry = HistoryInfoEntry.from_string(contents)
        assert not entry.discarded
        discarded_entry = entry.as_discarded()
        assert discarded_entry.discarded
        assert discarded_entry.as_discarded() is discarded_entry

    @need_no_fix
    def test_fix_not_needed(self, entry_str):
        """Check that correct entries need no fixing."""
        entry = HistoryInfoEntry.from_string(entry_str)
        assert not entry.needs_fixing
        assert entry.as_fixed() is entry
        assert not entry.misses_mandatory_fields

    @auto_fix
    def test_fix_entry(self, entry_str):
        """Check that correct entries need no fixing."""
        entry = HistoryInfoEntry.from_string(entry_str)
        fixed = entry.as_fixed()
        assert entry.needs_fixing
        assert not fixed.needs_fixing
        assert fixed != entry
        assert str(fixed) != str(entry)

    @cant_fix
    def test_cannot_be_fixed(self, entry_str):
        """Check complaints when trying to unfix a non-fixable entry."""
        entry = HistoryInfoEntry.from_string(entry_str)
        assert not entry.was_understood
        assert not entry.can_be_removed
        assert entry.as_fixed() is entry

    @auto_fix
    def test_format_problematic_fields_auto_fix(self, entry_str):
        """Check that no problematic fields are marked."""
        entry = HistoryInfoEntry.from_string(entry_str)
        # pylint: disable-next=magic-value-comparison
        assert 'fixable' in entry.format_problematic_fields()

    @cant_fix
    def test_format_problematic_fields_cant_fix(self, entry_str):
        """Check that no problematic fields are marked."""
        entry = HistoryInfoEntry.from_string(entry_str)
        problems = entry.format_problematic_fields()
        assert any(s in problems for s in ('edited', 'missing'))

    @need_no_fix
    def test_format_problematic_fields_no_fix(self, entry_str):
        """Check that no problematic fields are marked."""
        entry = HistoryInfoEntry.from_string(entry_str)
        assert not entry.format_problematic_fields()

    non_iso_time = parametrize_with_cases(
        'entry_str',
        cases=(
            CorrectEntry.case_german_datetime,
            ),
        )

    @non_iso_time
    def test_to_default_time(self, entry_str):
        """Check correct conversion of TIME format to the default."""
        entry = HistoryInfoEntry.from_string(entry_str)
        entry_default = entry.with_time_format('default')
        assert entry.needs_fixing
        assert not entry_default.needs_fixing
        assert entry_default != entry
        assert entry_default.time_format != entry.time_format
        assert entry_default.time_format is TimestampFormat.DEFAULT
        assert entry_default.timestamp == entry.timestamp
        assert str(entry_default) != str(entry)

    @parametrize(fmt=(f for f in TimestampFormat if f.writable))
    def test_with_time_format(self, fmt):
        """Check conversion to another time format."""
        entry = HistoryInfoEntry.from_string(CorrectEntry().case_no_notes())
        entry_edited = entry.with_time_format(fmt)
        assert entry.timestamp == entry_edited.timestamp
        if fmt is not TimestampFormat.DEFAULT:
            pytest.xfail(reason=('Not sure yet if we should really mark '
                                 'those as such. We currently do not'))
            assert entry_edited.needs_fixing

    empty_field = {
        'job_nums empty': tuple(),
        }
    invalid_field = {
        'folder_name': 1,
        'r_ref': 1+2j,
        'run_info': {1: 1},
        'tensor_nums': {1: 1},
        'tensor_nums negative': (-1,),
        'job_nums': (1.23,),
        'timestamp': {1: 1},
        **empty_field,
        }

    @parametrize('field,value', invalid_field.items(), ids=invalid_field)
    def test_entry_init_wrong_type(self, field, value, caplog, re_match):
        """Check complaints for a HistoryInfoEntry with invalid field type."""
        # We have to use caplog, as HistoryInfoEntry does not raise
        # exceptions for now, although we should eventually.
        dummy_str = CorrectEntry().case_no_notes()
        with caplog.at_level(logging.CRITICAL):
            dummy = HistoryInfoEntry.from_string(dummy_str)
        kwargs = {f.name: getattr(dummy, f.name)
                  for f in fields(dummy)
                  if f.init}
        field_name, *_ = field.split()
        kwargs[field_name] = value
        entry = HistoryInfoEntry(**kwargs)
        assert _MSG_NOT_UNDERSTOOD_PREFIX in caplog.text
        assert not entry.was_understood

        if not field in self.empty_field:
            return
        assert any(re_match(rf'^{_TAG[field_name]}\s*$', line)
                   for line in str(entry).splitlines())

    def test_discarded_entry_with_missing_notes(self, caplog):
        """Check correct formatting of a discarded entry without notes."""
        missing = cases_bookkeeper.case_missing_field('notes', caplog)
        entry = HistoryInfoEntry.from_string(missing)
        discarded = entry.as_discarded()
        assert _TAG['notes'] not in str(entry)
        assert _TAG['notes'] in str(discarded)


class TestHistoryEntryRaises:
    """Tests for HistoryInfoEntry conditions that raise exceptions."""

    @parametrize_with_cases('contents,exc',
                            cases=cases_bookkeeper,
                            has_tag=Tag.RAISES)
    def test_entry_from_string_raises(self, contents, exc):
        """Check complaints when reading a string entry."""
        with pytest.raises(exc):
            HistoryInfoEntry.from_string(contents)

    @parametrize_with_cases('contents,exc',
                            cases=cases_bookkeeper,
                            filter=has_tags(Tag.RAISES, Tag.ENTRY))
    def test_from_string_raises(self, contents, exc):
        """Check complaints when reading a info file with invalid entries."""
        with pytest.raises(exc):
            HistoryInfoEntry.from_string(contents)

    _invalid_time_fmt = {
        '_NOT_A VALID_TIME_FORMAT': ValueError,
        None: TypeError
        }

    @parametrize('fmt,exc', _invalid_time_fmt.items(), ids=_invalid_time_fmt)
    def test_wrong_time_format(self, fmt, exc):
        """Check complaints for an invalid TIME format."""
        entry = HistoryInfoEntry.from_string(CorrectEntry().case_no_notes())
        with pytest.raises(exc):
            entry.with_time_format(fmt)

    @parametrize(pos=(0, 1, 2))
    def test_from_string_unknown_field(self, pos):
        """Check complaints if an unknown filed is found."""
        entry_lines = CorrectEntry().case_notes().splitlines()
        entry_lines.insert(pos, '# UNKOWN    with some value')
        entry_str = '\n'.join(entry_lines)
        with pytest.raises(EntrySyntaxError):
            HistoryInfoEntry.from_string(entry_str)

    def test_cannot_parse(self, monkeypatch):
        """Check complaints if we cannot parse at all a string."""
        entry_str = CorrectEntry().case_notes()
        monkeypatch.setattr('viperleed.calc.bookkeeper.history._ENTRY_RE',
                            re.compile('ThisDoesNotMatchAtAll'))
        with pytest.raises(HistoryInfoError):
            HistoryInfoEntry.from_string(entry_str)


class TestHistoryInfoFile:
    """Collection of tests for reading/writing the history.info file."""

    def test_history_info_read_contents(self, after_calc_run):
        """Check that the history.info file is read correctly."""
        bookkeeper, *_ = after_calc_run
        history_info = bookkeeper.history_info
        actual_file = bookkeeper.cwd / HISTORY_INFO_NAME
        assert actual_file.exists()
        assert history_info.path == actual_file

    def test_history_info_entry_parsing(self, history_info_file):
        """Check correct parsing of a one-entry history.info file."""
        history_info, contents = history_info_file
        assert bool(history_info.last_entry) == bool(contents)

    def test_history_info_has_notes(self, history_info_file):
        """Check that the history.info file is read correctly."""
        history_info, contents = history_info_file
        has_notes = NOTES_TEST_CONTENT in contents
        last_entry = history_info.last_entry
        if last_entry and not isinstance(last_entry, PureCommentEntry):
            assert bool(last_entry.has_notes) == has_notes
        if has_notes:
            assert history_info.last_entry.notes == NOTES_TEST_CONTENT

    def test_history_info_discarded(self, history_info_file):
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

    def test_history_info_regenerate_from_entries(self, history_info_file):
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

    def test_history_info_discard_last_entry(self, history_info_file, caplog):
        """Check we can discard the last history.info entry."""
        history_info, *_ = history_info_file
        last_entry = history_info.last_entry
        if not last_entry:
            with pytest.raises(NoHistoryEntryError):
                history_info.discard_last_entry()
            return
        history_info.discard_last_entry()
        try:
            was_discarded = last_entry.discarded
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
        for tag in _TAG.values():
            n_entries = text.count(tag)
            if n_entries:
                break
        return n_entries

    def test_history_info_remove_last_entry(self, history_info_file):
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


class TestHistoryInfoRaises:
    """Tests for HistoryInfoFile conditions that raise exceptions."""

    def test_filenotfound_at_init(self):
        """Check complaints when initialized with a non-existing file."""
        with pytest.raises(FileNotFoundError):
            HistoryInfoFile('a_file_that_does_not_exist.zz_yy',
                            create_new=False)

    @parametrize_with_cases('contents,exc',
                            cases=cases_bookkeeper,
                            filter=has_tags(Tag.RAISES, Tag.HISTORY))
    def test_file_with_invalid_entry(self, contents, exc, make_history_file):
        """Check complaints when reading a info file with invalid entries."""
        info, *_ = make_history_file(contents)
        with pytest.raises(exc):
            info.read()
        # Disable as we really want to check the private member,
        # since the @property returns a default value of unset
        # pylint: disable-next=protected-access
        assert info._time_format is None
