"""Tests for module entry of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-09-17'
__license__ = 'GPLv3+'

from collections import defaultdict
from dataclasses import fields as data_fields
import logging

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases
from pytest_cases.filters import has_tag

from viperleed.calc.bookkeeper.history.constants import HISTORY_INFO_NAME
from viperleed.calc.bookkeeper.history.entry.entry import HistoryInfoEntry
from viperleed.calc.bookkeeper.history.entry.entry import PureCommentEntry
from viperleed.calc.bookkeeper.history.entry.entry import SyntaxErrorLogger
from viperleed.calc.bookkeeper.history.entry.enums import FaultyLabel
from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.field import FieldBase
from viperleed.calc.bookkeeper.history.entry.field import UnknownField
from viperleed.calc.bookkeeper.history.entry.field_collection import FieldList
from viperleed.calc.bookkeeper.history.entry.time_field import TimestampFormat
from viperleed.calc.bookkeeper.history.errors import FieldsScrambledError
from viperleed.calc.bookkeeper.history.errors import FixableSyntaxError
from viperleed.calc.bookkeeper.history.errors import HistoryInfoError

from .....helpers import CustomTestException
from .....helpers import exclude_tags
from .....helpers import not_raises
from ...tag import BookkeeperTag as Tag
from . import cases_entry


# For entries that contain funny stuff
MSG_NOT_UNDERSTOOD_PREFIX = f'{HISTORY_INFO_NAME}: Could not understand '
MSG_FAULTY_PREFIX = f'{HISTORY_INFO_NAME}: Faulty entry is'

CorrectEntry = cases_entry.CasesInfoEntryCorrect


auto_fix = parametrize_with_cases('entry_str',
                                  cases=cases_entry,
                                  has_tag=Tag.AUTO_FIX)
cant_fix = parametrize_with_cases(
    'entry_str',
    cases=cases_entry,
    has_tag=Tag.CANT_FIX,
    filter=exclude_tags(Tag.AUTO_FIX,
                        Tag.MULTI_ENTRY,
                        Tag.NEEDS_NO_FIX,
                        Tag.RAISES,
                        Tag.EMPTY)
    )
need_no_fix = parametrize_with_cases('entry_str',
                                     cases=cases_entry,
                                     has_tag=Tag.NEEDS_NO_FIX)
non_raising = parametrize_with_cases('entry_str',
                                     cases=cases_entry,
                                     filter=exclude_tags(Tag.RAISES))


@fixture(name='make_entry', scope='session')
def factory_make_entry():
    """Return an HistoryInfoEntry from a string or from fields."""
    def _make(as_str=None, **fields):
        if as_str is not None:
            return HistoryInfoEntry.from_string(as_str)
        return HistoryInfoEntry(**fields)
    return _make


class TestHistoryEntry:
    """Collection of tests for history.info entries."""

    def test_attribute_sorting(self, make_entry):
        """Check that the initialization attributes are sorted as expected."""
        entry = make_entry()
        attributes = FieldList(getattr(entry, f.name)
                               for f in data_fields(entry)
                               if f.init)
        with not_raises(FieldsScrambledError):
            attributes.check_sorted()

    @auto_fix
    def test_format_problematic_fields_auto_fix(self, entry_str, make_entry):
        """Check that no problematic fields are marked."""
        entry = make_entry(entry_str)
        formatted = entry.format_problematic_fields()
        assert FaultyLabel.FIXABLE.value in formatted

    @cant_fix
    def test_format_problematic_fields_cant_fix(self, entry_str,
                                                make_entry, subtests):
        """Check that no problematic fields are marked."""
        entry = make_entry(entry_str)
        problems = entry.format_problematic_fields()
        with subtests.test('edited or missing'):
            assert any(label.value in problems
                       for label in (FaultyLabel.EDITED, FaultyLabel.MISSING))

        # Check that all the mandatory fields are present
        mandatory = [t for t in FieldTag if FieldBase.for_tag(t).is_mandatory]
        with subtests.test('all mandatory'):
            assert all(t.value in problems for t in mandatory)

        # Check also that all the extra fields are present
        # pylint: disable-next=protected-access
        extras = [f for f in entry._raw_fields if isinstance(f, UnknownField)]
        with subtests.test('all extra lines'):
            assert all(str(e) in problems for e in extras)

    @need_no_fix
    def test_format_problematic_fields_no_fix(self, entry_str, make_entry):
        """Check that no problematic fields are marked."""
        entry = make_entry(entry_str)
        assert not entry.format_problematic_fields()

    def test_from_string_not_string(self):
        """Check complaints when from_string is called with a non-string."""
        with pytest.raises(TypeError):
            HistoryInfoEntry.from_string(tuple())

    _pure_comment = (
        '# this is just a comment',
        '',
        '    ',
        '\n'*3,
        )

    @parametrize('value', _pure_comment)
    def test_from_string_only_comments(self, value, make_entry):
        """Check correct parsing of a comment-only string."""
        entry = make_entry(value)
        assert isinstance(entry, PureCommentEntry)

    @parametrize(pos=(0, 1, 2))
    def test_from_string_unknown_field(self, pos, caplog, make_entry):
        """Check complaints if an unknown filed is found."""
        entry_lines = CorrectEntry().case_notes().splitlines()
        entry_lines.insert(pos, '# UNKOWN    with some value')
        entry_str = '\n'.join(entry_lines)
        make_entry(entry_str)
        assert FaultyLabel.EXTRA.value in caplog.text

    @parametrize_with_cases('entry_str',
                            cases=cases_entry.CasesInfoEntryCommented)
    def test_has_comments(self, entry_str, make_entry):
        """Check an entry that has no extra comments."""
        entry = make_entry(entry_str)
        assert entry.has_comments

    @need_no_fix
    def test_has_comments_none(self, entry_str, make_entry):
        """Check an entry that has no extra comments."""
        entry = make_entry(entry_str)
        assert not entry.has_comments

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
    # pylint: disable-next=too-many-arguments  # 3 fixtures
    def test_init_wrong_type(self, field, value, make_entry, caplog, re_match):
        """Check complaints for a HistoryInfoEntry with invalid field type."""
        # We have to use caplog, as HistoryInfoEntry does not raise
        # exceptions for now, although we should eventually.
        with caplog.at_level(logging.CRITICAL):
            dummy = make_entry(CorrectEntry().case_no_notes())
        kwargs = {f.name: getattr(dummy, f.name)
                  for f in data_fields(dummy)
                  if f.init}
        field_name, *_ = field.split()
        kwargs[field_name] = value
        entry = make_entry(**kwargs)
        assert MSG_NOT_UNDERSTOOD_PREFIX in caplog.text
        assert not entry.was_understood

        if field not in self.empty_field:
            return
        tag = getattr(dummy, field_name).tag
        assert any(re_match(rf'^{tag.value}\s*$', line)
                   for line in str(entry).splitlines())

    def test_no_raw_fields_to_replace(self, make_entry):
        """Check that, when not read from string, raw fields stay empty."""
        entry = make_entry()
        # pylint: disable-next=protected-access           # OK in tests
        assert not entry._raw_fields
        discarded = entry.as_discarded()   # Copies over the raw fields
        # pylint: disable-next=protected-access           # OK in tests
        assert not discarded._raw_fields

    non_iso_time = parametrize_with_cases(
        'entry_str',
        cases=(
            cases_entry.CasesInfoEntryAutoFixFields.case_german_datetime,
            ),
        )

    @non_iso_time
    def test_to_default_time(self, entry_str, make_entry):
        """Check correct conversion of TIME format to the default."""
        entry = make_entry(entry_str)
        entry_default = entry.with_time_format('default')
        assert entry.needs_fixing
        assert not entry_default.needs_fixing
        assert entry_default != entry
        assert entry_default.time_format != entry.time_format
        assert entry_default.time_format is TimestampFormat.DEFAULT
        assert entry_default.timestamp.value == entry.timestamp.value
        assert str(entry_default) != str(entry)

    _warn_issues = parametrize_with_cases(
        'entry_str',
        cases=cases_entry,
        filter=exclude_tags(Tag.RAISES, Tag.NO_ISSUES, Tag.NEEDS_NO_FIX),
        )

    @_warn_issues
    def test_warns_issues(self, entry_str, make_entry, caplog):
        """Check that entries with problems emit logging warnings."""
        try:
            entry = make_entry(entry_str)
        except ValueError as exc:
            # pylint: disable-next=magic-value-comparison
            assert 'multiple' in str(exc)
            return
        try:
            problems = entry.format_problematic_fields()
        except AttributeError:
            assert isinstance(entry, PureCommentEntry)
            return
        log = caplog.text
        # Only once irrespective of whether there are issues with
        # the entry at the level of the fields or of the whole entry.
        assert log.count(MSG_FAULTY_PREFIX) == 1
        assert problems
        assert log.count(problems) == 1

    _preserved = parametrize_with_cases(
        'entry_str',
        cases=cases_entry,
        filter=exclude_tags(Tag.RAISES, Tag.NOT_PRESERVED, Tag.MULTI_ENTRY),
        )

    @_preserved
    def test_whitespace_preserved(self, entry_str, make_entry):
        """Check that str(entry) does not add extra white space."""
        entry = make_entry(entry_str)
        as_string = str(entry)
        if as_string.startswith('\n') and not entry_str.startswith('\n'):
            # We always add one new-line character at the beginning
            # but it does not matter here. Checks concerning the
            # correctness of the character are done in test_history
            as_string = as_string[1:]
        assert as_string == entry_str

    @_preserved
    def test_from_to_string(self, entry_str, make_entry):                       # TODO: This test is somewhat slow!
        """Ensure that making a entry from str(entry) gives stable results."""
        entry = make_entry(entry_str)
        for _ in range(5):
            as_string = str(entry)
            entry = make_entry(as_string)
        if as_string.startswith('\n') and not entry_str.startswith('\n'):
            # We always add one new-line character at the beginning
            # but it does not matter here. Checks concerning the
            # correctness of the character are done in test_history
            as_string = as_string[1:]
        assert as_string == entry_str

    @parametrize(fmt=(f for f in TimestampFormat if f.writable))
    def test_with_time_format(self, fmt, make_entry):
        """Check conversion to another time format."""
        entry = make_entry(CorrectEntry().case_no_notes())
        entry_edited = entry.with_time_format(fmt)
        assert entry.timestamp.value == entry_edited.timestamp.value
        if fmt is not TimestampFormat.DEFAULT:
            pytest.xfail(reason=('Not sure yet if we should really mark '       # TODO
                                 'those as such. We currently do not'))
            assert entry_edited.needs_fixing

    def test_multiple_notes_no_extras(self, make_entry):
        """Check that an entry with repeated notes has no UnknownField."""
        base = CorrectEntry().case_no_notes()
        base, _ = base.rsplit('\n', 1)
        notes = f'{FieldTag.NOTES} This is a test\n   on multiple lines'
        entry_str = base + '\n' + '\n'.join(notes for _ in range(5))
        entry = make_entry(entry_str)
        assert not any(isinstance(f, UnknownField) for f in entry._raw_fields)


class TestHistoryEntryDiscard:
    """Collection of tests concerning discarding and discarded entries."""

    @parametrize(discard=(True, False))
    def test_discarded_at_init(self, discard, make_entry):
        """Check correct initialization of an entry with discarding info."""
        entry = make_entry(discarded_=discard)
        assert entry.is_discarded == discard
        assert not entry.has_notes  # Despite the .is_discarded state

    def test_discarded_with_missing_notes(self, make_entry):
        """Check correct formatting of a discarded entry without notes."""
        missing = cases_entry.case_missing_field(FieldTag.NOTES)
        entry = make_entry(missing)
        discarded = entry.as_discarded()
        assert FieldTag.NOTES.value not in str(entry)
        assert FieldTag.NOTES.value in str(discarded)
        assert discarded.is_discarded

    def test_discard_twice(self, make_entry):
        """Ensure correct behavior of discarding an entry twice."""
        entry = make_entry(CorrectEntry().case_no_notes())
        assert not entry.is_discarded
        discarded_entry = entry.as_discarded()
        assert discarded_entry.is_discarded
        assert discarded_entry.as_discarded() is discarded_entry


class TestHistoryEntryFix:
    """Collection of tests for fixing various format issues in entries."""

    auto_fix_whole_entry = parametrize_with_cases(
        'entry_str',
        cases=cases_entry,
        has_tag=Tag.AUTO_FIX_ENTRY
        )
    auto_fix_field_only = parametrize_with_cases(
        'entry_str',
        cases=cases_entry,
        has_tag=Tag.AUTO_FIX,
        filter=exclude_tags(Tag.AUTO_FIX_ENTRY)
        )

    @cant_fix
    def test_fix_cannot_be_fixed(self, entry_str, make_entry):
        """Check complaints when trying to fix a non-fixable entry."""
        entry = make_entry(entry_str)
        assert not entry.was_understood
        assert not entry.can_be_removed
        fixed = entry.as_fixed()
        if entry.needs_fixing:
            # Some form of fixing performed at the whole-entry level
            # pylint: disable-next=protected-access       # OK in tests
            assert not fixed._needs_fix_for_fields
            assert fixed is not entry
        else:
            assert fixed is entry

    @auto_fix_whole_entry
    def test_fix_entry_as_a_whole(self, entry_str, make_entry):
        """Check appropriate fixing of an entry with entry-level issues."""
        entry = make_entry(entry_str)
        assert entry.needs_fixing
        # Disable protected-access (W0212) as it's OK in tests
        assert entry._needs_fix_for_entry       # pylint: disable=W0212

        fixed = entry.as_fixed()
        assert not fixed.needs_fixing
        assert fixed != entry
        assert str(fixed) != str(entry)

    @auto_fix
    def test_fix_entry_whole_raises_unfixed_fields(self, entry_str,
                                                   make_entry):
        """Check complaints when fixing entry issues before field issues."""
        entry = make_entry(entry_str)
        assert entry.needs_fixing
        with pytest.raises(HistoryInfoError):
            entry.fix_entry_issues()

    @auto_fix_whole_entry
    def test_fix_failed(self, entry_str, make_entry, monkeypatch):
        """Check complaints when fixing an entry fails."""
        def _do_not_fix(*_):
            pass
        entry = make_entry(entry_str)
        fixers = (attr
                  for attr in dir(entry)
                  if attr.startswith('_do_fix_action_'))
        for fix_method in fixers:
            monkeypatch.setattr(HistoryInfoEntry, fix_method, _do_not_fix)
        assert entry.needs_fixing
        with pytest.raises(HistoryInfoError):
            entry.as_fixed()

    @auto_fix_field_only
    def test_fix_fields_only(self, entry_str, make_entry):
        """Check appropriate fixing of an entry with funny fields."""
        entry = make_entry(entry_str)
        assert entry.needs_fixing
        # Disable protected-access (W0212) as it's OK in tests
        assert entry._needs_fix_for_fields     # pylint: disable=W0212
        assert not entry._needs_fix_for_entry  # pylint: disable=W0212

        fixed = entry.as_fixed()
        assert not fixed.needs_fixing
        assert fixed != entry
        assert str(fixed) != str(entry)

    @need_no_fix
    def test_fix_not_needed(self, entry_str, make_entry):
        """Check that correct entries need no fixing."""
        entry = make_entry(entry_str)
        assert not entry.needs_fixing
        assert entry.as_fixed() is entry
        assert not entry.misses_mandatory_fields


class TestHistoryEntryOutdated:
    """Collection of tests for the .is_only_outdated property."""

    @auto_fix
    def test_auto_fixable(self, entry_str, make_entry, current_cases):
        """Check correct outdated state of fixable entries."""
        case = next(iter(current_cases.values())).func
        is_outdated = has_tag(Tag.OLD)
        entry = make_entry(entry_str)
        assert entry.is_only_outdated == is_outdated(case)

        fixed = entry.as_fixed()
        assert not fixed.is_only_outdated

    @need_no_fix
    def test_no_fix_needed(self, entry_str, make_entry):
        """Check that non-faulty entries are up to date."""
        entry = make_entry(entry_str)
        assert not entry.is_only_outdated

    @cant_fix
    def test_unfixable(self, entry_str, make_entry):
        """Check that faulty entries are considered up to date."""
        entry = make_entry(entry_str)
        assert not entry.is_only_outdated

    @parametrize_with_cases('entry_str', cases=cases_entry.case_missing_field)
    def test_missing(self, entry_str, make_entry):
        """Check that entries with missing fields are not outdated."""
        entry = make_entry(entry_str)
        assert not entry.is_only_outdated


class TestHistoryEntryRaises:
    """Tests for HistoryInfoEntry conditions that raise exceptions."""

    def test_cannot_parse(self, make_entry, monkeypatch):
        """Check complaints if we cannot parse at all a string."""
        entry_str = CorrectEntry().case_notes()
        field_rgx = ('viperleed.calc.bookkeeper.history.entry'
                     '.field.FieldBase.rgx_pattern')
        monkeypatch.setattr(field_rgx, 'ThisDoesNotMatchAtAll')
        with pytest.raises(HistoryInfoError):
            make_entry(entry_str)

    @parametrize_with_cases('contents,exc',
                            cases=cases_entry,
                            has_tag=Tag.RAISES)
    def test_from_string_raises(self, contents, exc, make_entry):
        """Check complaints when reading a string entry."""
        with pytest.raises(exc):
            make_entry(contents)

    def test_wrong_field_type(self):
        """Check complaints when initializing with an unexpected field type."""
        with pytest.raises(TypeError):
            HistoryInfoEntry(UnknownField('abcd'))

    _invalid_time_fmt = {
        '_NOT_A VALID_TIME_FORMAT': ValueError,
        None: TypeError
        }

    @parametrize('fmt,exc', _invalid_time_fmt.items(), ids=_invalid_time_fmt)
    def test_wrong_time_format(self, fmt, exc, make_entry):
        """Check complaints for an invalid TIME format."""
        entry = make_entry(CorrectEntry().case_no_notes())
        with pytest.raises(exc):
            entry.with_time_format(fmt)


class TestPureCommentEntry:
    """Collection of tests for a comment-only entry."""

    def test_cannot_remove(self):
        """Ensure PureCommentEntry can never be removed."""
        entry = PureCommentEntry()
        assert not entry.can_be_removed

    def test_str(self):
        """Check that __str__ produces the expected result."""
        value = 'abcd'
        entry = PureCommentEntry(value)
        assert str(entry) is value


class TestSyntaxErrorLogger:
    """Tests for the SyntaxErrorLogger context manager."""

    @fixture(name='logger', scope='session')
    def fixture_logger(self):
        """Return a SyntaxErrorLogger."""
        def _make(todos=None, for_field=None):
            if todos is None:
                todos = defaultdict(set)
            return SyntaxErrorLogger(todos, for_field=for_field)
        return _make

    def test_logger_unhashable(self, logger):
        """Check correct storage of TODOs with an unhashable field."""
        field = UnknownField({1: 11})
        # Notice that we purposely do not use an IdentitySet,
        # which would be able to handle non-hashable items.
        todos = defaultdict(set)
        # Notice that the order of the contexts is important here.
        # Swapping the two may make the logger see a pytest.Failed
        # exception that is not handled correctly.
        with not_raises(FixableSyntaxError), logger(todos, for_field=field):
            raise FixableSyntaxError('fixable')
        assert todos

    def test_logger_propagates_exception(self, logger):
        """Check that custom exceptions are not swallowed."""
        reason = 'custom'
        with pytest.raises(CustomTestException, match=reason), logger():
            raise CustomTestException(reason)
