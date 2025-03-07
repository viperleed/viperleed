"""Tests for module notes_field of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-09-01'
__license__ = 'GPLv3+'

import functools
from dataclasses import dataclass
import operator

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.field import EmptyField
from viperleed.calc.bookkeeper.history.entry.field import FieldBase
from viperleed.calc.bookkeeper.history.entry.field import MissingField
from viperleed.calc.bookkeeper.history.entry.field import UnknownField
from viperleed.calc.bookkeeper.history.entry.notes_field import NotesField
from viperleed.calc.bookkeeper.history.entry.notes_field import _DISCARDED
from viperleed.calc.lib.dataclass_utils import set_frozen_attr

from .....helpers import not_raises


@fixture(name='make_notes')
def fixture_make_notes(make_field_factory):
    """Return a factory of NotesField instances."""
    return make_field_factory(NotesField)


@dataclass
class NotesCaseInfo:
    """Expected values of various attributes of a note case."""

    value: str = None        # The value, excluding DISCARDED
    discarded: bool = False  # Whether the note is DISCARDED
    n_lines: int = None      # The number of lines in self.value
    bool_: bool = True       # The truth value of the note
    str_ : str = None        # The expected result of str(notes)


class CasesNotesField:
    """Cases for testing NotesField."""

    @staticmethod
    def make_notes_and_info(value, make_notes, info=None):
        """Return a NotesField and a NotesCaseInfo."""
        if info is None:
            info = NotesCaseInfo()
        if info.value is None:
            info.value = tuple(value.splitlines())
        if info.str_ is None:
            info.str_ = value
        if info.n_lines is None:
            info.n_lines = len(info.value)
        if not info.str_.startswith(FieldTag.NOTES.value):
            tag = str(FieldTag.NOTES)
            if not info.str_:
                tag = tag.rstrip()
            info.str_ = f'{tag}{info.str_}'
        return make_notes(value), info

    def make_discarded(self, value, make_notes, info=None, newline=False):
        """Return a discarded NotesField and a NotesCaseInfo."""
        value_str = '' if value is EmptyField else value
        if isinstance(value_str, tuple):
            value_str = '\n'.join(value)
        newline = newline and value_str
        discarded = value_str + '\n' if newline else value_str
        notes, info = self.make_notes_and_info(discarded + _DISCARDED,
                                               make_notes, info)
        if isinstance(value, tuple):
            info.value = value
        elif not value or value is EmptyField:
            info.value = EmptyField
        else:
            info.value = (value,)
        info.discarded = True
        return notes, info

    def case_discarded(self, make_notes):
        """Return an empty, discarded NotesField."""
        _, info = self.case_empty(make_notes)
        info.str_ = None  # Auto-compute
        return self.make_discarded('', make_notes, info)

    @parametrize(pos=range(5))
    def case_discarded_not_at_end(self, pos, make_notes):
        """Return a discarded note where DISCARDED is in the middle."""
        value = tuple(f'line {i}' for i in range(6))
        lines = list(value)
        lines.insert(pos, f'  \t{_DISCARDED}  ')
        info = NotesCaseInfo()
        info.value = value
        info.str_ = '\n'.join(line.rstrip() for line in lines)
        info.discarded = True
        return self.make_notes_and_info('\n'.join(lines), make_notes, info)

    def case_comment(self, make_notes):
        """Return a one-line NotesField with comments in the value."""
        return self.make_notes_and_info('Test note # plus some comment',
                                        make_notes)

    def case_empty(self, make_notes):
        """Return an empty NotesField."""
        info = NotesCaseInfo()
        info.bool_ = False
        info.value = EmptyField
        info.n_lines = 0
        notes, info = self.make_notes_and_info('', make_notes, info)
        return notes, info

    def case_missing(self):
        """Return a missing NotesField."""
        info = NotesCaseInfo()
        info.bool_ = False
        info.value = MissingField
        info.str_ = ''
        return NotesField(), info

    @parametrize(nlines=(2, 3, 5, 7))
    @parametrize(discarded=(True, False))
    def case_multiline(self, nlines, discarded, make_notes):
        """Return a NotesField with multiple lines."""
        value = '\n'.join(f'Line {i+1}' for i in range(nlines))
        make = (self.make_notes_and_info if not discarded
                else functools.partial(self.make_discarded, newline=True))
        notes, info = make(value, make_notes)
        info.value = tuple(value.splitlines())
        info.n_lines = len(info.value)
        return notes, info

    def case_trailing_spaces(self, make_notes):
        """Return a multi-line NotesField with trailing spaces."""
        lines = [' '*2*i + f'Line {i}' + ' '*i for i in range(5)]
        # Add trailing newlines, removed by FieldBase by not NotesField
        lines.extend(' '*i for i in range(5))
        lines_stripped = tuple(line.rstrip() for line in lines)
        notes, info = self.make_notes_and_info('\n'.join(lines), make_notes)
        info.str_ = str(FieldTag.NOTES) + '\n'.join(lines_stripped)
        info.value = lines_stripped
        return notes, info

    def case_one_line(self, make_notes):
        """Return a single-line NotesField."""
        return self.make_notes_and_info('Initial note', make_notes)

    def case_one_line_discarded(self, make_notes):
        """Return a single-line, discarded NotesField."""
        notes, info = self.case_one_line(make_notes)
        info.str_ = None  # Auto-compute again
        return self.make_discarded(notes.value, make_notes,
                                   info=info, newline=True)


all_notes = parametrize_with_cases('notes,info', cases=CasesNotesField)


class TestNotesField:
    """Tests for NotesField instances."""

    @all_notes
    def test_init(self, notes, info):
        """Check simple initialization of notes."""
        assert notes.value == info.value
        assert notes.is_discarded == info.discarded
        assert not notes.has_comments
        assert notes.tag is FieldTag.NOTES

    _add = {  # operand, number of added lines
        'empty notes': (NotesField(''), 1),
        'empty extra': (UnknownField(EmptyField), 1),
        'extra': (UnknownField('Additional note'), 1),
        'missing notes': (NotesField(), 0),
        'notes field': (NotesField('Additional note'), 1),
        'string': ('Additional note', 1),
        'string empty': ('', 1),
        'string multi-line': ('Additional note\n'*3, 4),
        # And the same stuff, but with _DISCARDED (where sensible)
        'empty discarded': (NotesField('').as_discarded(), 0),
        'extra discarded': (UnknownField(_DISCARDED), 0),
        'notes field discarded': (
            NotesField('Additional note').as_discarded(),
            1
            ),
        'string discarded only': (_DISCARDED, 0),
        'string multi-line discarded': ('Additional note\n'*3 + _DISCARDED, 3),
        }

    @all_notes
    @parametrize('extra,n_lines_extra', _add.values(), ids=_add)
    # pylint: disable-next=too-many-arguments        # OK with fixtures
    def test_add(self, notes, info, extra, n_lines_extra, subtests):
        """Check expected outcome of adding lines to notes."""
        ori_notes = notes
        notes += extra
        extra_str = str(extra).replace(notes.tag.value, '').strip()
        with subtests.test('Immutable'):
            assert notes is not ori_notes
        with subtests.test('Preserve discarded'):
            extra_discarded = extra_str.endswith(_DISCARDED)
            notes_discarded = notes.is_discarded
            assert notes_discarded == ori_notes.is_discarded or extra_discarded
            assert notes_discarded == info.discarded or extra_discarded
        with subtests.test('Nr. lines'):
            n_lines_expect = info.n_lines or 0
            n_lines_expect += n_lines_extra
            needs_one_more_line = (
                # Only one line is considered specially for _DISCARDED
                ori_notes.is_discarded and extra_discarded
                or
                # Adding to empty non-discarded gives
                # "Notes:...\n", not "Notes:..."
                ori_notes.is_empty and not ori_notes.is_discarded
                )
            if needs_one_more_line:
                n_lines_expect += 1
            try:
                assert len(notes.value) == n_lines_expect
            except TypeError:
                assert extra.is_missing
                assert notes.is_missing or notes.is_empty
        with subtests.test('Value OK'):
            expect = (tuple() if ori_notes.is_missing
                      else (('',) if not ori_notes else ori_notes.value))
            try:
                extra_str = '' if extra.is_empty else extra.value
            except AttributeError:  # Already a string
                extra_str = extra
            if extra_str is not MissingField:
                expect += ((extra_str,) if isinstance(extra_str, str)
                           else extra_str)  # a tuple
            assert notes.value == expect or MissingField

    @all_notes
    # pylint: disable-next=unused-argument       # info, from decorator
    def test_add_multiple(self, notes, info):
        """Check correct behavior for adding multiple lines at once."""
        non_discarded = (v for v, _ in self._add.values()
                         if not str(v).rstrip().endswith(_DISCARDED))
        # First, only those that are not DISCARDED,
        # and should not change the discarded state
        new_notes = sum(non_discarded, notes)
        assert new_notes is not notes
        assert new_notes.is_discarded == notes.is_discarded

        # Then all of them
        all_added = (v for v, _ in self._add.values())
        new_notes_discarded = sum(all_added, notes)
        assert new_notes_discarded is not notes
        assert new_notes_discarded.is_discarded

    @parametrize_with_cases('notes,_', cases=CasesNotesField)
    @parametrize_with_cases('other,__', cases=CasesNotesField)
    def test_add_not_commutative(self, notes, other, _, __):
        """Check that a+b != b+a."""
        compare = getattr(
            operator,
            'eq' if notes == other or notes.is_missing or other.is_missing
            else 'ne'
            )
        assert compare(notes + other, other + notes)

    _non_notes_add = {k: v
                      for k, (v, _) in _add.items()
                      if not isinstance(v, NotesField)}

    @all_notes
    @parametrize(other=_non_notes_add.values(), ids=_non_notes_add)
    # pylint: disable-next=unused-argument       # info, from decorator
    def test_add_not_commutative_non_notes(self, notes, other, info):
        """Check that a+b != b+a."""
        try:
            other_missing = other.is_missing
        except AttributeError:
            other_missing = False
        try:
            # pylint: disable-next=protected-access       # OK in tests
            other_str = other._get_string_value()
        except AttributeError:  # Already a string
            other_str = other
        equal = (
            notes == other
            or notes.is_missing
            or other_missing
            # pylint: disable-next=protected-access   # OK in tests
            or (notes._get_string_value() == other_str)
            )
        compare = getattr(operator, 'eq' if equal else 'ne')
        assert compare(notes + other, other + notes)

    _add_invalid = {
        'not string': 123,
        'unsupported field': FieldBase.from_string('# R REF 0.123'),
        }

    @all_notes
    @parametrize(extra=_add_invalid.values(), ids=_add_invalid)
    # pylint: disable-next=unused-argument       # info, from decorator
    def test_add_raises(self, notes, info, extra):
        """Check complaints when trying to add an invalid type."""
        with pytest.raises(TypeError):
            notes += extra

    @all_notes
    def test_bool(self, notes, info):
        """Check expected result of bool(notes)."""
        assert bool(notes) == info.bool_

    @all_notes
    def test_str(self, notes, info):
        """Check the result of str(notes)."""
        assert str(notes) == info.str_

    @all_notes
    # pylint: disable-next=unused-argument       # info, from decorator
    def test_str_preserved(self, notes, info):
        """Check that making notes from str(notes) preserves str."""
        as_str = str(notes)
        notes_again = NotesField.from_string(as_str)
        assert str(notes_again) == as_str

    @all_notes
    def test_as_discarded(self, notes, info):
        """Check outcome of marking a note as discarded."""
        discarded_version = notes.as_discarded()
        str_discarded = str(discarded_version)
        assert notes.is_discarded == info.discarded
        assert discarded_version.is_discarded
        assert any(line.endswith(_DISCARDED)
                   for line in str_discarded.splitlines())
        if notes.is_discarded:
            assert discarded_version is notes
        else:
            assert discarded_version is not notes

    @all_notes
    # pylint: disable-next=unused-argument  # info, from decorator
    def test_check_not_empty(self, notes, info):
        """Check there never are complaints, even for empty notes."""
        with not_raises(Exception):
            # pylint: disable-next=protected-access       # OK in tests
            notes._check_not_empty()

    @parametrize_with_cases(
        'notes, _',
        cases=(CasesNotesField.case_empty, CasesNotesField.case_discarded),
        )
    def test_discarded_empty(self, notes, _):
        """Check string version of a discarded empty note."""
        expect_str = 'Notes: DISCARDED'
        assert str(notes.as_discarded()) == expect_str

        # Now add some notes and check again
        extra = 'Some text'
        notes += extra
        # expect_str = f'Notes: {extra}\nDISCARDED\n'
        lines = (f'{_DISCARDED}' if notes.is_discarded else f'\n{extra}',
                 f'{extra}' if notes.is_discarded else f'{_DISCARDED}')
        expect_str = 'Notes: ' + '\n'.join(lines)
        expect_value = (() if notes.is_discarded else ('',)) + (extra,)
        discarded = notes.as_discarded()
        assert notes.value == expect_value
        assert discarded.value == expect_value
        assert str(discarded) == expect_str

    _discard_midway = {
        'one line, start': f'  {_DISCARDED} other text',
        'one line, midway': f'some {_DISCARDED} text',
        'one line, end': f'this line was {_DISCARDED}',
        'more lines': f'first\nsecond\nthis is a {_DISCARDED} one\nlast',
        }

    @parametrize(value=_discard_midway.values(), ids=_discard_midway)
    def test_discarded_midway(self, value, make_notes):
        """Check correct identification of DISCARDED somewhere in a line."""
        notes = make_notes(value)
        assert not notes.is_discarded
        assert notes.needs_fixing

        notes_fixed = notes.as_fixed()
        with not_raises(Exception):
            notes_fixed.check_value()
        assert notes_fixed.is_discarded
        assert not notes_fixed.needs_fixing

        n_lines_ori = len(notes.value)
        n_lines_fix = len(notes_fixed.value)
        assert n_lines_fix - n_lines_ori in {0, 1}

    empty_notes = parametrize_with_cases('notes, _',
                                         cases=CasesNotesField.case_empty)

    @empty_notes
    def test_empty_notes_acceptable(self, notes, _):
        """Check that an empty notes filed is not marked as problematic."""
        assert notes.value is EmptyField
        assert notes.was_understood
        assert not notes.needs_fixing

    _invalid = {
        'int': 1,
        'dict': {1: 11},
        'dict with string keys': {'1': 11},
        'tuple of non-string': ('1', 2, '3'),
        }

    @parametrize(value=_invalid.values(), ids=_invalid)
    def test_find_discarded_line_raises(self, value, make_notes):
        """Check complaints for invalid notes values."""
        notes = make_notes(value)
        with pytest.raises(TypeError):
            # pylint: disable-next=protected-access       # OK in tests
            notes._find_discarded_line()

    @all_notes
    # pylint: disable-next=unused-argument  # info, from decorator
    def test_update_discarded(self, notes, info):
        """Check correct identification of DISCARDED state."""
        # pylint: disable-next=protected-access           # OK in tests
        value = notes._get_string_value()
        if value:
            value += '\n'
        value += _DISCARDED
        set_frozen_attr(notes, 'value', value)
        # pylint: disable-next=protected-access           # OK in tests
        notes._cleanup_str_value()  # Store as tuple
        # pylint: disable-next=protected-access           # OK in tests
        notes._update_discarded()
        assert notes.is_discarded


class TestNotesFieldClassVar:
    """Tests for the class-level attributes of NotesField."""

    def test_mandatory_flag(self):
        """Check that notes are mandatory."""
        assert NotesField.is_mandatory

    def test_tag(self):
        """Check the correct tag."""
        assert NotesField.tag is FieldTag.NOTES
