"""Tests for module notes_field of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-09-01'
__license__ = 'GPLv3+'

import functools
from dataclasses import dataclass

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
            info.value = value
        if info.str_ is None:
            info.str_ = value.rstrip()
        if not info.str_.startswith(FieldTag.NOTES.value):
            info.str_ = f'{FieldTag.NOTES}{info.str_}'.rstrip()
        return make_notes(value), info

    def make_discarded(self, value, make_notes, info=None, newline=False):
        """Return a discarded NotesField and a NotesCaseInfo."""
        value_str = '' if value is EmptyField else value
        newline = newline and value_str
        discarded = value_str + '\n' if newline else value_str
        notes, info = self.make_notes_and_info(discarded + _DISCARDED,
                                               make_notes, info)
        info.value = value
        info.discarded = True
        return notes, info

    def case_discarded(self, make_notes):
        """Return an empty, discarded NotesField."""
        notes, info = self.case_empty(make_notes)
        info.str_ = None  # Auto-compute
        return self.make_discarded(notes.value, make_notes, info)

    def case_comment(self, make_notes):
        """Return a one-line NotesField with comments in the value."""
        return self.make_notes_and_info('Test note # plus some comment',
                                        make_notes)

    def case_empty(self, make_notes):
        """Return an empty NotesField."""
        info = NotesCaseInfo()
        info.bool_ = False
        info.value = EmptyField
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
        return make(value, make_notes)

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

    _add = {
        'empty extra': NotesField(''),
        'extra': UnknownField('Additional note'),
        'missing extra': NotesField(),
        'notes field': NotesField('Additional note'),
        'string': 'Additional note',
        }

    @all_notes
    @parametrize(extra=_add.values(), ids=_add)
    def test_add_notes(self, notes, info, extra, subtests):
        """Check expected outcome of adding lines to notes."""
        ori_notes = notes
        notes += extra
        with subtests.test('Immutable'):
            assert notes is not ori_notes
        with subtests.test('Preserve discarded'):
            assert notes.is_discarded == ori_notes.is_discarded
            assert notes.is_discarded == info.discarded
        with subtests.test('Value OK'):
            expect = '' if not ori_notes else ori_notes.value
            try:
                extra_str = ('' if extra.is_empty or extra.is_missing
                             else extra.value)
            except AttributeError:  # Already a string
                extra_str = extra
            if expect and extra_str:
                expect += '\n'
            expect += extra_str
            assert notes.value == expect

    @all_notes
    # pylint: disable-next=unused-argument  # info, from decorator
    def test_add_multiple(self, notes, info):
        """Check correct behavior for adding multiple lines at once."""
        new_notes = sum(self._add.values(), notes)
        assert new_notes is not notes
        assert new_notes.is_discarded == notes.is_discarded

    _add_invalid = {
        'not string': 123,
        'unsupported field': FieldBase.from_string('# R REF 0.123'),
        }

    @all_notes
    @parametrize(extra=_add_invalid.values(), ids=_add_invalid)
    # pylint: disable-next=unused-argument  # info, from decorator
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
    def test_as_discarded(self, notes, info):
        """Check outcome of marking a note as discarded."""
        discarded_version = notes.as_discarded()
        assert discarded_version.is_discarded
        assert str(discarded_version).endswith(_DISCARDED)
        assert notes.is_discarded == info.discarded
        if notes.is_discarded:
            assert discarded_version is notes
        else:
            assert discarded_version is not notes

    def test_empty_notes_acceptable(self, make_notes):
        """Check that an empty notes filed is not marked as problematic."""
        notes, *_ = CasesNotesField().case_empty(make_notes)
        assert notes.value is EmptyField
        assert notes.was_understood
        assert not notes.needs_fixing

    @all_notes
    # pylint: disable-next=unused-argument  # info, from decorator
    def test_check_not_empty(self, notes, info):
        """Check there never are complaints, even for empty notes."""
        with not_raises(Exception):
            # pylint: disable-next=protected-access       # OK in tests
            notes._check_not_empty()

    @all_notes
    # pylint: disable-next=unused-argument  # info, from decorator
    def test_update_discarded(self, notes, info):
        """Check correct identification of DISCARDED state."""
        # pylint: disable-next=protected-access           # OK in tests
        value = notes._get_string_value() + _DISCARDED
        set_frozen_attr(notes, 'value', value)
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
