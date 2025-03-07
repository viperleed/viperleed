"""Test field_collection module of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-09-01'
__license__ = 'GPLv3+'

from random import shuffle

import pytest
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.bookkeeper.history.entry.enums import DuplicateType
from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.field import FieldBase
from viperleed.calc.bookkeeper.history.entry.field import UnknownField
from viperleed.calc.bookkeeper.history.entry.field_collection import FieldList
from viperleed.calc.bookkeeper.history.entry.list_of_int_field import (
    JobIdsField,
    RunInfoField,
    TensorNumsField,
    )
from viperleed.calc.bookkeeper.history.entry.notes_field import NotesField
from viperleed.calc.bookkeeper.history.entry.rfactor_field import RRefField
from viperleed.calc.bookkeeper.history.entry.rfactor_field import RSuperField
from viperleed.calc.bookkeeper.history.entry.string_field import FolderField
from viperleed.calc.bookkeeper.history.entry.string_field import JobNameField
from viperleed.calc.bookkeeper.history.entry.time_field import TimestampField
from viperleed.calc.bookkeeper.history.errors import FieldsScrambledError

from .....helpers import not_raises


def shuffled(sequence):
    """Return an in-place shuffled version of sequence."""
    shuffle(sequence)
    return sequence


def twice(obj):
    """Return obj twice."""
    return obj, obj


MockField = type('MockField', (), {'tag': 1})


class CasesFieldList:
    """Test cases for FieldList."""

    @staticmethod
    def _make_field_list(*values):
        """Return a FieldList and the values used to make it."""
        return FieldList(*values), values

    def case_empty(self):
        """Return an empty FieldList."""
        return self._make_field_list()

    def case_one_field(self):
        """Return a FieldList with a single item."""
        return self._make_field_list(UnknownField('abc'))

    def case_no_duplicates(self):
        """Return a FieldList with multiple fields and no duplicates."""
        return self._make_field_list(
            UnknownField('abc'),
            UnknownField('def'),
            RRefField('0.123'),
            TensorNumsField('None'),
            )

    def case_duplicates(self):
        """Return a FieldList with two identical fields."""
        return self._make_field_list(
            UnknownField('abc'),
            UnknownField('def'),
            RRefField('0.123'),
            TensorNumsField('None'),
            TensorNumsField('None'),
            )

    def case_non_hashable(self):
        """Return a FieldList with a non-hashable field."""
        return self._make_field_list(UnknownField({}))


all_fields = parametrize_with_cases('case_', CasesFieldList)


class TestFieldListSequence:
    """Tests for the methods that make FieldList class a MutableSequence."""

    @all_fields
    def test_init(self, case_):
        """Check initialization of FieldList."""
        fields, values = case_
        assert len(fields) == len(values)
        assert tuple(fields) == values
        assert bool(fields) == bool(values)

    def test_init_generator(self):
        """Check initialization from a single generator."""
        nitems = 2
        fields = (UnknownField('abc') for _ in range(nitems))
        field_list = FieldList(fields)
        assert len(field_list) == nitems

    def test_init_invalid(self):
        """Check initialization with invalid item types."""
        with pytest.raises(TypeError, match='Can only contain'):
            FieldList('123')

    @all_fields
    def test_contains(self, case_):
        """Check result of __contains__."""
        fields, values = case_
        try:
            first = values[0]
        except IndexError:  # Empty
            return
        assert first in fields
        assert FieldBase.for_tag(first.tag)(first) in fields

    _not_there = (
        UnknownField('not present'),
        {},
        'abc',
        )

    @all_fields
    @parametrize(not_there=_not_there)
    def test_contains_false(self, case_, not_there):
        """Check that a different object is not in a FieldList."""
        fields, *_ = case_
        assert not_there not in fields

    @all_fields
    def test_delitem_first(self, case_):
        """Check deletion of the first item."""
        fields, values = case_
        if not fields:
            return
        with not_raises(Exception):
            del fields[0]
        assert fields == values[1:]

    @all_fields
    def test_delitem_rest(self, case_):
        """Check deletion of all but the first item."""
        fields, values = case_
        with not_raises(Exception):
            del fields[1::]
        assert fields == values[:1]

    _eq = {
        'empty FieldList': ((), FieldList(), True),
        'empty list': ((), [], True),
        'empty tuple': ((), (), True),
        'empty, different list': ((), [1, 2], False),
        'empty, different tuple': ((), ('1', '2'), False),
        'not empty, FieldList': (
            (TensorNumsField(), UnknownField(), NotesField()),
            FieldList(TensorNumsField(), UnknownField(), NotesField()),
            True,
            ),
        'duplicates': (
            *twice((TensorNumsField(), UnknownField(),
                    NotesField(), NotesField())),
            True,
            ),
        'different': (
            (UnknownField(), TensorNumsField(), TensorNumsField()),
            (UnknownField(), TensorNumsField(), UnknownField()),
            False,
            ),
        'wrong type': ((UnknownField(1),), 1, False),
        }

    @parametrize('values,other,expect', _eq.values(), ids=_eq)
    def test_equal(self, values, other, expect):
        """Check equality of FieldList(*values) and `other`."""
        fields = FieldList(*values)
        equal = fields == other
        assert equal is expect

    @all_fields
    def test_getitem(self, case_):
        """Check correct item access."""
        fields, values = case_
        if fields:
            assert fields[0] == values[0]
        with not_raises(Exception):
            _ = fields[1:10:2]
        with pytest.raises(IndexError):
            _ = fields[len(fields)]

    _index = {
        'extra': (UnknownField('2'), 1),
        'standard': (NotesField('3'), 2),
        'duplicate': (UnknownField('1'), 0),
        'non-hashable': (UnknownField({4}), 4),
        }

    @parametrize('value,expect', _index.values(), ids=_index)
    def test_index(self, value, expect):
        """Check expected result of indexing a field."""
        fields = FieldList(
            UnknownField('1'),
            UnknownField('2'),
            NotesField('3'),
            UnknownField('2'),
            UnknownField({4}),  # not hashable
            )
        result = fields.index(value)
        assert result == expect

    _index_raises = {
        'not present': NotesField('this one is not there'),
        'not hashable': UnknownField({1}),
        }

    @all_fields
    @parametrize(find=_index_raises.values(), ids=_index_raises)
    def test_index_raises(self, case_, find):
        """Check complaints when trying to index a non-existent field."""
        fields, *_ = case_
        with pytest.raises(ValueError):
            fields.index(find)

    _index_start_stop = {
        'start only': (3,),
        'start and stop': (1, -5),
        }

    @all_fields
    @parametrize(start_stop=_index_start_stop.values(), ids=_index_start_stop)
    def test_index_start_stop(self, case_, start_stop):
        """Check indexing when start or stop values are given."""
        fields, *_ = case_
        with pytest.raises(ValueError, match='not in list'):
            fields.index(NotesField('this is not there'), *start_stop)

    @all_fields
    def test_insert(self, case_):
        """Check correct item insertion."""
        fields, *_ = case_
        insert = UnknownField()
        fields.insert(0, insert)
        assert fields[0] is insert

    def test_insert_raises(self):
        """Check complaints when inserting a non-field."""
        fields = FieldList()
        with pytest.raises(TypeError, match='Can only contain'):
            fields.insert(0, 1.35)

    @all_fields
    def test_setitem(self, case_):
        """Check correct behavior of setting new items."""
        fields, *_ = case_
        set_one, *set_three = (UnknownField(f'extra_{i+1}') for i in range(4))
        if fields:
            fields[-1] = set_one
            assert fields[-1] is set_one
        fields[0:3] = set_three
        assert fields[0:3] == set_three

    _set_invalid = {
        'one': (0, 'set'),
        'two': (slice(0, 2), ('set', UnknownField('set'))),
        'generator': (slice(0, 3), ('set' for _ in range(3))),
        }

    @parametrize('index,value', _set_invalid.values(), ids=_set_invalid)
    def test_setitem_raises(self, index, value):
        """Check complaints when setting non-field items."""
        initial = (UnknownField('test'),)
        fields = FieldList(*initial)
        with pytest.raises(TypeError, match='Can only contain'):
            fields[index] = value
        assert fields == initial

    _repr = {
        'empty': ((), 'FieldList()'),
        'items': ((UnknownField('1'), TensorNumsField('2')),
                  '''\
FieldList(
    UnknownField(value='1'),
    TensorNumsField(value='2'),
    )'''),
        }

    @parametrize('values,expect', _repr.values(), ids=_repr)
    def test_repr(self, values, expect):
        """Check expected outcome of __repr__."""
        fields = FieldList(*values)
        assert repr(fields) == expect
        assert str(fields) == expect


class TestFieldListMethods:
    """Tests for non-MutableSequence methods of FieldList."""

    _sorted = {
        'empty': (),
        'no unknown': (
            TensorNumsField('1, 2'),
            RRefField('1.23'),
            ),
        'all': tuple(FieldBase.for_tag(t)() for t in FieldTag),
        'with extras': (
            UnknownField('one'),
            TensorNumsField('1, 2'),
            UnknownField('two'),
            UnknownField('three'),
            RRefField('1.23'),
            UnknownField('four'),
            ),
        }

    @parametrize(values=_sorted.values(), ids=_sorted)
    def test_check_sorted(self, values):
        """Check correct identification of sort order."""
        fields = FieldList(*values)
        with not_raises(FieldsScrambledError):
            fields.check_sorted()

    def test_check_sorted_raises(self):
        """Check complaints for non-sorted fields."""
        fields = FieldList(
            UnknownField('one'),
            RRefField('1.23'),
            UnknownField('two'),
            UnknownField('three'),
            TensorNumsField('1, 2'),
            UnknownField('four'),
            )
        with pytest.raises(FieldsScrambledError):
            fields.check_sorted()

    @all_fields
    def test_copy(self, case_):
        """Check correct duplication of a FieldList."""
        fields, *_ = case_
        fields_copy = fields.copy()
        assert isinstance(fields_copy, FieldList)
        assert fields == fields_copy
        assert fields is not fields_copy
        assert all(f_copy is f_ori  # Copy is shallow
                   for f_copy, f_ori in zip(fields_copy, fields))

    _duplicates = {
        'empty': ((), FieldBase, (DuplicateType.NONE, ())),
        'not present': ([UnknownField('1'), UnknownField('1')],
                        NotesField('test'),
                        (DuplicateType.NONE, ())),
        'not duplicate': ([UnknownField('1'), NotesField()], UnknownField,
                          (DuplicateType.NONE, ())),
        'different': (
            [UnknownField('1'), UnknownField('2')], UnknownField,
            (DuplicateType.DIFFERENT, (UnknownField('1'), UnknownField('2'))),
            ),
        'identical': (
            [UnknownField('1'), UnknownField('1')], UnknownField,
            (DuplicateType.IDENTICAL, (UnknownField('1'), UnknownField('1'))),
            ),
        'notes identical': (
            [NotesField() for _ in range(5)], NotesField,
            (DuplicateType.IDENTICAL, tuple(NotesField() for _ in range(5))),
            ),
        'notes different': (
            [NotesField(str(i)) for i in range(5)], NotesField,
            (DuplicateType.NOTES, tuple(NotesField(str(i)) for i in range(5))),
            ),
        }

    @parametrize('values,field,expect', _duplicates.values(), ids=_duplicates)
    def test_duplicates(self, values, field, expect):
        """Check correct identification of duplicates."""
        fields = FieldList(*values)
        result = fields.duplicates(field)
        assert result == expect

    _type = {
        'empty': ((), NotesField, ()),
        'not there': ((UnknownField(),), NotesField('123'), ()),
        'one': ((UnknownField('1'),), UnknownField, (UnknownField('1'),)),
        'multiple': (
            (UnknownField('1'), UnknownField('2'), NotesField('123')),
            UnknownField, (UnknownField('1'), UnknownField('2'))
            ),
        'hashable non field': ((), '123', ()),
        }

    def test_find_exact_index_raises(self):
        """Check complaints for index look-up of a non-existing value."""
        fields = FieldList()
        with pytest.raises(ValueError):
            # pylint: disable-next=protected-access       # OK in tests
            fields._find_exact_index(UnknownField())

    @parametrize('values,field,selected', _type.values(), ids=_type)
    def test_first_by_type(self, values, field, selected):
        """Check selection of the first filed with a given type."""
        try:
            expect = selected[0]
        except IndexError:
            expect = None
        fields = FieldList(*values)
        assert fields.first_by_type(field) == expect

    def test_has_only_pure_comments(self):
        """Check the has_only_pure_comments property."""
        unknown_field = UnknownField()
        fields = FieldList(unknown_field)
        assert fields.has_only_pure_comments
        fields.append(NotesField(''))
        assert not fields.has_only_pure_comments

    # Here some initialization items for testing insert_sorted:
    _insert_sort_items = {
        'unsorted': (TensorNumsField(),
                     UnknownField(),
                     NotesField(),
                     RRefField()),  # Purposely at the wrong position
        'unknown first': (UnknownField(),
                          RRefField()),
        'only optional': (UnknownField(),
                          UnknownField(),
                          JobNameField(),
                          RunInfoField(),
                          UnknownField(),
                          RRefField(),
                          RSuperField(),
                          UnknownField()),
        'unknown only': (UnknownField(),
                         UnknownField()),
        'no unknown': (TensorNumsField(),
                       JobIdsField(),
                       JobIdsField(),
                       JobIdsField(),
                       JobIdsField(),
                       JobNameField()),
        }
    # And the items to insert
    _insert_sort = {  # items, new field, expected insertion index
        'same type, right order': ('unsorted', TensorNumsField(''), 1),
        'new type, last': ('unsorted', RSuperField(''), 4),
        'same type, wrong order': ('unsorted', RRefField(''), 4),
        'new type, in between': ('unsorted', JobIdsField(''), 2),
        'insert missing tensors': ('unknown first', TensorNumsField(), 1),
        'insert missing notes': ('unknown first', NotesField(), 2),
        'insert missing folder': ('only optional', FolderField(), 8),
        'insert missing time': ('only optional', TimestampField(), 5),
        'insert notes at end': ('only optional', NotesField(), 8),
        'insert after unknown': ('unknown only', TensorNumsField(), 2),
        'insert unknown': ('unknown first', UnknownField(), 1),
        'insert unknown into unknown': ('unknown only', UnknownField(), 2),
        'insert unknown at end': ('no unknown', UnknownField(), 6),
        'insert same value': ('no unknown', JobIdsField(), 5),
        }

    @parametrize('items_key,new_field,expect',
                 _insert_sort.values(),
                 ids=_insert_sort)
    def test_insert_sorted(self, items_key, new_field, expect):
        """Check that insertion of a `new_field` works as expected."""
        fields = FieldList(*self._insert_sort_items[items_key])
        fields.insert_sorted(new_field)
        # pylint: disable-next=protected-access           # OK in tests
        new_ind = fields._find_exact_index(new_field)
        assert new_ind == expect
        # pylint: disable-next=magic-value-comparison
        if items_key != 'unsorted':
            with not_raises(FieldsScrambledError):
                fields.check_sorted()

    _insert_indices = {
        # We can't really use all permutations, as they would be
        # 9! == 362880 tests. Too many. Pick a few representative
        # ones, plus have the test permute them again every time
        # so we may spot problematic ones.
        (0, 1, 2, 3, 4, 5, 6, 7, 8),
        (0, 3, 4, 1, 7, 8, 6, 2, 5),
        (0, 3, 5, 2, 8, 1, 7, 6, 4),
        (0, 6, 5, 1, 8, 4, 7, 3, 2),
        (0, 7, 8, 1, 4, 5, 3, 6, 2),
        (1, 0, 2, 5, 6, 7, 4, 3, 8),
        (1, 0, 5, 8, 4, 3, 2, 6, 7),
        (1, 2, 5, 8, 4, 6, 0, 7, 3),
        (1, 4, 2, 3, 6, 0, 5, 8, 7),
        (1, 4, 2, 8, 3, 7, 0, 6, 5),
        (1, 5, 3, 6, 8, 7, 0, 2, 4),
        (2, 1, 4, 6, 0, 5, 7, 8, 3),
        (2, 5, 8, 1, 4, 7, 3, 0, 6),
        (3, 1, 4, 2, 0, 5, 6, 7, 8),
        (3, 4, 7, 0, 8, 6, 2, 5, 1),
        (3, 5, 2, 8, 1, 4, 6, 0, 7),
        (3, 5, 4, 0, 7, 8, 1, 2, 6),
        (4, 0, 5, 7, 8, 6, 3, 1, 2),
        (4, 1, 3, 0, 5, 7, 2, 6, 8),
        (4, 2, 3, 1, 8, 6, 7, 0, 5),
        (4, 2, 5, 7, 1, 8, 0, 6, 3),
        (4, 3, 6, 8, 0, 5, 2, 1, 7),
        (4, 5, 0, 6, 2, 8, 1, 3, 7),
        (4, 6, 2, 5, 0, 8, 3, 1, 7),
        (4, 6, 2, 7, 5, 1, 8, 0, 3),
        (4, 6, 8, 3, 7, 1, 5, 2, 0),
        (5, 1, 8, 4, 6, 0, 7, 2, 3),
        (5, 2, 3, 6, 0, 8, 7, 4, 1),
        (5, 3, 8, 4, 6, 1, 7, 0, 2),
        (5, 6, 1, 2, 8, 7, 3, 4, 0),
        (5, 6, 7, 8, 3, 4, 0, 2, 1),
        (5, 7, 3, 6, 0, 1, 4, 8, 2),
        (6, 1, 7, 0, 2, 4, 3, 5, 8),
        (6, 4, 8, 0, 2, 3, 1, 5, 7),
        (7, 0, 2, 5, 6, 3, 8, 4, 1),
        (7, 8, 3, 5, 1, 2, 4, 0, 6),
        (8, 0, 3, 1, 7, 4, 2, 6, 5),
        (8, 0, 6, 3, 1, 7, 4, 5, 2),
        (0, 7, 2, 6, 5, 4, 3, 8, 1),
        (2, 0, 7, 1, 6, 3, 4, 5, 8),
        (2, 1, 7, 6, 5, 4, 8, 3, 0),
        (2, 6, 7, 4, 1, 0, 3, 8, 5),
        (3, 1, 0, 7, 8, 5, 2, 6, 4),
        (5, 3, 2, 6, 7, 8, 4, 0, 1),
        (5, 6, 7, 1, 0, 8, 2, 4, 3),
        (5, 7, 3, 1, 0, 2, 8, 4, 6),
        (6, 0, 5, 4, 2, 8, 7, 1, 3),
        (6, 2, 1, 7, 3, 5, 4, 0, 8),
        (6, 2, 5, 8, 3, 1, 0, 7, 4),
        (6, 4, 8, 2, 1, 5, 3, 7, 0),
        (6, 4, 8, 7, 0, 1, 2, 3, 5),
        (7, 1, 2, 4, 5, 3, 8, 6, 0),
        (7, 5, 0, 2, 4, 1, 8, 6, 3),
        (7, 5, 2, 3, 8, 0, 1, 6, 4),
        (7, 6, 5, 4, 3, 0, 1, 2, 8),
        (7, 8, 5, 4, 6, 1, 0, 2, 3),
        (8, 2, 0, 5, 6, 7, 4, 1, 3),
        (8, 5, 7, 4, 2, 6, 0, 1, 3),
        (8, 7, 6, 0, 3, 2, 4, 5, 1),
        }

    def _check_insert_one_at_a_time(self, items, indices):
        """Check that insertion of fields one at a time is sorted."""
        fields = FieldList()
        sorted_items = [items[i] for i in indices]
        for item in sorted_items:
            fields.insert_sorted(item)
        try:
            fields.check_sorted()
        except FieldsScrambledError:
            reason = (f'Wrong sorting with {indices}.\nInserted\n   '
                      + ',\n   '.join(type(f).__name__ for f in sorted_items)
                      + f'\nAnd got {fields}')
            pytest.fail(reason)

    @parametrize(indices=_insert_indices)
    def test_insert_sorted_one_at_a_time(self, indices):
        """Check that unsorted insertion of items gives a sorted list."""
        items = [FieldBase.for_tag(t)() for t in FieldTag
                 if t is not FieldTag.UNKNOWN]
        if len(items) != len(indices):
            raise ValueError(
                f'Indices {indices} are out of date. Probably the number of '
                'FieldTag items has changed. Please update indices so that '
                f'len(indices) == {len(items)}'
                )
        indices = list(indices)
        self._check_insert_one_at_a_time(items, indices)
        for _ in range(100):  # Try 100 times to generate a new one
            shuffle(indices)
            if tuple(indices) not in self._insert_indices:
                self._check_insert_one_at_a_time(items, indices)
                return

    _insert_sort_raises = {
        'has no .tag': '123',
        'not a field': MockField(),
        }

    @parametrize(field=_insert_sort_raises.values(), ids=_insert_sort_raises)
    def test_insert_sorted_raises(self, field):
        """Check complaints when trying to insert an invalid field."""
        fields = FieldList()
        with pytest.raises(TypeError, match='is not a field'):
            fields.insert_sorted(field)

    @parametrize('values,field,expect', _type.values(), ids=_type)
    def test_select_type(self, values, field, expect):
        """Check outcome for selecting only fields of a certain type."""
        fields = FieldList(*values)
        selected = fields.select_type(field)
        assert selected == expect

    def test_select_type_raises(self):
        """Check complaints when trying to select a non-field."""
        fields = FieldList()
        with pytest.raises(ValueError):
            fields.select_type([])

    _sort = {
        'no extra, sorted': twice(
            tuple(FieldBase.for_tag(t)() for t in FieldTag
                  if t is not FieldTag.UNKNOWN),
            ),
        'no extra, scrambled': (
            shuffled(list(FieldBase.for_tag(t)() for t in FieldTag
                          if t is not FieldTag.UNKNOWN)),
            list(FieldBase.for_tag(t)() for t in FieldTag
                if t is not FieldTag.UNKNOWN),
            ),
        'with extra, sorted': twice(
            tuple(
                f for i, t in enumerate(FieldTag)
                for f in (FieldBase.for_tag(t)(str(i)), UnknownField(str(i)))
                ),
            ),
        'with extra, scrambled': (
            (NotesField('1'), UnknownField('2'), TensorNumsField('3'),
             UnknownField('4')),
            (TensorNumsField('3'), UnknownField('2'), NotesField('1'),
             UnknownField('4')),
            ),
        # NB: these examples are a bit unusual, as unknown fields
        # following notes are usually merged with the notes themselves
        'duplicates, sorted': twice((TensorNumsField('1'), UnknownField('2'),
                                     TensorNumsField('3'), UnknownField('4'),
                                     NotesField('5'), UnknownField('6'))),
        'duplicates, scrambled': (
            (NotesField('5'), TensorNumsField('1'), NotesField('7'),
             UnknownField('2'), TensorNumsField('3'), UnknownField('4'),
             UnknownField('6')),
            (TensorNumsField('1'), TensorNumsField('3'), NotesField('5'),
             UnknownField('2'), NotesField('7'), UnknownField('4'),
             UnknownField('6')),
            ),
        }

    @parametrize('values,expect', _sort.values(), ids=_sort)
    def test_sort(self, values, expect):
        """Check correct in-place sorting."""
        fields = FieldList(*values)
        fields.sort()
        assert fields == expect

    _remove = {
        'none': (),
        'first': (UnknownField('123'),),
        'two': (UnknownField('123'), UnknownField('123')),
        'three': (UnknownField('123'), UnknownField('123'),
                  UnknownField('123')),
        'two, different, sorted differently': (
            TensorNumsField('789'),
            UnknownField('123'),
            ),
        'all of them' : (
            UnknownField('456'),
            UnknownField('123'),
            UnknownField('123'),
            TensorNumsField('789'),
            UnknownField('123'),
            ),
        }

    @parametrize(to_remove=_remove.values(), ids=_remove)
    def test_remove_fields(self, to_remove):
        """Check outcome of removal of fields."""
        items = (
            UnknownField('123'),
            UnknownField('456'),
            UnknownField('123'),
            UnknownField('123'),
            TensorNumsField('789'),
            )
        fields = FieldList(*items)
        len_before = len(fields)
        fields.remove_fields(*to_remove)
        assert len(fields) == len_before - len(to_remove)

        # Check that the right items were removed
        counts_expect = {f: items.count(f) for f in items}
        for item in to_remove:
            counts_expect[item] -= 1
        counts_after = {f: fields.count(f) for f in items}
        assert counts_after == counts_expect

    _remove_invalid = {
        'all not there': (TensorNumsField(),),
        'second not there': (UnknownField(), TensorNumsField()),
        'first not there': (TensorNumsField(), UnknownField()),
        'not a field': (UnknownField(), '123',),
        'twice the same, but it is there only once': (UnknownField(),) * 2,
        }

    @parametrize(invalid=_remove_invalid.values(), ids=_remove_invalid)
    def test_remove_fields_raises(self, invalid):
        """Check complaints when removing invalid fields."""
        fields = FieldList(UnknownField())
        fields_before = fields.copy()
        with pytest.raises(ValueError):
            fields.remove_fields(*invalid)
        assert fields == fields_before

    _replace = {  # index to replace, new field
        # Notice that the indices here correspond to the FIRST
        # occurrence of the field that will be replaced.
        'unique before and after': (0, TensorNumsField('abc')),
        'duplicate before, not after': (0, TensorNumsField('abc')),
        'duplicate before and after': (0, TensorNumsField('789')),
        'more duplicates': (1, TensorNumsField('789')),
        }

    @parametrize('ind,new_field', _replace.values(), ids=_replace)
    def test_replace(self, ind, new_field):
        """Check correct replacement of a field."""
        fields = FieldList(
            UnknownField('123'),
            UnknownField('456'),
            UnknownField('123'),
            TensorNumsField('789'),
            )
        old_field = fields[ind]
        fields_before = fields.copy()
        fields.replace(old_field, new_field)
        assert any(after is not before
                   for before, after in zip(fields_before, fields))
        assert fields[ind] is new_field
        assert not any(f is new_field
                       for i, f in enumerate(fields)
                       if i != ind)

    def test_replace_not_there(self):
        """Check complaints when trying to replace a non-present field."""
        fields = FieldList()
        with pytest.raises(ValueError, match='not in'):
            fields.replace(UnknownField(), TensorNumsField())
