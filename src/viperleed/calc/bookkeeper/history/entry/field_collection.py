"""Module field_collection of viperleed.calc.bookkeeper.history.entry.

Defines base classes for collections of FieldBase objects.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-07-31'
__license__ = 'GPLv3+'

from collections import defaultdict
from collections.abc import MutableSequence
from copy import deepcopy
from typing import Generator

from viperleed.calc.bookkeeper.history.entry.enums import DuplicateType
from viperleed.calc.bookkeeper.history.entry.enums import FieldTag
from viperleed.calc.bookkeeper.history.entry.field import FieldBase
from viperleed.calc.bookkeeper.history.entry.field import UnknownField
from viperleed.calc.bookkeeper.history.entry.notes_field import NotesField
from viperleed.calc.bookkeeper.history.errors import FieldsScrambledError
from viperleed.calc.lib.sequence_utils import conditional_sort


class FieldList(MutableSequence):
    """A sorted sequence of fields with fast access to its members.

    Fast access is guaranteed as long as all the items are hashable
    and there are no duplicates.
    """

    def __init__(self, *fields):
        """Initialize instance."""
        self._seq = list(fields)
        if len(self._seq) == 1 and isinstance(self._seq[0], Generator):
            self._seq = list(self._seq[0])
        self._check_item_types(*self._seq)
        self._maps = {
            'by_tag': {},
            'by_type': defaultdict(list),
            'by_index': {},
            # value_to_index is the inverse of by_index, but is filled
            # only if there's no duplicates and if all are hashable
            'value_to_index': {},
            }
        self._make_maps()

    def __bool__(self):
        """Return whether there is anything in this FieldList."""
        return bool(self._seq)

    def __contains__(self, value):
        """Return whether `value` is in this FieldList."""
        _map = self._maps['value_to_index']
        if not _map:  # Some item is non-hashable or duplicate
            return super().__contains__(value)
        try:
            return value in _map
        except TypeError:  # value non-hashable
            return super().__contains__(value)

    def __delitem__(self, index):
        """Remove the item(s) at `index`."""
        del self._seq[index]
        self._make_maps()

    def __eq__(self, other):
        """Return whether self == other."""
        try:
            same_length = len(self) == len(other)
        except TypeError:  # No __len__
            return NotImplemented
        if not same_length:
            return NotImplemented
        if any(self_i != other_i for self_i, other_i in zip(self, other)):
            return NotImplemented
        return True

    def __getitem__(self, index):
        """Return the items(s) at `index`."""
        try:
            return self._maps['by_index'][index]
        except (KeyError, TypeError, AttributeError):
            # Not there, not hashable, no _maps yet
            return self._seq[index]

    def __len__(self):
        """Return the number of items in this FieldList."""
        return len(self._seq)

    def __repr__(self):
        """Return a literal string version of this FieldList."""
        sep = '\n    '
        items = sep.join(f'{repr(item)},' for item in self)
        items = f'{sep}{items}{sep}' if items else items
        return f'{type(self).__name__}({items})'

    def __setitem__(self, index, value):
        """Assign one or more `index` as `value`."""
        if isinstance(index, slice):
            value = tuple(value)  # In case it's an iterator
            self._check_item_types(*value)
        else:
            self._check_item_types(value)
        self._seq[index] = value
        self._make_maps()

    @property
    def has_only_pure_comments(self):
        """Return whether all fields in this sequence are comment-only."""
        types = tuple(t for t, v in self._maps['by_type'].items() if v)
        return types == (UnknownField,)

    def check_sorted(self):
        """Raise FieldsScrambledError if items are not sorted as expected."""
        cls = type(self)
        # Use only known fields, and skip duplicates
        first_items = [duplicates[0]
                       for t, duplicates in self._maps['by_type'].items()
                       if t is not UnknownField and duplicates]
        sorted_items = cls(*first_items)
        sorted_items.sort()
        if first_items != sorted_items:
            raise FieldsScrambledError(
                'Lines are not in the expected order. Should be\n'
                + '\n'.join(str(f) for f in sorted_items)
                )

    def copy(self):
        """Return a shallow copy of this FieldList."""
        cls = type(self)
        return cls(*self._seq)

    def duplicates(self, field_type):
        """Return duplicate fields of a given type.

        Parameters
        ----------
        field_type : type(FieldBase) or FieldBase
            The field type to be checked.

        Returns
        -------
        duplicate_reason : DuplicateType
            Which kind of duplication there is.
        duplicates : tuple
            The duplicate fields.
        """
        if isinstance(field_type, FieldBase):
            field_type = type(field_type)
        by_type = self.select_type(field_type)
        # pylint: disable-next=magic-value-comparison  # Clear enough
        if len(by_type) < 2:
            return DuplicateType.NONE, tuple()
        first = by_type[0]
        if all(f == first for f in by_type[1:]):
            return DuplicateType.IDENTICAL, by_type
        duplicate_reason = (DuplicateType.NOTES if field_type is NotesField
                            else DuplicateType.DIFFERENT)
        return duplicate_reason, by_type

    def first_by_type(self, field_type):
        """Return the first item of type `field_type` or None."""
        by_type = self.select_type(field_type)
        return by_type[0] if by_type else None

    def index(self, value, start=0, stop=None):
        """Return the first index of `value`."""
        _map = self._maps['value_to_index']
        if start or stop or not _map:
            args = (start, stop) if stop is not None else (start,)
            return self._seq.index(value, *args)
        try:
            return _map[value]
        except (KeyError, TypeError):
            raise ValueError(f'{value} not in list') from None

    def insert(self, index, value):
        """Insert `value` at `index`."""
        self._check_item_types(value)
        self._seq.insert(index, value)
        self._make_maps()

    def insert_sorted(self, value):
        """Insert `value` in the right sort order."""
        try:
            index = self._find_sorted_insert_index(value)
        except IndexError:  # Should go last
            self.append(value)
        else:
            self.insert(index, value)

    def remove_fields(self, *fields):
        """Remove one or more fields from this FieldList."""
        # Store some backup values in case any fields are invalid
        backup_seq = self._seq.copy()
        backup_maps = deepcopy(self._maps)
        for field in fields:
            # We do the "expensive" item-by-item removal in order to
            # ensure that the value-to-index maps are always up to
            # date. This prevents us from having to (i) fiddle with
            # keeping track of possible repeated indices in case of
            # duplicates, (ii) doing the previous in case the state
            # of duplication changes as a result of removal of some
            # items, and (iii) making sure we remove the right fields
            # (the map is out of date as soon as one item is removed).
            try:
                self.remove(field)
            except ValueError:
                self._seq = backup_seq
                self._maps = backup_maps
                raise

    def replace(self, old_field, new_field):
        """Replace (the first occurrence of) an item with another one."""
        index = self.index(old_field)
        self[index] = new_field

    def select_type(self, field_type):
        """Return a tuple of fields with a given type."""
        if isinstance(field_type, FieldBase):
            field_type = type(field_type)
        map_ = self._maps['by_type']
        try:
            return tuple(map_[field_type])
        except TypeError as exc:  # No KeyError with defaultdict
            raise ValueError(f'Unkown field type {field_type}.') from exc

    def sort(self):
        """Sort this FieldList in place, keeping UnknownField unmoved."""
        _map = {t: ind for ind, t in enumerate(FieldTag)}
        self._seq = conditional_sort(
            self._seq,
            skip=lambda field: isinstance(field, UnknownField),
            key=lambda field: _map[field.tag],
            )
        self._make_maps()

    def _check_item_types(self, *fields):
        """Raise TypeError unless all items are FieldBase."""
        if not all(isinstance(f, FieldBase) for f in fields):
            raise TypeError(
                'Can only contain FieldBase objects. Found '
                + ', '.join(f'{type(i).__name__}' for i in fields)
                )

    def _find_exact_index(self, value):
        """Find the index of the item that is the same object as `value`."""
        try:
            return next(i for i, f in enumerate(self) if f is value)
        except StopIteration:
            raise ValueError(f'{value} not in {self}') from None

    def _find_sorted_insert_index(self, value):
        """Return an index for inserting `value` in a sorted manner."""
        by_tag = self._maps['by_tag']
        try:
            with_same_tag = by_tag[value.tag]
        except (KeyError, AttributeError) as exc:  # Not a field
            raise TypeError(f'{repr(value)} is not a field') from exc

        if with_same_tag:  # Goes last among those with the same tag
            previous_item = with_same_tag[-1]
            return self._find_exact_index(previous_item) + 1
        # None yet with the same tag. Find an appropriate insertion
        # index by skimming through the rest of the items, but always
        # insert after UnknownField rather than before.
        index = self._find_sorted_insert_index_new_type(value)
        while self[index].tag is FieldTag.UNKNOWN:
            index += 1
        return index

    def _find_sorted_insert_index_new_type(self, value):
        """Return an index for inserting a `value` with a new type in order."""
        by_tag = self._maps['by_tag']
        # Find out the first item in self with the next tag. Find first
        # the next 'bucket', that is, the first item in the by_tag map
        # coming after this value that has at least one item. Notice
        # that we do not iterate over by_tag.items(), as it has also
        # non-tagged fields, which we want to always leave at the end.
        # The correct bucket is the first one coming after the tag of
        # this one, but may not be the one right after it, as we may
        # be missing entirely items of the "next tag" too. Also, make
        # sure to skip FieldTag.UNKNOWN fields in this search, as they
        # may be anywhere and "do not matter".
        if value.tag is FieldTag.UNKNOWN:
            raise IndexError    # Goes last
        tags = [t for t in FieldTag if t is not FieldTag.UNKNOWN]
        start = tags.index(value.tag)
        potential_indices = []  # items: (bucket_tag, insertion_index)
        try:
            next_tag = next(t for t in tags[start:] if by_tag[t])
        except StopIteration:   # May need to go last
            potential_indices.append((None, len(self)))
        else:
            next_items = by_tag[next_tag]
            potential_indices.append((
                next_tag,
                self._find_exact_index(next_items[0]),
                ))

        # Now do the same with the previous tag
        try:
            previous_tag = next(t for t in reversed(tags[:start]) if by_tag[t])
        except StopIteration:   # May need to go first
            potential_indices.append((None, 0))
        else:
            previous_items = by_tag[previous_tag]
            potential_indices.append((
                previous_tag,
                self._find_exact_index(previous_items[-1]) + 1,
                ))
        if all(t for t, _ in potential_indices):
            # Take the one that is closest tag-wise
            _, index = min((abs(start - tags.index(t)), ind)
                           for t, ind in potential_indices)
            return index
        try:
            return next(ind for t, ind in potential_indices if t)
        except StopIteration:  # Last resort. Place it last
            raise IndexError from None

    def _make_maps(self):
        """Populate self._maps."""
        self._make_non_indexed_maps()
        self._make_indexed_maps()

    def _make_non_indexed_maps(self):
        """Make all maps that are not index based."""
        self._make_map_by_type()
        self._make_map_by_tag()

    def _make_indexed_maps(self):
        """Map indices to elements and viceversa for fast access."""
        self._maps['by_index'] = map_ = dict(enumerate(self._seq))
        try:
            to_index = {v: k for k, v in map_.items()}
        except TypeError:  # Some non-hashable
            to_index = {}
        if to_index and len(to_index) != len(map_):  # Duplicates
            to_index = {}
        self._maps['value_to_index'] = to_index

    def _make_map_by_tag(self):
        """Populate the internal mapping that stores fields by their tag."""
        # This is somewhat similar to the by_type one, but helps with
        # some searching operations. Rather than using a defaultdict,
        # make all of the tags, so they come already sorted correctly.
        map_ = {tag: [] for tag in FieldTag}
        map_[None] = []  # Non-tagged FieldBase always at the end
        # Notice that we DO NOT use "for field in self"
        # as that one may access an outdated map!
        for field in self._seq:
            map_[field.tag].append(field)
        self._maps['by_tag'] = map_

    def _make_map_by_type(self):
        """Populate the internal mapping that stores fields by their type."""
        map_ = defaultdict(list)
        # Notice that we DO NOT use "for field in self"
        # as that one may access an outdated map!
        for field in self._seq:
            map_[type(field)].append(field)
        self._maps['by_type'] = map_
