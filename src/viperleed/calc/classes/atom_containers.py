"""Module atom_collections of viperleed.tleemlib.classes.

Contains definitions of collections (i.e., containers) of Atoms.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-06'
__license__ = 'GPLv3+'

from abc import abstractmethod
from collections.abc import Container, Iterable, MutableSequence, Sequence

from viperleed.calc.classes import atom


class AtomContainerError(Exception):
    """Base exception for containers of Atom objects."""


class AtomListOutOfDateError(AtomContainerError):
    """AtomList needs to update_atoms_map."""

    def __init__(self, *args):
        """Initialize instance."""
        if not args:
            args = ('Call .update_atoms_map',)
        super().__init__(*args)


class DuplicateAtomsError(AtomContainerError):
    """There are multiple atoms with the same num."""

    def __init__(self, *args):
        """Initialize instance."""
        if not args:
            args = ('Multiple atoms with the same .num',)
        super().__init__(*args)


# pylint: disable=too-few-public-methods
# Disabled because AtomCollection is essentially only meant to
# be the base class for collections of Atom objects. For now,
# the only reasonable methods are n_atoms and __contains__.
class AtomContainer(Container):
    """Base class for all other collections of Atom objects."""

    @property
    @abstractmethod
    def n_atoms(self):
        """Return the number of atoms in this AtomList."""
# pylint: enable=too-few-public-methods


class AtomList(AtomContainer, MutableSequence):
    """A sequence of Atom objects.

    Attributes
    ----------
    strict : bool
        Whether this AtomList only allows unique .num values. This
        attribute can temporarily be turned off if atom numbers are
        made unique later. After .num(s) are unique, it is advisable
        to restore strict to True and to call update_atoms_map.
        A strict AtomList has some of its operations optimized
        for O(1) access.
    """

    def __init__(self, *atoms, strict=True):
        """Initialize instance."""
        if atoms and len(atoms) == 1 and isinstance(atoms[0], Iterable):
            atoms = atoms[0]
        self._atoms = list(atoms)
        if not all(isinstance(at, atom.Atom) for at in self):
            raise TypeError(
                f'{type(self).__name__}: only Atom objects allowed'
                )
        self._map = {at.num: at for at in self._atoms}
        self.strict = strict  # True == allow only unique num
        self._sort_map = {}   # Atom: index, for restoring sorting

        if self.strict and len(self._map) != self.n_atoms:
            raise DuplicateAtomsError

    def __contains__(self, value):
        """Return whether an atom is in this list."""
        if not isinstance(value, atom.Atom):
            return False
        try:
            atom_in_self = self.get(value.num)
        except KeyError:
            return False
        return value is atom_in_self

    def __delitem__(self, index):
        """Remove atom(s) at index (or slice)."""
        try:
            to_remove = self._atoms[index]
        except IndexError:
            raise IndexError('list assignment index out of range') from None
        del self._atoms[index]
        if not isinstance(index, slice):
            to_remove = [to_remove]
        for atom_ in to_remove:
            del self._map[atom_.num]

    def __eq__(self, other):
        """Return whether this list is equal to another one."""
        try:
            equal = self.n_atoms == len(other)
        except TypeError:
            equal = False
        if not equal:
            return NotImplemented
        equal = all(at_self is at_other
                    for at_self, at_other in zip(self, other))
        if not equal:
            return NotImplemented
        return True

    def __iter__(self):
        """Return an iterator of atoms in this AtomList."""
        # Reimplemented as the Sequence ABC uses a for and __getitem__
        return iter(self._atoms)

    def __getitem__(self, index):
        """Return the Atom(s) at index (or slice)."""
        return self._atoms[index]

    def __len__(self):
        """Return the number of atoms in this list."""
        return len(self._atoms)

    def __repr__(self):
        """Return a string representation of this atom list."""
        _repr_ = f'{type(self).__name__}('
        _repr_ += ', '.join(str(at) for at in self)
        return _repr_ + f', strict={self.strict})'

    def __reversed__(self):
        """Return a reverse iterator of atoms in this AtomList."""
        return reversed(self._atoms)

    def __setitem__(self, index, value):
        """Assign the atom(s) at index (or slice) to value."""
        # When a slice, we expect an iterable, otherwise a single
        # element. For the case of a slice, make sure we consume a
        # possible iterator right away, as we'll need it twice
        if isinstance(value, Sequence):
            new_atoms = value
        elif isinstance(value, Iterable):
            new_atoms = tuple(value)
        elif not isinstance(index, slice):
            new_atoms = (value,)

        if isinstance(index, slice):
            replaced = set(self._atoms[index])
        else:
            replaced = {self._atoms[index]}

        try:
            _invalid = {at for at in new_atoms
                        if at not in replaced and at.num in self._map}
        except AttributeError:
            raise TypeError(
                f'{type(self).__name__}: only Atom objects allowed'
                ) from None
        if _invalid and self.strict:
            raise DuplicateAtomsError(
                f'Atom(s) {_invalid} have the same .num as others '
                f'already present in this {type(self).__name__}'
                )
        self._map.update((at.num, at) for at in new_atoms)
        if not isinstance(index, slice):
            new_atoms = new_atoms[0]
        self._atoms[index] = new_atoms

    @property
    def n_atoms(self):
        """Return the number of atoms in this AtomList."""
        return len(self)

    def clear(self):
        """Remove all atoms."""
        self._atoms.clear()
        self._map = {}
        self._sort_map = {}

    def copy(self):
        """Return a copy of this AtomList."""
        return self.__class__(*self, strict=self.strict)

    def count(self, value):
        """Return the number of occurrences of value."""
        try:
            present = self._map[value.num] is value
        except (KeyError, AttributeError):
            present = False
        if self.strict:
            return 1 if present else 0
        if not present:
            return 0
        return self._atoms.count(value)

    _no_default = object()

    def get(self, num, default=_no_default):
        """Return the atom with a given number."""
        if not self.strict:
            raise NotImplementedError(
                f'Cannot use .get for non-strict {type(self).__name__}. You '
                f'can use, e.g., next(at for at in atoms if at.num=={num}) '
                '(for the first occurrence)'
                )
        try:
            return self._map[num]
        except KeyError:
            pass
        if default is not self._no_default:
            return default
        if self.n_atoms != len(self._map):
            raise AtomListOutOfDateError
        raise KeyError(num)

    def index(self, value, *args):
        """Return the index of an atom."""
        return self._atoms.index(value, *args)

    def insert(self, index, value):
        """Insert a single atom at a given index."""
        if not isinstance(value, atom.Atom):
            raise TypeError(
                f'{type(self).__name__}: only Atom objects allowed'
                )
        if value.num in self._map and self.strict:
            raise DuplicateAtomsError(
                f'{value} has the same .num as antoher one '
                f'already present in this {type(self).__name__}'
                )
        self._atoms.insert(index, value)
        self._map[value.num] = value

    def restore_sorting(self, sort_map=None):
        """Restore a sorting order (or the last one stored)."""
        if not sort_map:
            sort_map = self._sort_map
        if not sort_map:
            return
        try:
            self.sort(key=sort_map.get)
        except (TypeError, KeyError) as exc:
            raise RuntimeError(
                'Impossible to restore sorting as some elements '
                'have been added since save_sorting() was called'
                ) from exc

    def reverse(self):
        """Reverse in place."""
        self._atoms.reverse()

    def save_sorting(self):
        """Store the current sorting so it can be later restored."""
        self._sort_map = dict(zip(self, range(self.n_atoms)))
        return self._sort_map

    def sort(self, *, key=None, reverse=False):
        """Sort this list of atoms."""
        self._atoms.sort(key=key, reverse=reverse)

    def update_atoms_map(self):
        """Refresh the {num: atom} map.

        Typically called after atom numbers have been modified.

        Raises
        ------
        DuplicateAtomsError
            If multiple atoms in self have the same .num.
        """
        self._map = {at.num: at for at in self._atoms}
        if len(self) != len(self._map):
            raise DuplicateAtomsError
