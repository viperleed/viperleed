"""Tests for module viperleed.calc.classes.atom_containers."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-06'
__license__ = 'GPLv3+'

from operator import attrgetter

import pytest
from pytest_cases import fixture, parametrize, parametrize_with_cases

from viperleed.calc.classes.atom import Atom
from viperleed.calc.classes.atom_containers import AtomList
from viperleed.calc.classes.atom_containers import AtomListOutOfDateError
from viperleed.calc.classes.atom_containers import DuplicateAtomsError


class CasesAtomList:
    """Collection of pytest cases for testing AtomList objects."""

    def case_one_atom(self):
        """Return a tuple with single Atom."""
        return (Atom('H', (0, 0, 0), 1, None),)

    def case_three_atoms(self):
        """Return three distinct atoms."""
        return (*self.case_one_atom(),
                Atom('O', (0.1, 0.3, 0.4), 2, None),
                Atom('C', (0.5, 0.1, 0.8), 3, None))


THREE_ATOMS = CasesAtomList().case_three_atoms


@fixture(name='make_atomlist', scope='session')
def factory_make_atomlist():
    """Return an AtomList from a sequence of atoms."""
    def _make(*atoms, **kwargs):
        return AtomList(*atoms, **kwargs)
    return _make


class TestAtomList:
    """Collection of tests for AtomList class."""

    @parametrize_with_cases('atoms', cases=CasesAtomList)
    def test_creation(self, make_atomlist, atoms):
        """Check that creation of an AtomList is successful."""
        atom_list = make_atomlist(atoms)
        assert atom_list.n_atoms == len(atoms)
        assert atom_list == atoms

    def test_clear(self, make_atomlist):
        """Check successful emptying of an AtomList."""
        atom_list = make_atomlist(*THREE_ATOMS())
        atom_list.clear()
        assert not atom_list.n_atoms

    _strict_counts = {True: (1, 1, 0),
                      False: (2, 1, 0)}

    @parametrize('strict,expected_counts', _strict_counts.items(),
                 ids=(f'strict={k}' for k in _strict_counts))
    def test_count(self, strict, expected_counts, make_atomlist):
        """Ensure atom counts are as expected."""
        first, second, _ = all_atoms = THREE_ATOMS()
        atoms = first, second
        if not strict:
            atoms += (first,)
        atom_list = make_atomlist(*atoms, strict=strict)
        counts = tuple(atom_list.count(at) for at in all_atoms)
        assert counts == expected_counts

    def test_does_not_contain_non_atom(self, make_atomlist):
        """Check that a non-Atom is never in an AtomList."""
        atom_list = make_atomlist(THREE_ATOMS())
        assert 10 not in atom_list

    @parametrize(other=(10, [1, 2, 5]))
    def test_not_equal(self, other, make_atomlist):
        """Check that a non-Atom never compares equal to an AtomList."""
        atom_list = make_atomlist(THREE_ATOMS())
        assert atom_list != other

    def test_delitem(self, make_atomlist):
        """Check the successful deletion of an Atom."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(*atoms)
        del atom_list[1]
        assert len(atom_list) == 2
        assert atoms[1] not in atom_list

    def test_get(self, make_atomlist):
        """Check correct return and exceptions in AtomList.get()."""
        first, second, _ = THREE_ATOMS()
        atom_list = make_atomlist(second, first)
        assert atom_list.get(1) is first
        assert atom_list.get(2) is second
        with pytest.raises(KeyError):
            assert atom_list.get(3)

    def test_insert(self, make_atomlist):
        """Check the successful insertion of an Atom."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(*atoms[:2])
        atom_list.insert(1, atoms[2])
        assert atom_list.n_atoms == 3
        assert atom_list[1] == atoms[-1]

    def test_restore_sort(self, make_atomlist):
        """Check that sorting is restored successfully."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(atoms)
        atom_list.save_sorting()
        atom_list.sort(key=attrgetter('el'))
        atom_list.restore_sorting()
        assert atom_list == atoms

    def test_setitem(self, make_atomlist):
        """Check successful assignment of an AtomList item."""
        first, second, third = THREE_ATOMS()
        atom_list = make_atomlist(first, third)
        atom_list[1] = second
        assert atom_list[1] is second

    def test_setitem_slice(self, make_atomlist):
        """Check successful assignment of an AtomList item."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(atoms[0])
        atom_list[1:1] = atoms[1:]
        assert atom_list == atoms

    def test_setitem_genexpr(self, make_atomlist):
        """Check successful assignment of an AtomList item."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(atoms[0])
        atom_list[1:1] = (at for at in atoms[1:])
        assert atom_list == atoms

    def test_update_atoms_map(self, make_atomlist):
        """Ensure {num: Atom} map is created successfully."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(*atoms)
        expected_nums = [at.num for at in atom_list]

        # pylint: disable=protected-access
        # There's no public interface to directly access the _map
        assert all(n in atom_list._map for n in expected_nums)
        for atom in atoms:
            atom.num += 1
        atom_list.update_atoms_map()
        assert all(n+1 in atom_list._map for n in expected_nums)
        assert len(atom_list._map) == 3

    def test_sort(self, make_atomlist):
        """Test successful sorting of an AtomList."""
        first, second, third = THREE_ATOMS()
        atom_list = make_atomlist(second, first, third)
        atom_list.sort(key=attrgetter('el'))
        assert atom_list == [third, first, second]

    def test_reverse(self, make_atomlist):
        """Test successful reversing of an AtomList."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(*atoms)
        atom_list.reverse()
        assert atom_list == atoms[::-1]


class TestAtomListRaises:
    """Collection of AtomList tests for paths that that raise exceptions."""

    def test_delitem_indexerror(self, make_atomlist):
        """Check complaints about atom duplication during initialization."""
        atom_list = make_atomlist()
        with pytest.raises(IndexError):
            del atom_list[0]

    @parametrize_with_cases('atoms', cases=CasesAtomList)
    def test_duplicate_atoms_init(self, make_atomlist, atoms):
        """Check complaints about atom duplication during initialization."""
        with pytest.raises(DuplicateAtomsError):
            make_atomlist(*atoms, *atoms)

    def test_duplicate_atoms_insert(self, make_atomlist):
        """Check complaints about atom duplication during initialization."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(atoms)
        with pytest.raises(DuplicateAtomsError):
            atom_list.insert(0, atoms[0])

    def test_duplicate_atoms_setitem(self, make_atomlist):
        """Check complaints about atom duplication during initialization."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(atoms)
        with pytest.raises(DuplicateAtomsError):
            atom_list[2] = atoms[0]

    @parametrize_with_cases('atoms', cases=CasesAtomList)
    def test_duplicate_atoms_update(self, make_atomlist, atoms):
        """Check complaints about atom duplication during map updating."""
        atom_list = make_atomlist(*atoms, *atoms, strict=False)
        with pytest.raises(DuplicateAtomsError):
            atom_list.update_atoms_map()

    def test_get_non_strict_raises(self, make_atomlist):
        """Check that calling get() on a non-strict AtomList raises."""
        atom_list = make_atomlist(*THREE_ATOMS(), strict=False)
        with pytest.raises(NotImplementedError):
            atom_list.get(10)

    def test_get_out_of_date(self, make_atomlist):
        """Check that calling .get() on an outdated AtomList raises."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(*atoms, *atoms, strict=False)
        atom_list.strict = True
        with pytest.raises(AtomListOutOfDateError):
            assert atom_list.get(10)

    def test_restore_sorting_added_atom(self, make_atomlist):
        """Check that we cannot restore sorting if atoms are added."""
        atoms = THREE_ATOMS()
        atom_list = make_atomlist(atoms[0])
        atom_list.save_sorting()
        atom_list.extend(atoms[1:])
        with pytest.raises(RuntimeError):
            atom_list.restore_sorting()

    def test_typerror_init(self, make_atomlist):
        """Check that initializing with non-Atoms raises."""
        with pytest.raises(TypeError):
            make_atomlist(None)

    def test_typerror_insert(self, make_atomlist):
        """Check that initializing with non-Atoms raises."""
        atom_list = make_atomlist(THREE_ATOMS())
        with pytest.raises(TypeError):
            atom_list.insert(0, None)

    def test_typerror_setitem(self, make_atomlist):
        """Check that initializing with non-Atoms raises."""
        atom_list = make_atomlist(THREE_ATOMS())
        with pytest.raises(TypeError):
            atom_list[0] = None
