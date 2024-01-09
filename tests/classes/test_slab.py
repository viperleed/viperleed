"""Tests for viperleed.tleedmlib.classes.slab.

Created on 2023-07-28

@author: Alexander M. Imre (@amimre)
@author: Michele Riva (@michele-riva)
"""

from copy import deepcopy
import operator
from random import shuffle
from pathlib import Path
import sys

import numpy as np
import pytest
from pytest_cases import fixture, fixture_ref
from pytest_cases import parametrize, parametrize_with_cases

VPR_PATH = str(Path(__file__).resolve().parents[3])
if VPR_PATH not in sys.path:
    sys.path.append(VPR_PATH)

# pylint: disable=wrong-import-position
# Cannot do anything about it until we make viperleed installable
from viperleed.tleedmlib.base import pairwise
from viperleed.tleedmlib.base import NonIntegerMatrixError, SingularMatrixError
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.rparams import Rparams
from viperleed.tleedmlib.classes.slab import Slab
from viperleed.tleedmlib.classes.slab import slab_errors as err
from viperleed.tleedmlib.classes.slab import surface_slab
from viperleed.tleedmlib.classes.sym_entity import SymPlane

from .. import poscar_slabs
from ..helpers import not_raises
# pylint: enable=wrong-import-position

CasePOSCARSlabs = poscar_slabs.CasePOSCARSlabs
todo = pytest.mark.skip('To be implemented')


@fixture(name='shuffle_slab', scope='session')
def make_shuffled_slab():
    """Shuffle the atoms of a slab at random."""
    def _shuffle(slab):
        slab.atlist.strict = False  # Silence duplicate errors
        shuffle(slab.atlist)
        slab.atlist.strict = True   # Now all atoms should be unique
        slab.atlist.update_atoms_map()
    return _shuffle


class TestAtomTransforms:
    """Test simple transformations of the atoms of a slab."""

    def test_mirror(self, manual_slab_3_atoms):
        """Test the expected outcome of mirroring atoms of a simple slab."""
        slab = manual_slab_3_atoms
        mirrored_slab = deepcopy(slab)
        symplane = SymPlane((0, 0), (0, 1), abt=slab.ab_cell.T)
        mirrored_slab.mirror_atoms(symplane)
        assert all(at.is_same_xy(mir_at)
                   for at, mir_at in zip(slab, reversed(mirrored_slab)))

    def test_mirror_on_slanted_cell(self, make_poscar):
        """Test the expected outcome of mirroring atoms of a slanted slab."""
        slab, *_ = make_poscar(poscar_slabs.SLAB_36C_cm)
        slab.create_sublayers(0.1)
        sym_plane = SymPlane((0, 0), (0, 1), abt=slab.ab_cell.T)
        assert slab.is_mirror_symmetric(sym_plane, eps=1e-6)

    def test_rotation_symmetric_raises(self, ag100):
        """Check that no symmetry checks can be performed without sublayers."""
        slab, *_ = ag100
        slab.sublayers.clear()
        with pytest.raises(err.MissingSublayersError):
            slab.is_rotation_symmetric((0, 0), 4, 1e-3)

    def test_180_rotation(self, manual_slab_3_atoms):
        """Test the expected outcome of rotating atoms of a simple slab."""
        slab = manual_slab_3_atoms
        rotated_slab = deepcopy(slab)
        rotated_slab.rotate_atoms(order=2)
        assert all(at.is_same_xy(rot_at)
                   for at, rot_at in zip(slab, reversed(rotated_slab)))


class TestAtomsAndElements:
    """Collection of tests for atom additions/removals."""

    def test_empty_slab(self):
        """Check that an empty slab has no atoms, layers, etc..."""
        slab = Slab()
        assert not slab.atlist
        assert not slab.elements
        assert not slab.layers
        assert slab.planegroup == 'unknown'

    def test_add_one_atom_n_elements(self):
        """Check that adding one atom to a slab updates elements correctly."""
        slab = Slab()
        new_atom = Atom('C', (0, 0, 0), 1, slab)
        slab.atlist.append(new_atom)
        slab.update_element_count()
        assert new_atom.el in slab.elements
        assert slab.n_per_elem[new_atom.el] == 1

    @todo
    def test_atlist_is_not_list(self):                                          # TODO: Should consider various situations to make sure that no Slab method messes with the atlist
        """TODO"""

    def test_chemelem_upon_element_mix_changed(self, ag100):
        """Check that chemelem are changed when ELEMENT_MIX is."""
        slab, rpars, *_ = ag100
        rpars.ELEMENT_MIX = {'Ag': {'Fe', 'Co'}}
        slab.full_update(rpars)
        assert slab.chemelem == {'Fe', 'Co'}

    def test_remove_one_atom_n_elements(self, ag100):
        """Check that removing one atom updates elements correctly."""
        slab, *_ = ag100
        n_ag_atoms = slab.n_per_elem['Ag']
        slab.atlist.pop()
        slab.update_element_count()
        assert slab.n_per_elem['Ag'] == n_ag_atoms - 1

    remove_ats = {  # fe_atom_num_to_remove, is_bulk_atom
        'non bulk': (19, False),
        'bulk': (33, True),
        }

    @parametrize(make_bulk=(True, False))
    @parametrize('to_remove,in_bulk', remove_ats.values(), ids=remove_ats)
    def test_update_atom_numbers(self, make_bulk, in_bulk, to_remove, subtests):
        """Check correct removal of an Fe atom, with and without bulk."""
        fe3o4, rpars, *_ = CasePOSCARSlabs().case_poscar_fe3o4_001_cod()
        elems = 'Fe', 'O'
        n_fe, n_oxygen = (fe3o4.n_per_elem[k] for k in elems)
        if make_bulk:
            bulk = fe3o4.make_bulk_slab(rpars)
            n_fe_bulk, n_oxygen_bulk = (bulk.n_per_elem[k] for k in elems)
            old_at_nrs_bulk = [at.num for at in bulk]
        fe3o4.atlist.remove(fe3o4.atlist.get(to_remove))
        fe3o4.update_element_count()
        if make_bulk and in_bulk:
            bulk.atlist.remove(bulk.atlist.get(to_remove))
            bulk.update_element_count()

        fe3o4.update_atom_numbers()
        with subtests.test('nr. atoms'):
            assert fe3o4.n_per_elem['Fe'] == n_fe - 1
            assert fe3o4.n_per_elem['O'] == n_oxygen
        new_at_nrs = [at.num for at in fe3o4]
        with subtests.test('atom.num'):
            assert new_at_nrs == list(range(1, n_fe + n_oxygen))
        if not make_bulk:
            return
        with subtests.test('nr. atoms, bulk'):
            assert bulk.n_per_elem['O'] == n_oxygen_bulk
            assert bulk.n_per_elem['Fe'] == n_fe_bulk - (1 if in_bulk else 0)
        with subtests.test('atom.num, bulk'):
            assert all(at.num in new_at_nrs for at in bulk)
        with subtests.test('atom.num, bulk, changed'):
            assert [at.num for at in bulk] != old_at_nrs_bulk

@todo
class TestBulk3DOperations:
    """Tests for 3D symmetry operations of BulkSlab objects."""

    def test_get_candidate_layer_periods(self):
        """TODO"""

    def test_bulk_screw_symmetric(self):                                        # TODO: especially consider the cases where an should be on an axis but it is away, and those where it should be at an n-fold position but it isn't.
        """TODO"""

    def test_bulk_glide_symmetric(self):                                        # TODO: especially consider the cases where an should be on an plane but it is away, and those where it should be at a glide-symmetric position but it isn't.
        """TODO"""


class TestBulkDetectAndExtraBulk:
    """Collection of tests for adding bulk units to slabs."""

    @todo
    def test_detect_bulk(self):                                                 # TODO: also check that rp and sl are unchanged if it fails
        """TODO"""

    def test_detect_bulk_raises(self, ag100):
        """Check complaints for invalid inputs."""
        slab, rpars, *_ = ag100
        rpars.BULK_LIKE_BELOW = -0.5
        with pytest.raises(ValueError):
            slab.detect_bulk(rpars)

    @todo
    def test_with_extra_bulk_units(self):                                       # TODO: also check the number of bulk layers
        """TODO"""

    @todo
    def test_with_double_thickness_once(self):
        """TODO"""

    @todo
    def test_with_double_thickness_twice(self):
        """TODO"""


class TestBulkRepeat:
    """Collection of test for bulk-repeat finding and returning."""

    @todo
    def test_identify(self):
        """TODO"""

    def test_identify_raises_without_bulkslab(self, ag100):                     # TODO: also other raises
        """Check complaints when called without a bulk slab."""
        slab, *_ = ag100
        slab.bulkslab = None
        with pytest.raises(err.MissingBulkSlabError):
            slab.identify_bulk_repeat(.1)

    @todo
    def test_get(self):
        """TODO"""


@todo
class TestBulkUcell:
    """Tests concerning reduction of bulk unit cell and C vector."""

    def test_get_min_c(self):                                                   # TODO: also various raises
        """TODO"""

    def test_ensure_min_c(self):
        """TODO"""

    def test_apply_bulk_ucell_reduction(self):                                  # TODO: separately ab only and c_vec. Both with recenter==True/False
        """TODO"""

    def test_minimal_bulk_ab(self):                                             # TODO: SurfaceSlab
        """TODO"""


class TestCoordinates:
    """Collection of tests for Cartesian/fractional coordinates."""

    cart_from_frac = {  # update_origin, expected cartpos of 1st atom
        'no update': (False, (0.3, 0.8, -1.5)),
        'update': (True, (0.3, 0.8, 0)),
        }

    @parametrize('update_,expect', cart_from_frac.values(), ids=cart_from_frac)
    def test_cartesian_from_fractional(self, update_, expect,
                                       manual_slab_3_atoms):
        """Check correct update of Cartesian atom coordinates."""
        slab = manual_slab_3_atoms
        atom = slab.atlist[0]
        atom.pos = np.array([0.1, 0.2, 0.3])
        slab.update_cartesian_from_fractional(update_origin=update_)
        assert atom.cartpos == pytest.approx(expect)

    collapse_cart = {  # new_cartpos, expected cartpos after collapsing
        'large values': ((5, 6, 3), (2, 2, -2)),
        'small eps': ([1 - 1e-9, 2 + 1e-9, -3 - 1e-15], (1, 2, -3)),
        }

    @parametrize('new_cart,expect', collapse_cart.values(), ids=collapse_cart)
    def test_collapse_cartesian(self, new_cart, expect, manual_slab_3_atoms):
        """Check that Cartesian coordinates are correctly collapsed."""
        slab = manual_slab_3_atoms
        atom = slab.atlist[0]
        atom.cartpos = np.array(new_cart)
        slab.collapse_cartesian_coordinates()
        assert atom.cartpos == pytest.approx(expect, abs=1e-8)

    def test_collapse_fractional(self, manual_slab_3_atoms):                    # TODO: use both methods. Find especially cases that are 'problematic' with %1.0: e.g., 1-1e-8, 1-1e-9, 1-1e-15
        """Check that fractional coordinates are correctly collapsed."""
        slab = manual_slab_3_atoms
        atom = slab.atlist[0]
        atom.pos = np.array([1.1, 2.2, -3.3])
        slab.collapse_fractional_coordinates()
        assert atom.pos == pytest.approx([0.1, 0.2, 0.7])

    def test_fractional_from_cartesian(self, manual_slab_3_atoms):
        """Check correct update of fractional atom coordinates."""
        slab = manual_slab_3_atoms
        atom = slab.atlist[0]
        atom.cartpos = np.array([0.3, 0.8, -1.5])
        slab.update_fractional_from_cartesian()
        assert atom.pos == pytest.approx([0.1, 0.2, 0.3])


@todo
class TestDuplicateAtoms:
    """Tests for checking detection and removal of duplicate atoms."""

    @parametrize(info=poscar_slabs.WITH_DUPLICATE_ATOMS)
    def test_remove_duplicates(self, info, make_poscar):                        # TODO: n_atoms, raises, others? check method
        """Check correct removal of duplicate atoms."""
        slab, rpars, *_ = make_poscar(info)
        n_atoms_before = slab.n_atoms
        with not_raises(err.SlabError):
            slab.remove_duplicate_atoms(rpars.SYMMETRY_EPS,
                                        rpars.SYMMETRY_EPS.z)
        assert n_atoms_before > slab.n_atoms
        with not_raises(err.AtomsTooCloseError):
            slab.check_atom_collisions()

    @parametrize(info=poscar_slabs.WITH_DUPLICATE_ATOMS)
    def test_with_duplicate_atoms(self, info, make_poscar):
        """Check that POSCARs with duplicates are handled correctly."""
        slab, *_ = make_poscar(info)
        with pytest.raises(err.AtomsTooCloseError):
            slab.check_atom_collisions()

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_without_duplicates(self, args):
        """Check that POSCARs without duplicates are handled correctly."""
        slab, *_ = args
        with not_raises(err.AtomsTooCloseError):
            slab.check_atom_collisions()


class TestEquivalence:
    """Collection of tests for the is_equivalent method."""

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_equivalent_copy(self, args, shuffle_slab):
        """Check that a slab is equivalent to its deepcopy."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        shuffle_slab(slab_copy)  # Make sure order does not matter
        assert slab.is_equivalent(slab_copy, eps=1e-3)

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_equivalent_translated_2d(self, args, shuffle_slab):
        """Check equivalence by translating in plane by a unit vector."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        shuffle_slab(slab_copy)  # Make sure order does not matter
        for atom in slab_copy:
            atom.pos[:2] += np.random.randint(-10, 10, size=2)
        slab_copy.update_cartesian_from_fractional()
        assert slab.is_equivalent(slab_copy, eps=1e-3)

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_slab_translated_not_equivalent(self, args):
        """Check that a translated slab is not equivalent."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        for atom in slab_copy:
            atom.pos[0] -= 0.1
        slab_copy.update_cartesian_from_fractional()
        assert not slab.is_equivalent(slab_copy, eps=1e-10)

    def test_non_slab_not_equivalent(self, ag100):
        """Check equivalence with a non-slab object."""
        slab, *_ = ag100
        assert not slab.is_equivalent('not a slab')

    def test_more_layers_not_equivalent(self, ag100):
        """Check (in)equivalence of two slabs with different nr. of layers."""
        slab, rpars, *_ = ag100
        thicker, _ = slab.with_extra_bulk_units(rpars, 1)
        assert not slab.is_equivalent(thicker)

    def test_more_atoms_not_equivalent(self, ag100):
        """Check (in)equivalence of two slabs with different nr. of atoms."""
        slab, *_ = ag100
        two_by_one = np.diag((2, 1))
        assert not slab.is_equivalent(slab.make_supercell(two_by_one))


# pylint: disable-next=protected-access
@pytest.mark.skipif(not surface_slab._HAS_ASE, reason='No ASE module')
class TestFromAse:
    """Collection of tests for the from_ase class method."""

    def test_no_ase(self):
        """Check complaints when no ase module is present."""
        # pylint: disable=protected-access
        surface_slab._HAS_ASE = False
        with pytest.raises(ModuleNotFoundError):
            Slab.from_ase(None)
        surface_slab._HAS_ASE = True

    def test_not_an_ase_atoms(self):
        """Check complaints for the wrong type."""
        with pytest.raises(TypeError):
            Slab.from_ase('invalid')


class TestMakeBulkSlab:
    """Collection of tests for the make_bulk_slab method."""

    _invalid = {
        'no layers': (Slab(), err.MissingLayersError),
        'one layer': (fixture_ref('manual_slab_3_atoms'),
                      err.TooFewLayersError),
        }

    @todo
    def test_valid(self):                                                       # TODO: also LOG of warning if a_bulk > b_bulk
        """TODO"""

    @parametrize('slab,exc', _invalid.values(), ids=_invalid)
    def test_invalid(self, slab, exc):
        """Check complaints for invalid conditions when making bulk."""
        with pytest.raises(exc):
            slab.make_bulk_slab(Rparams())


class TestProperties:
    """Collection of tests for various @property of Slab."""

    def test_thickness(self, ag100):
        """Check expected slab thickness."""
        slab, *_ = ag100
        assert slab.thickness == pytest.approx(10.18233, abs=1e-4)

    def test_vacuum_gap(self, ag100):
        """Check expected thickness of the vacuum gap."""
        slab, *_ = ag100
        assert slab.vacuum_gap == pytest.approx(10.18233, abs=1e-4)


# This class currently contains NO TESTS, as the previous tests were
# somewhat faulty. The purpose of restore_ori_state is to:
#    (i) convert the current positions, vibrations and occupations
#        into VIBROCC offsets, and
#   (ii) fully clear the displacements of all atoms
# The tests were originally set up by @amimre, who found a bug that
# concerns runs where multiple refcalc-search pairs exists in RUN.
# This bug is documented in Issue #107.
class TestRestoreOristate:
    """Collection of tests for reverting a slab to its ref-calc state."""

    @pytest.fixture(name='slab_and_copy')
    def fixture_slab_and_copy(self, ag100_with_displacements_and_offsets):
        """Return a Ag(100) slab and its deepcopy."""
        slab, *_ = ag100_with_displacements_and_offsets
        return slab, deepcopy(slab)

    @staticmethod
    def check_displacements_equal(atom1, atom2, which, element, subtests):
        """Check that two atoms have the same requested displacements."""
        _test_str = f'{atom1}, {element}, '
        displ_1 = atom1.displacements[which]
        displ_2 = atom2.displacements[which]
        with subtests.test(_test_str + 'total'):
            assert displ_1[element] == pytest.approx(displ_2[element])

        offs_1 = displ_1.vibrocc_offset
        offs_2 = displ_2.vibrocc_offset
        with subtests.test(_test_str + 'vibrocc'):
            assert offs_1[element] == pytest.approx(offs_2[element])

        if which == 'geo':  # Check also offset
            offs_1 = displ_1.displacements_offset
            offs_2 = displ_2.displacements_offset
            with subtests.test(_test_str + 'geo offset'):
                assert offs_1[element] == pytest.approx(offs_2[element])


class TestRevertUnitCell:
    """Tests for reverting the unit cell of a slab."""

    @todo
    def test_revert_unit_cell(self):                                            # TODO: Probably best to pick a few random operations and make sure that reverting one+rest, a few+rest, or all of them at once gives the same result. This should include unit cell as well as all atom frac and cart coordinates
        """TODO"""

    @staticmethod
    def check_identical(slab, other, subtests):
        """Check that slab and other are identical."""
        with subtests.test('ucell'):
            assert slab.ucell == pytest.approx(other.ucell)
        atom_pairs = zip(slab, other)
        with subtests.test('atom elements'):
            assert all(at.el == at2.el for at, at2 in atom_pairs)
        with subtests.test('atom pos'):
            assert all(np.allclose(at.pos, at2.pos) for at, at2 in atom_pairs)
        with subtests.test('atom cartpos'):
            assert all(np.allclose(at.cartpos, at2.cartpos)
                       for at, at2 in atom_pairs)

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_one_operation(self, args, subtests):
        """Check correct result of reverting one unit-cell operation."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        slab.rotate_unit_cell(6)
        slab.revert_unit_cell()
        self.check_identical(slab, slab_copy, subtests)

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_few_operation(self, args, subtests):
        """Same as above, but reverting a few operations."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        slab.rotate_unit_cell(6)
        slab.rotate_unit_cell(4)
        slab.rotate_unit_cell(8)                                                # TODO: so far, this is the only type of operation we have built in
        slab.revert_unit_cell()
        self.check_identical(slab, slab_copy, subtests)

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_nothing_to_undo(self, args, subtests):                             # TODO: both by having nothing to undo, and by passing as many as there are operations. Check especially by manually translating atoms out of the base cell. – @michele-riva: not sure what you mean by this
        """Check that reverting with no operations does nothing."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        slab.revert_unit_cell()
        self.check_identical(slab, slab_copy, subtests)

    def test_raises(self, ag100):
        """Check complaints for invalid operation types."""
        slab, *_ = ag100
        slab.ucell_mod.append(('invalid', 'invalid'))
        with pytest.raises(RuntimeError):
            slab.revert_unit_cell()


@todo
class TestSlabLayers:
    """Collection of tests concerning slab (sub)layers."""

    def test_bulk_layers(self):
        """TODO"""

    def test_create_layers(self):                                               # TODO: check also logging with cuts that (do not) create empty layers
        """Check that layers are created correctly."""

    def test_create_sublayers(self):                                            # TODO: also test if this works fine excluding the second sort-by-element run
        """Check that sublayers are created correctly."""

    def test_full_update_with_layers_defined(self):                             # TODO: test correct behaviour for (i) coords initially outside the unit cell, and (ii) no topat_ori_z available
        """TODO"""

    def test_interlayer_spacing(self):                                          # TODO: also raises.
        """TODO"""

    def test_slab_lowocc_sublayer(self):
        """TODO"""


class TestSlabRaises:
    """Collection of tests for diverse exception-raising conditions."""

    def test_ab_in_plane_raises(self):
        """Check complaints when accessing an undefined unit cell."""
        with pytest.raises(err.InvalidUnitCellError):
            Slab().check_a_b_in_plane()

    _sublayers = {  # attr_to_clear, exception
        'no atoms': ('atlist', err.EmptySlabError),
        'no elements': ('n_per_elem', err.MissingElementsError),
        }

    @parametrize('attr,exc', _sublayers.values(), ids=_sublayers)
    def test_create_sublayers_raises(self, attr, exc, ag100):
        """Check complaints when creating sublayers."""
        slab, *_ = ag100
        getattr(slab, attr).clear()
        with pytest.raises(exc):
            slab.create_sublayers()

    def test_from_slab_not_a_slab(self):
        """Check complaints when from_slab is passed a non-Slab."""
        with pytest.raises(TypeError):
            Slab.from_slab('invalid')

    def test_interlayer_spacing_few_layers(self, manual_slab_3_atoms):
        """Check complaints when there's not enough layers."""
        slab = manual_slab_3_atoms
        with pytest.raises(err.TooFewLayersError):
            _ = slab.smallest_interlayer_spacing

    _props = {
        'ab_cell': err.InvalidUnitCellError,
        'fewest_atoms_sublayer': err.MissingSublayersError,
        'smallest_interlayer_spacing': err.MissingLayersError,
        }

    @parametrize('attr_name,exc', _props.items(), ids=_props)
    def test_property_raises(self, attr_name, exc):
        """Check complaints when accessing a property of an empty slab."""
        with pytest.raises(exc):
            _ = getattr(Slab(), attr_name)


class TestSorting:
    """Collection of tests for slab sorting."""

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_element_sort(self, args, shuffle_slab):
        """Check correct element-based sorting of a Slab."""
        slab, *_ = args
        shuffle_slab(slab)
        slab.sort_by_element()
        element_index_orig_order = [slab.elements.index(at.el) for at in slab]
        assert element_index_orig_order == sorted(element_index_orig_order)

    def test_element_sort_raises_with_outdated_elements(self, ag100):
        """Ensure sorting complains when elements are outdated."""
        slab, *_ = ag100
        new_atom = slab.atlist[0].duplicate()
        new_atom.el = 'C'
        with pytest.raises(err.SlabError):
            slab.sort_by_element()

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_sort_original(self, args, shuffle_slab):
        """Check correct sorting of a Slab to original atom numbers."""
        slab, *_ = args
        original_atlist = deepcopy(slab.atlist)
        shuffle_slab(slab)
        slab.sort_original()
        assert all(all(at1.pos == at2.pos) and at1.el == at2.el
                   for at1, at2 in zip(slab, original_atlist))

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    def test_simple_sort_by_z(self, args, shuffle_slab):
        """Check correct outcome of sorting by z coordinates."""
        slab, *_ = args
        shuffle_slab(slab)
        slab.sort_by_z()
        assert all(at1.pos[2] <= at2.pos[2] for at1, at2 in pairwise(slab))

    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar)
    @parametrize(bottom_to_top=(True, False))
    def test_z_sort(self, args, bottom_to_top):
        """Check successful sorting of atoms by out-of-plane position."""
        slab, *_ = args
        slab.sort_by_z(bottom_to_top=bottom_to_top)
        _ordered = operator.ge if bottom_to_top else operator.le                # TODO: swap when flipping .cartpos
        assert all(_ordered(at1.cartpos[2], at2.cartpos[2])
                   for at1, at2 in pairwise(slab))


class TestSuperAndSubCell:
    """Collection of tests for creation of larger and smaller slab versions."""

    invalid = {
        'not integer': (np.eye(2)*0.5, NonIntegerMatrixError),
        'singular': (np.eye(2)*0, SingularMatrixError),
        'non-2x2': (np.eye(3), ValueError),
        }

    @todo
    def test_supercell_valid(self):                                             # TODO: diagonal and non-diagonal (for some weird basis cell?). I think the old version was failing under some non-diagonal situations. Explicitly test the two removed update_origin.
        """TODO"""

    @parametrize('matrix,exc', invalid.values(), ids=invalid)
    def test_supercell_invalid(self, matrix, exc, ag100):
        """Check complaints for invalid input to make_subcell."""
        slab, *_ = ag100
        with pytest.raises(exc):
            slab.make_supercell(matrix)

    sub_invalid = {
        **invalid,
        'not a subcell': (np.diag((2, 2)), ValueError),
        }

    @todo
    def test_subcell_valid(self):
        """TODO"""

    @parametrize('matrix,exc', sub_invalid.values(), ids=sub_invalid)
    def test_subcell_invalid(self, matrix, exc, ag100):
        """Check complaints for invalid input to make_subcell."""
        slab, rpars, *_ = ag100
        with pytest.raises(exc):
            slab.make_subcell(rpars, matrix)


class TestUnitCellTransforms:
    """Test simple transformations of the unit cell of a slab."""

    def test_angle_between_ucell_and_coord_sys_zero(self, manual_slab_3_atoms):
        """Check correct identification of the rotation of the a vector."""
        slab = manual_slab_3_atoms
        assert slab.angle_between_ucell_and_coord_sys == pytest.approx(0)

    def test_angle_between_ucell_and_coord_sys_30(self):
        """Check correct identification of the rotation of the a vector."""
        slab = Slab()
        slab.ucell = np.array([[np.cos(np.pi/6), np.sin(np.pi/6), 0],
                               [np.cos(np.pi/3*4), np.sin(np.pi/3*4), 0],
                               [0, 0, 1]]).T
        assert slab.angle_between_ucell_and_coord_sys == pytest.approx(30)

    _scalings = {  # Applied, expected result
        'scalar': ((2,),  np.diag((6, 8, 10))),
        'vector': ((1/3, 3.14, 0.1), np.diag((1, 12.56, 0.5)))
        }
    _scale_invalid = {
        'too few': (tuple(), TypeError),
        'two': ((1, 2), TypeError),
        'too many': ((1, 2, 3, 4), TypeError),
        'not a number': (('invalid',), TypeError),
        'singular': ((0, 1, 1), ValueError),
        }

    @parametrize('scaling,expected', _scalings.values(), ids=_scalings)
    def test_apply_scaling(self, scaling, expected, manual_slab_3_atoms):
        """Check expected outcome of scaling the unit cell."""
        slab = manual_slab_3_atoms
        fractionals = [at.pos.copy() for at in slab]
        slab.apply_scaling(*scaling)
        assert slab.ucell == pytest.approx(expected)
        assert all(np.allclose(at.pos, ori_pos)
                   for at, ori_pos in zip(slab, fractionals))

    @parametrize('scaling,exc', _scale_invalid.values(), ids=_scale_invalid)
    def test_apply_scaling_raises(self, scaling, exc, ag100):
        """Check complaints for invalid arguments to apply_scaling."""
        slab, *_ = ag100
        with pytest.raises(exc):
            slab.apply_scaling(*scaling)

    _invalid_matrix = {
        'not a matrix': ('invalid', ValueError),
        'not orthogonal': (np.identity(3)+1, ValueError),
        }

    @parametrize('matrix,exc', _invalid_matrix.values(), ids=_invalid_matrix)
    def test_matrix_transform_raises(self, matrix, exc, ag100):
        """Check complaints for unexpected inputs."""
        slab, *_ = ag100
        with pytest.raises(exc):
            slab.apply_matrix_transformation(matrix)

    def test_project_c_to_z(self, make_poscar, subtests):
        """Check that the c vector is parallel to the z axis after projection."""
        slab, *_ = make_poscar(poscar_slabs.SLAB_36C_cm)
        slab_copy = deepcopy(slab)
        slab.project_c_to_z()
        with subtests.test('c along z'):
            assert np.allclose(slab.ucell.T[2][:2], 0)
        with subtests.test('atom positions'):
            assert all(  # atoms 4,5 are wrapped around
                np.allclose(at.cartpos, at_copy.cartpos)
                for at, at_copy in zip(slab.atlist[~4:5], slab_copy.atlist[~4:5])
                )

    def test_rotation_on_trigonal_slab(self, manual_slab_1_atom_trigonal):
        """Test application of a rotation to a trigonal slab."""
        rot_15 = [[0.96592583, -0.25881905,  0.],
                  [0.25881905,  0.96592583,  0.],
                  [0.,  0.,  1.]]
        expected_cell = [[0.96592583,  0.25881905,  0.],
                         [-2.99808654, 2.30249368,  0.],
                         [0.44828774,  2.1906707,   3.]]
        expected_atom_cartpos = [-1.86064664,  1.88257645]
        slab = manual_slab_1_atom_trigonal
        slab.apply_matrix_transformation(rot_15)
        assert np.allclose(slab.ucell.T, expected_cell)
        assert np.allclose(slab.atlist[0].cartpos[:2], expected_atom_cartpos)

    def test_ucell_array(self, manual_slab_3_atoms):
        """Check that the array shape of the unit cell is as expected."""
        assert manual_slab_3_atoms.ucell.shape == (3, 3)


@todo
class TestUnitCellReduction:
    """Tests for minimization of various bits of the unit cell."""

    def test_ab_cell_minimization(self):                                        # TODO: Use also "Sb on Si(111)" case from Max Buchta
        """TODO"""

    def test_ab_cell_already_minimal(self):
        """TODO"""


@todo
def test_contains():
    """TODO"""


@todo
def test_check_ab_in_plane():
    """TODO"""


@todo
def test_nearest_neighbors():
    """TODO"""


def test_ucell_ori_after_clear_symmetry_and_ucell_history(ag100):
    """Ensure modifications to the unit cell don't go into ucell_ori."""
    slab, *_ = ag100
    assert slab.celltype == 'unknown'
    ucell_ori = slab.ucell.copy()
    slab.clear_symmetry_and_ucell_history()
    slab.ucell *= 3
    assert np.allclose(ucell_ori, slab.ucell_ori)


def test_translation_symmetry_different_species(make_poscar):
    """Check that an MgO slab is not equivalent upon Mg->O translation."""
    slab, rpars, *_ = make_poscar(poscar_slabs.SLAB_MgO)
    mg_to_oxygen = slab.ab_cell.T[0] / 2
    assert not slab.is_translation_symmetric(mg_to_oxygen, rpars.SYMMETRY_EPS)


@pytest.mark.xfail(reason='Issue #140')
def test_layer_cutting_for_slab_with_incomplete_bulk_layer(make_poscar):
    """Test for issue #140."""
    slab, rpars, *_ = make_poscar(poscar_slabs.SLAB_Cu2O_111)
    rpars.BULK_LIKE_BELOW = 0.35
    slab.detect_bulk(rpars)
    assert slab.smallest_interlayer_spacing >= 1.0
