"""Tests for viperleed.calc.classes.slab."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-07-28'
__license__ = 'GPLv3+'

from collections import Counter
from copy import deepcopy
import logging
import operator
from random import shuffle

import numpy as np
import pytest
from pytest_cases import fixture, fixture_ref
from pytest_cases import parametrize, parametrize_with_cases
from scipy.spatial.distance import cdist as euclid_distance

from viperleed.calc.classes.atom import Atom
from viperleed.calc.classes.atom_containers import AtomList
from viperleed.calc.classes.rparams import LayerCuts
from viperleed.calc.classes.rparams import Rparams
from viperleed.calc.classes.slab import BulkSlab
from viperleed.calc.classes.slab import Slab
from viperleed.calc.classes.slab import errors as err
from viperleed.calc.classes.slab import surface_slab
from viperleed.calc.classes.sym_entity import SymPlane
from viperleed.calc.files.parameters.errors import InconsistentParameterError
from viperleed.calc.lib import leedbase
from viperleed.calc.lib.coordinates import add_edges_and_corners
from viperleed.calc.lib.coordinates import collapse
from viperleed.calc.lib.itertools_utils import pairwise
from viperleed.calc.lib.matrix import NonIntegerMatrixError
from viperleed.calc.lib.matrix import SingularMatrixError

from ...helpers import exclude_tags, not_raises
from .. import cases_ase, poscar_slabs
from ..tags import CaseTag as Tag

CasePOSCARSlabs = poscar_slabs.CasePOSCARSlabs
todo = pytest.mark.skip('To be implemented')
if_ase = pytest.mark.skipif(
    not surface_slab._HAS_ASE,  # pylint: disable=protected-access
    reason='No ASE module'
    )
with_bulk_repeat = parametrize_with_cases('args', cases=CasePOSCARSlabs,
                                          has_tag=Tag.BULK_PROPERTIES)
with_layers = parametrize_with_cases(
    'args',
    cases=CasePOSCARSlabs.case_layer_info_poscar
    )
infoless_poscar = parametrize_with_cases(
    'args',
    cases=CasePOSCARSlabs.case_infoless_poscar
    )


@fixture(name='shuffle_slab', scope='session')
def make_shuffled_slab():
    """Shuffle the atoms of a slab at random."""
    def _shuffle(slab):
        slab.atlist.strict = False  # Silence duplicate errors
        shuffle(slab.atlist)
        slab.atlist.strict = True   # Now all atoms should be unique
        slab.atlist.update_atoms_map()
    return _shuffle


@fixture(name='check_identical')
def factory_check_identical(subtests):
    """Check that slab and other are identical."""
    def _get_frac_and_cart(slab):
        frac = np.array([at.pos for at in slab])
        cart = frac.dot(slab.ucell.T)
        return frac, cart

    def check_identical(slab, other):
        with subtests.test('n_atoms'):
            assert slab.n_atoms == other.n_atoms
        with subtests.test('ucell'):
            assert slab.ucell == pytest.approx(other.ucell)
        with subtests.test('atom elements'):
            assert all(at.el == at2.el for at, at2 in zip(slab, other))

        # Compare coordinates including replicas for edge/corner atoms
        eps = 1e-4
        frac_slab, cart_slab = _get_frac_and_cart(slab)
        frac_other, cart_other = _get_frac_and_cart(other)
        cart_other, _ = add_edges_and_corners(cart_other, frac_other,
                                              (eps,)*3, other.ucell.T)
        frac_other = cart_other.dot(np.linalg.inv(other.ucell.T))
        with subtests.test('atom pos'):
            frac_dist = euclid_distance(frac_slab, frac_other).min(axis=1)
            assert all(frac_dist <= eps)
        with subtests.test('atom cartpos'):
            cart_dist = euclid_distance(cart_slab, cart_other).min(axis=1)
            assert all(cart_dist <= eps)
    return check_identical


@fixture(name='with_one_thick_bulk', scope='session')
def factory_with_one_thick_bulk():
    """Prepare slab and rpars to have a one-layer thick bulk slab."""
    def _make(slab, rpars, info):
        cut = info.bulk_properties.bulk_like_below
        rpars.LAYER_CUTS.update_from_sequence([cut])
        rpars.N_BULK_LAYERS = 1
        rpars.BULK_REPEAT = None
        slab.create_layers(rpars)

        # Shift the bottom atom to c=0 to prevent back-folding
        # of the top ones when creating the thick bulk
        c_shift = slab.remove_vacuum_at_bottom(rpars)
        slab.make_bulk_slab(rpars, recenter=False)
        return c_shift
    return _make


class TestABInPlane:
    """Tests for checking that the first two unit vectors have no z."""

    def test_not_in_plane(self, ag100):
        """Check complaints when one unit vector has z component."""
        slab, *_ = ag100
        slab.ucell.T[0, 2] = 1
        with pytest.raises(err.InvalidUnitCellError):
            slab.check_a_b_in_plane()

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_successful(self, args):
        """Check that a slab has no z component of a & b."""
        slab, *_ = args
        with not_raises(err.InvalidUnitCellError):
            slab.check_a_b_in_plane()

    def test_undefined_cell_raises(self):
        """Check complaints when accessing an undefined unit cell."""
        with pytest.raises(err.InvalidUnitCellError):
            Slab().check_a_b_in_plane()


class TestAtlistIsNotList:
    """Test for methods that may potentially transform .atlist into a list."""

    def _check_not_list(self, slab):
        """Check that slab.atlist is not a list."""
        assert isinstance(slab.atlist, AtomList)
        assert not isinstance(slab.atlist, list)
        assert not isinstance(slab.atlist, tuple)

    def test__add_one_bulk_cell(self, ag100):
        """Check that _add_one_bulk_cell does not mess with atlist."""
        slab, rpars, *_ = ag100
        # Easier to use with_extra_bulk_units which calls the
        # _add_one_bulk_cell method rather than fiddling with
        # all the input arguments to _add_one_bulk_cell
        thicker, _ = slab.with_extra_bulk_units(rpars, 5)
        self._check_not_list(slab)
        self._check_not_list(thicker)

    def test_create_layers(self, ag100):
        """Check that create_layers does not mess with atlist."""
        slab, rpars, *_ = ag100
        slab.create_layers(rpars)
        self._check_not_list(slab)

    @if_ase
    @parametrize_with_cases('args', cases=cases_ase)
    def test_from_ase(self, args):
        """Check that from_ase does not mess with atlist."""
        slab = Slab.from_ase(args[0])
        self._check_not_list(slab)

    def test_make_bulk_slab(self, ag100):
        """Check that make_bulk_slab does not mess with atlist."""
        slab, rpars, *_ = ag100
        bulk = slab.make_bulk_slab(rpars)
        self._check_not_list(slab)
        self._check_not_list(bulk)

    def test_remove_duplicate_atoms(self, ag100):
        """Check that remove_duplicate_atoms does not mess with atlist."""
        slab, rpars, *_ = ag100
        slab.atlist[0].duplicate()
        slab.remove_duplicate_atoms(rpars.SYMMETRY_EPS, rpars.SYMMETRY_EPS.z)
        self._check_not_list(slab)


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

    def test_mirror_on_slanted_cell(self, subtests):
        """Test the expected outcome of mirroring atoms of a slanted slab."""
        slab, *_ = CasePOSCARSlabs().case_poscar_36carbon_atoms_cm()
        slab.create_sublayers(0.1)
        with subtests.test('valid mirror plane'):
            sym_plane = SymPlane((0, 0), (0, 1), abt=slab.ab_cell.T)
            assert slab.is_mirror_symmetric(sym_plane, eps=1e-6)
        with subtests.test('not a glide'):
            sym_plane = SymPlane((0, 0), (1, 0), abt=slab.ab_cell.T,
                                 ty='glide')
            assert not slab.is_mirror_symmetric(sym_plane, eps=1e-6)

    def test_rotation_symmetric_raises(self, ag100):
        """Check that no symmetry checks can be performed without sublayers."""
        slab, *_ = ag100
        slab.sublayers = ()
        with pytest.raises(err.MissingSublayersError):
            slab.is_rotation_symmetric((0, 0), 4, 1e-3)

    def test_180_rotation(self, manual_slab_3_atoms):
        """Test the expected outcome of rotating atoms of a simple slab."""
        slab = manual_slab_3_atoms
        rotated_slab = deepcopy(slab)
        rotated_slab.rotate_atoms(order=2)
        assert all(at.is_same_xy(rot_at)
                   for at, rot_at in zip(slab, reversed(rotated_slab)))

    def test_translate_c_clears_layers(self, ag100):
        """Check that no (sub)layers are present after a c shift."""
        slab, *_ = ag100
        slab.create_sublayers(0.1)
        assert slab.layers
        slab.translate_atoms_c(0.15)
        assert not slab.layers
        assert not slab.sublayers


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
    def test_update_atom_numbers(self, make_bulk, in_bulk,
                                 to_remove, subtests):
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
        with subtests.test('no. atoms'):
            assert fe3o4.n_per_elem['Fe'] == n_fe - 1
            assert fe3o4.n_per_elem['O'] == n_oxygen
        new_at_nrs = [at.num for at in fe3o4]
        with subtests.test('atom.num'):
            assert new_at_nrs == list(range(1, n_fe + n_oxygen))
        if not make_bulk:
            return
        with subtests.test('no. atoms, bulk'):
            assert bulk.n_per_elem['O'] == n_oxygen_bulk
            assert bulk.n_per_elem['Fe'] == n_fe_bulk - (1 if in_bulk else 0)
        with subtests.test('atom.num, bulk'):
            assert all(at.num in new_at_nrs for at in bulk)
        with subtests.test('atom.num, bulk, changed'):
            assert [at.num for at in bulk] != old_at_nrs_bulk


class TestAtomsInMultipleCells:
    """Check correct behaviour of .has_atoms_in_multiple_c_cells()."""

    def test_all_at_edges(self):
        """Check correct outcome when all atoms are at cell edges."""
        slab, *_ = CasePOSCARSlabs().case_poscar_mgo()
        edge_atoms = [slab.atlist.get(num) for num in (1, 4, 17, 19)]
        other_atom = slab.atlist.get(12)
        for atom in (other_atom, *edge_atoms):
            atom.pos[2] += 3
        assert slab.has_atoms_in_multiple_c_cells()
        # Keep only atoms at c edges
        slab.atlist.clear()
        slab.atlist.extend(edge_atoms)
        assert not slab.has_atoms_in_multiple_c_cells()

    @parametrize_with_cases('args', cases=CasePOSCARSlabs)
    def test_all_in_first(self, args):
        """Check that all atoms belong to the same cell."""
        slab, *_ = args
        # Notice that they all are because we full_update when
        # creating the case, which collapses to the base cell
        assert not slab.has_atoms_in_multiple_c_cells()

    @parametrize(shift_cell=(1, -8, 27))
    def test_all_in_same(self, ag100, shift_cell):
        """Check that all atoms belong to the same cell (although not base)."""
        slab, *_ = ag100
        for atom in slab:
            atom.pos[2] = atom.pos[2] + shift_cell
        slab.update_cartesian_from_fractional()
        assert not slab.has_atoms_in_multiple_c_cells()

    def test_in_multiple_cells(self, ag100):
        """Check correct detection of atoms in multiple cells."""
        slab, *_ = ag100
        for atom in slab:
            shift_cell = np.random.randint(-5, 5)
            while not shift_cell:
                shift_cell = np.random.randint(-5, 5)
            atom.pos[2] = atom.pos[2] + shift_cell
        slab.update_cartesian_from_fractional()
        assert slab.has_atoms_in_multiple_c_cells()

    def test_in_same_except_edges(self):
        """Ensure that cell-belonging checks discard atoms at edges."""
        slab, *_ = CasePOSCARSlabs().case_poscar_mgo()
        atom_at_edge = slab.atlist.get(1)  # First atom at c=1e-5
        atom_at_edge.pos[2] += 7
        slab.update_cartesian_from_fractional()
        assert not slab.has_atoms_in_multiple_c_cells()

    def test_raises_with_no_atoms(self):
        """Check complaints if there's no atoms in a slab."""
        with pytest.raises(err.EmptySlabError):
            Slab().has_atoms_in_multiple_c_cells()


class TestBulk3DOperations:
    """Tests for 3D symmetry operations of BulkSlab objects."""

    @staticmethod
    def move_atom(slab, atom_index, distance):
        """Move the atom at a given index by distance in a random direction."""
        atom = slab.atlist.get(atom_index)
        direction = np.random.randn(2)
        direction /= np.linalg.norm(direction)
        atom.cartpos[:2] += distance * direction
        slab.update_fractional_from_cartesian()

    @parametrize_with_cases('bulk', cases=poscar_slabs, has_tag=Tag.BULK)
    def test_get_candidate_layer_periods(self, bulk):
        """Check correct identification of sublayer periods."""
        slab, rpars, info = bulk
        if not slab.sublayers:
            slab.create_sublayers(rpars.SYMMETRY_EPS.z)
        periods = slab.get_candidate_layer_periods(rpars.SYMMETRY_EPS.z)
        assert periods == info.bulk.periods

    @parametrize_with_cases('bulk', cases=poscar_slabs.case_double_bulk)
    def test_get_candidate_layer_periods_raises(self, bulk):
        """Check correct identification of sublayer periods."""
        slab, rpars, *_ = bulk
        with pytest.raises(err.MissingSublayersError):
            slab.get_candidate_layer_periods(rpars.SYMMETRY_EPS.z)

    def test_get_candidate_layer_periods_thin(self, ag100):
        """Check that a slab with few layers has no candidate periods."""
        slab, rpars, *_ = ag100
        bulk = slab.make_bulk_slab(rpars)
        assert not bulk.get_candidate_layer_periods(rpars.SYMMETRY_EPS.z)

    @fixture(name='bulk_slab_with_glide')
    def fixture_bulk_slab_with_glide(self):
        """Return a very basic BulkSlab with a 3D glide plane."""
        slab = BulkSlab()
        slab.ucell = 4.0*np.diag((1, 1, 2))
        slab.atlist.extend((
            Atom('Fe', np.array([0.0, 0, 0.0]), 1, slab),
            Atom('Fe', np.array([0.5, 0, 0.5]), 2, slab),
            Atom('C', np.array([0.25, 0, 0.0]), 3, slab),  # on glide
            Atom('C', np.array([0.25, 0, 0.5]), 4, slab),  # on glide
            ))
        slab.update_element_count()
        slab.update_cartesian_from_fractional()
        eps = 0.1
        slab.create_sublayers(eps)
        ab_cell = slab.ab_cell.T
        glide = SymPlane(np.array([0.25, 0]), ab_cell[1], ab_cell)
        return slab, glide, eps

    @pytest.mark.xfail(reason='Known bug. Will be fixed in better-symmetry',
                       strict=False)  # Sometimes they succeed
    @parametrize(atoms_to_move=((1,), (3,), (1, 3)))
    def test_bulk_glide_symmetric_atom_moved(self, atoms_to_move,
                                             bulk_slab_with_glide):
        """Check identification of bulk glide plane with atoms a bit off."""
        slab, glide, eps = bulk_slab_with_glide
        for atom_ind in atoms_to_move:
            self.move_atom(slab, atom_ind, 0.99*eps)
        assert slab.is_bulk_glide_symmetric(glide, 2, eps)

    fe3o4_bulk = poscar_slabs.CaseBulkSlabs.case_fe3o4_bulk

    @parametrize_with_cases('bulk', cases=fe3o4_bulk)
    def test_bulk_glide_symmetric(self, bulk):
        """Check correct identification of bulk glide plane."""
        slab, rpars, *_ = bulk
        ab_cell = slab.ab_cell.T
        glide = SymPlane(np.array([0, 0]), ab_cell[0] + ab_cell[1], ab_cell)
        assert slab.is_bulk_glide_symmetric(glide, 3, rpars.SYMMETRY_EPS)

    @parametrize_with_cases('bulk', cases=fe3o4_bulk)
    def test_bulk_screw_symmetric(self, bulk, subtests):
        """Check correct identification of bulk screw axis."""
        slab, rpars, *_ = bulk
        kwargs = {'order': 4, 'sublayer_period': 3, 'eps': rpars.SYMMETRY_EPS}
        with subtests.test('Move atom that is transformed by 4-fold screw'):
            # Atom no. 49 (oxygen) is not on the 4-fold screw
            self.move_atom(slab, 49, 0.99*rpars.SYMMETRY_EPS)
            assert slab.is_bulk_screw_symmetric(**kwargs)
        with subtests.test('Move atom at 4-fold screw'):
            # Atom no. 33 (tetrahedral Fe) is on the 4-fold screw
            self.move_atom(slab, 33, 0.99*rpars.SYMMETRY_EPS)
            try:
                assert slab.is_bulk_screw_symmetric(**kwargs)
            except AssertionError:
                pytest.xfail('Known bug. Will be fixed in better-symmetry')
            else:
                raise AssertionError('XPASS: This test should FAIL unless '
                                     'this is the better-symmetry branch')
        with subtests.test('Move both atoms'):
            self.move_atom(slab, 33, 0.99*rpars.SYMMETRY_EPS)
            self.move_atom(slab, 49, 0.99*rpars.SYMMETRY_EPS)
            try:
                assert slab.is_bulk_screw_symmetric(**kwargs)
            except AssertionError:
                pytest.xfail('Known bug. Will be fixed in better-symmetry')
            else:
                raise AssertionError('XPASS: This test should FAIL unless '
                                     'this is the better-symmetry branch')


class TestBulkDetect:
    """Collection of tests for detecting bulk below a cutoff."""

    @staticmethod
    def prepare_to_detect(slab, rpars, bulk_like_below):
        """Prepare slab and rpars to detect_bulk."""
        rpars.BULK_LIKE_BELOW = bulk_like_below
        rpars.BULK_REPEAT = None
        slab.create_layers(rpars)
        slab.create_sublayers(rpars.SYMMETRY_EPS.z)

    @with_bulk_repeat
    def test_detect_bulk(self, args, subtests):
        """Test function for detecting bulk cuts and distances."""
        slab, rpars, info = args
        bulk_info = info.bulk_properties
        self.prepare_to_detect(slab, rpars, bulk_info.bulk_like_below)

        bulk_cuts, bulk_dist = slab.detect_bulk(rpars)
        with subtests.test('no. bulk cuts'):
            assert len(bulk_cuts) == len(bulk_info.bulk_cuts)
        with subtests.test('bulk cut values'):
            assert np.allclose(bulk_cuts, bulk_info.bulk_cuts, atol=1e-4)
        with subtests.test('distance between bulk layers'):
            assert bulk_dist == pytest.approx(bulk_info.bulk_dist, abs=1e-4)

    @with_bulk_repeat
    def test_detect_bulk_fail_leaves_rp_sl_unchanged(self, args,
                                                     check_identical):
        """Test function for detecting bulk cuts and distances."""
        slab, rpars, *_ = args
        self.prepare_to_detect(slab, rpars, 0.01)

        sl_copy, rp_copy = deepcopy(slab), deepcopy(rpars)
        with pytest.raises((err.TooFewLayersError, err.NoBulkRepeatError)):
            slab.detect_bulk(rpars)

        # Check that slab and rpars are unchanged
        check_identical(slab, sl_copy)

        assert rpars.BULK_LIKE_BELOW == rp_copy.BULK_LIKE_BELOW
        assert str(rpars.LAYER_CUTS) == str(rp_copy.LAYER_CUTS)
        assert rpars.BULK_REPEAT == rp_copy.BULK_REPEAT
        assert np.allclose(rpars.SUPERLATTICE, rp_copy.SUPERLATTICE)

    def test_detect_bulk_raises_negative_bulk_like_below(self, ag100):
        """Check complaints for negative BULK_LIKE_BELOW."""
        slab, rpars, *_ = ag100
        rpars.BULK_LIKE_BELOW = -0.5
        with pytest.raises(ValueError):
            slab.detect_bulk(rpars)

    def test_detect_bulk_raises_with_defined_bulk_repeat(self, ag100):
        """Check complaints for existing BULK_REPEAT."""
        slab, rpars, *_ = ag100
        rpars.BULK_LIKE_BELOW = 0.65
        rpars.BULK_REPEAT = np.array([1, 1, 1])
        with pytest.raises(ValueError):
            slab.detect_bulk(rpars)

    @pytest.mark.xfail(reason='Issue #140')
    def test_layer_cutting_for_slab_with_incomplete_bulk_layer(self,
                                                               make_poscar):
        """Test for issue #140."""
        slab, rpars, *_ = make_poscar(poscar_slabs.SLAB_Cu2O_111)
        rpars.BULK_LIKE_BELOW = 0.35
        slab.detect_bulk(rpars)
        assert slab.smallest_interlayer_gap >= 1.0


class TestBulkRepeat:
    """Collection of test for bulk-repeat finding and returning."""

    @with_bulk_repeat
    def test_get(self, args):
        """Test get_bulk_repeat method."""
        slab, rpars, info = args
        bulk_info = info.bulk_properties
        TestBulkDetect.prepare_to_detect(slab, rpars,
                                         bulk_info.bulk_like_below)
        slab.detect_bulk(rpars)
        repeat_vector = slab.get_bulk_repeat(rpars)
        repeat_vector_bulk = slab.bulkslab.get_bulk_repeat(rpars)
        atol = float(0.2*rpars.SYMMETRY_EPS)
        assert repeat_vector == pytest.approx(bulk_info.bulk_repeat, abs=atol)
        assert repeat_vector_bulk == pytest.approx(repeat_vector)

    def test_get_raises_no_bulk_layers(self, ag100):
        """Check complaints when get_bulk_repeat is called without bulk."""
        slab, rpars, *_ = ag100
        rpars.BULK_REPEAT = None
        slab.layers = ()
        with pytest.raises(err.TooFewLayersError) as exc:
            slab.get_bulk_repeat(rpars)
        assert exc.match('no bulk')

    def test_get_raises_only_bulk_layers(self, ag100):
        """Check complaints when get_bulk_repeat is called without non-bulk."""
        slab, rpars, *_ = ag100
        slab.make_bulk_slab(rpars)
        bulk_twice = slab.bulkslab.with_double_thickness()
        slab_with_ony_bulk = slab.from_slab(bulk_twice)
        rpars.BULK_REPEAT = None
        rpars.N_BULK_LAYERS = 2
        slab_with_ony_bulk.create_layers(rpars)
        with pytest.raises(err.TooFewLayersError) as exc:
            slab_with_ony_bulk.get_bulk_repeat(rpars)
        assert exc.match('only bulk')

    def test_get_returns_stored_bulk_repeat(self, ag100):
        """Test get_bulk_repeat returns stored bulk repeat if available."""
        slab, rpars, *_ = ag100
        rpars.BULK_REPEAT = np.array([1, 2, 3])
        assert np.allclose(slab.get_bulk_repeat(rpars), np.array([1, 2, 3]))

    @with_bulk_repeat
    def test_get_returns_z_vector(self, args):
        """Test get_bulk_repeat gives a z-only vector without BULK_REPEAT."""
        slab, rpars, info = args
        rpars.BULK_REPEAT = None
        atol = float(0.2*rpars.SYMMETRY_EPS)
        assert np.allclose(slab.get_bulk_repeat(rpars),
                           [0, 0, info.bulk_properties.bulk_repeat[2]],
                           atol=atol)

    @with_bulk_repeat
    def test_identify(self, args, subtests):
        """Tests identify_bulk_repeat method."""
        slab, rpars, info = args
        bulk_info = info.bulk_properties
        rpars.BULK_REPEAT = None
        slab.make_bulk_slab(rpars)
        repeat_vector = slab.identify_bulk_repeat(eps=rpars.SYMMETRY_EPS,
                                                  epsz=rpars.SYMMETRY_EPS.z)
        atol = float(0.2*rpars.SYMMETRY_EPS)
        with subtests.test('z component'):
            assert repeat_vector[2] == pytest.approx(bulk_info.bulk_repeat[2],
                                                     abs=atol)
        # Compare unit-cell-collapsed versions of the bulk repeat
        # vector to get the in-plane components right. This way we
        # maintain the test independent of the exact atom in the
        # sublayers that is used to determine the repeat vector.
        # E.g., in orthorhombic cells one may get a (also acceptable)
        # opposite sign for the in-plane components.
        collapsed_repeat, _ = collapse(repeat_vector, slab.ucell.T)
        collapsed_expected, _ = collapse(bulk_info.bulk_repeat, slab.ucell.T)
        with subtests.test('collapsed'):
            assert collapsed_repeat == pytest.approx(collapsed_expected,
                                                     abs=atol)

    _invalid_cuts = {
        'too few layers': 0.41,
        'no repeat vector': 0.32,
        'elements mismatched': 0.30,
        }

    @parametrize(cut=_invalid_cuts.values(), ids=_invalid_cuts)
    def test_identify_raises(self, cut, with_one_thick_bulk):
        """Check complaints if there's not enough layers for identification."""
        slab, rpars, info = CasePOSCARSlabs().case_poscar_tio2_small()
        info.bulk_properties.bulk_like_below = cut
        with_one_thick_bulk(slab, rpars, info)
        with pytest.raises(err.NoBulkRepeatError):
            slab.identify_bulk_repeat(eps=rpars.SYMMETRY_EPS,
                                      epsz=rpars.SYMMETRY_EPS.z)

    def test_identify_raises_without_bulkslab(self, ag100):
        """Check complaints when called without a bulk slab."""
        slab, *_ = ag100
        slab.bulkslab = None
        with pytest.raises(err.MissingBulkSlabError):
            slab.identify_bulk_repeat(.1)


class TestBulkUcell:
    """Tests concerning reduction of bulk unit cell and C vector."""

    @with_bulk_repeat
    def test_get_min_c(self, args, with_one_thick_bulk):
        """Test that get_minimal_c_vector works as expected."""
        slab, rpars, info = args
        # strip out existing layers and make a single "fake" cut at
        # BULK_LIKE_BELOW. This gives us a bulk slab with (at least)
        # two bulk layers, which is what we need to test get_min_c.
        with_one_thick_bulk(*args)
        bulk_slab = slab.bulkslab
        bulk_slab.create_sublayers(rpars.SYMMETRY_EPS.z)
        # z_periodic=False because we use the original unit cell
        min_c = bulk_slab.get_minimal_c_vector(rpars.SYMMETRY_EPS,
                                               rpars.SYMMETRY_EPS.z,
                                               z_periodic=False)
        atol = float(0.2*rpars.SYMMETRY_EPS)
        assert min_c == pytest.approx(info.bulk_properties.bulk_repeat,
                                      abs=atol)

    @with_bulk_repeat
    def test_ensure_min_c(self, args, with_one_thick_bulk, subtests):
        """Test that ensure_minimal_c_vector works as expected."""
        slab, rpars, info = args
        with_one_thick_bulk(*args)
        bulk_slab = slab.bulkslab
        bulk_slab.ensure_minimal_c_vector(rpars)
        c_vec = bulk_slab.c_vector
        bulk_info = info.bulk_properties
        with subtests.test('bulk unit-cell c'):
            atol = float(0.2*rpars.SYMMETRY_EPS)
            assert c_vec == pytest.approx(bulk_info.bulk_repeat, abs=atol)
        with subtests.test('no. bulk atoms'):
            assert bulk_slab.n_atoms == bulk_info.n_bulk_atoms

    def test_ensure_min_c_raises_with_vacuum_below(self, ag100,
                                                   with_one_thick_bulk):
        """Check complaints of ensure_minimal_c_vector when there's vacuum."""
        c_shift = with_one_thick_bulk(*ag100)
        slab, rpars, *_ = ag100
        cart_shift = c_shift * slab.c_vector[2]
        cart_shift *= -1                                                        # TODO: remove with .cartpos[2] flip -- Issue #174
        bulk = slab.bulkslab
        for atom in bulk:
            atom.cartpos[2] += cart_shift  # Shift back upwards
        bulk.update_fractional_from_cartesian()
        with pytest.raises(err.SlabError) as exc:
            bulk.ensure_minimal_c_vector(rpars, z_periodic=False)
        assert exc.match(r'.*vacuum')

    @with_bulk_repeat
    def test_get_min_c_raises_already_minimal(self, args, with_one_thick_bulk):
        """Test that get_minimal_c_vector raises AlreadyMinimalError."""
        slab, rpars, *_ = args
        with_one_thick_bulk(*args)
        bulk_slab = slab.bulkslab
        bulk_slab.ensure_minimal_c_vector(rpars)
        with pytest.raises(err.AlreadyMinimalError):
            bulk_slab.get_minimal_c_vector(rpars.SYMMETRY_EPS)

    transform = {
        '2x2': np.diag((2, 2)),
        }
    recenter = {'recenter': True,
                'no recenter': False}

    @staticmethod
    def _check_centred(slab):
        """Ensure that atoms in slab are centred around cell midpoint in z."""
        atom_c_frac = [at.pos[2] for at in slab]
        midpos = (max(atom_c_frac) + min(atom_c_frac)) / 2
        assert midpos == pytest.approx(0.5, abs=0.02)

    @fixture(name='check_reduced_correctly')
    def factory_check_reduced_correctly(self, subtests):
        """Assert correct application of apply_bulk_cell_reduction."""
        def _check(reduced, expected_ucell, info, recenter):
            with subtests.test('ucell'):
                assert reduced.ucell == pytest.approx(expected_ucell)
            with subtests.test('no. atoms'):
                assert reduced.n_atoms == info.bulk_properties.n_bulk_atoms
            if recenter:
                with subtests.test('atoms are centred'):
                    self._check_centred(reduced)
        return _check

    @parametrize('recenter', recenter.values(), ids=recenter)
    @parametrize('transform', transform.values(), ids=transform)
    @with_bulk_repeat
    def test_apply_bulk_ucell_reduction_ab(self, recenter, transform, args,
                                           check_reduced_correctly):
        """Test that apply_bulk_cell_reduction works as expected for ab."""
        slab, rpars, info = args
        slab.make_bulk_slab(rpars)
        original_ab_cell = slab.bulkslab.ab_cell.T.copy()
        supercell = slab.make_supercell(transform)
        supercell_bulk = supercell.make_bulk_slab(rpars)
        supercell_bulk.apply_bulk_cell_reduction(eps=rpars.SYMMETRY_EPS,
                                                 new_ab_cell=original_ab_cell,
                                                 recenter=recenter)
        check_reduced_correctly(supercell_bulk, slab.bulkslab.ucell,
                                info, recenter)

    @parametrize('recenter', recenter.values(), ids=recenter)
    @with_bulk_repeat
    def test_apply_bulk_ucell_reduction_c(self, recenter, args,
                                          check_reduced_correctly):
        """Test that apply_bulk_cell_reduction works as expected for c."""
        slab, rpars, info = args
        bulkslab = slab.make_bulk_slab(rpars)
        original_ucell = bulkslab.ucell
        thick_slab = bulkslab.with_double_thickness()
        thick_slab.apply_bulk_cell_reduction(eps=rpars.SYMMETRY_EPS,
                                             new_c_vec=original_ucell.T[2],
                                             recenter=recenter,
                                             z_periodic=True)
        check_reduced_correctly(thick_slab, original_ucell, info, recenter)

    @with_bulk_repeat
    def test_apply_bulk_ucell_reduction_nothing(self, args,
                                                check_reduced_correctly):
        """Test that apply_bulk_cell_reduction does nothing."""
        slab, rpars, info = args
        bulkslab = slab.make_bulk_slab(rpars)
        slab_copy = deepcopy(bulkslab)
        slab_copy.apply_bulk_cell_reduction(eps=rpars.SYMMETRY_EPS)
        check_reduced_correctly(slab_copy, bulkslab.ucell, info, True)

    @with_bulk_repeat
    def test_minimal_bulk_ab(self, args):
        """Test surface_slab method ensure_minimal_bulk_ab_cell().

        This is a test for the surface slab method only. The bulk slab method
        get_minimal_ab_cell is tested in TestBulkSlab.test_minimal_ab_cell.

        Parameters
        ----------
        args : tuple
            Items are Slab, Rparams and TestInfo objects.

        Returns
        -------
        None.
        """
        slab, rpars, info = args
        slab.make_bulk_slab(rpars)
        min_ab_cell = info.bulk_properties.bulk_ucell[:2, :2]
        large_bulk = slab.make_supercell(np.diag((2, 2))).make_bulk_slab(rpars)
        slab.bulkslab = large_bulk
        rpars.superlattice_defined = False  # Avoid consistency checks
        slab.ensure_minimal_bulk_ab_cell(rpars)
        # Notice that it is enough to check that the current
        # cell is a rotated version of the minimal one
        transform = slab.bulkslab.ab_cell.dot(np.linalg.inv(min_ab_cell))
        assert np.linalg.det(transform) == pytest.approx(1.0)

    def test_inconsistent_superlattice(self, ag100):
        """Check complaints if an invalid superlattice is given."""
        slab, rpars, *_ = ag100
        slab.make_bulk_slab(rpars)
        large_bulk = slab.make_supercell(np.diag((2, 2))).make_bulk_slab(rpars)
        slab.bulkslab = large_bulk
        rpars.SUPERLATTICE = np.arange(4).reshape(2, 2)
        rpars.superlattice_defined = True
        with pytest.raises(InconsistentParameterError):
            slab.ensure_minimal_bulk_ab_cell(rpars)


class TestContains:
    """Collection of tests for the __contains__ method."""

    def test_contains_atom(self, ag100, subtests):
        """Check that an Atom is/isn't part of a slab."""
        slab, *_ = ag100
        with subtests.test('Contains first atom'):
            assert slab.atlist[0] in slab
        with subtests.test('Does not contain another atom'):
            assert Atom('C', (0, 0, 0), 1, slab) not in slab

    def test_contains_element(self, ag100, subtests):
        """Check that a chemical element is/isn't part of a slab."""
        slab, *_ = ag100
        with subtests.test('Contains Ag'):
            assert 'Ag' in slab
        with subtests.test('Does not contain other element'):
            assert 'Fe' not in slab
        with subtests.test('Does not contain invalid element'):
            assert 'invalid' not in slab

    def test_contains_layer(self, ag100, subtests):
        """Check that a Layer is/isn't part of a slab."""
        slab, rpars, *_ = ag100
        slab.create_layers(rpars)
        lay = slab.layers[0]
        with subtests.test('Contains first sublayer'):
            assert lay in slab
        with subtests.test('Does not contain another SubLayer'):
            assert lay.__class__(slab, 0) not in slab
        with subtests.test('Does not contain sublayer after clearing'):
            slab.layers = ()
            assert lay not in slab

    def test_contains_site(self, ag100, subtests):
        """Check that a Site is/isn't part of a slab."""
        slab, *_ = ag100
        site = slab.sitelist[0]
        with subtests.test('Contains first site'):
            assert site in slab
        with subtests.test('Does not contain another site'):
            assert site.__class__('C', 'non_existing') not in slab

    def test_contains_sublayer(self, ag100, subtests):
        """Check that a Sublayer is/isn't part of a slab."""
        slab, rpars, *_ = ag100
        slab.create_sublayers(rpars.SYMMETRY_EPS.z)
        sub_lay = slab.sublayers[0]
        with subtests.test('Contains first sublayer'):
            assert sub_lay in slab
        with subtests.test('Does not contain another SubLayer'):
            assert sub_lay.__class__(slab, 0) not in slab
        with subtests.test('Does not contain sublayer after clearing'):
            slab.sublayers = ()
            assert sub_lay not in slab

    def test_does_not_contain(self, ag100):
        """Check that an object is not part of a slab."""
        slab, *_ = ag100
        assert object() not in slab


class TestCoordinates:
    """Collection of tests for Cartesian/fractional coordinates."""

    cart_from_frac = {  # update_origin, expected cartpos of 1st atom
        'no update': (False, (0.3, 0.8, -1.5)),                                 # TODO: Issue #174
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
        'large values': ((5, 6, 3), (2, 2, -2)),                                # TODO: Issue #174
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

    def test_collapse_fractional(self, manual_slab_3_atoms):
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
        atom.cartpos = np.array([0.3, 0.8, -1.5])                               # TODO: Issue #174
        slab.update_fractional_from_cartesian()
        assert atom.pos == pytest.approx([0.1, 0.2, 0.3])


class TestDuplicateAtoms:
    """Tests for checking detection and removal of duplicate atoms."""

    @parametrize(info=poscar_slabs.WITH_DUPLICATE_ATOMS)
    def test_remove_duplicates(self, info, make_poscar):
        """Check correct removal of duplicate atoms."""
        slab, rpars, *_ = make_poscar(info)
        n_atoms_before = slab.n_atoms
        with not_raises(err.SlabError):
            slab.remove_duplicate_atoms(rpars.SYMMETRY_EPS,
                                        rpars.SYMMETRY_EPS.z)
        assert n_atoms_before > slab.n_atoms
        with not_raises(err.AtomsTooCloseError):
            slab.check_atom_collisions()

    def test_check_collisions_raises_without_layers(self, ag100):
        """Check complaints if layers are not available."""
        slab, *_ = ag100
        slab.layers = ()
        with pytest.raises(err.MissingLayersError):
            slab.check_atom_collisions()

    @parametrize(info=poscar_slabs.WITH_DUPLICATE_ATOMS)
    def test_with_duplicate_atoms(self, info, make_poscar):
        """Check that POSCARs with duplicates are handled correctly."""
        slab, *_ = make_poscar(info)
        with pytest.raises(err.AtomsTooCloseError):
            slab.check_atom_collisions()

    @infoless_poscar
    def test_without_duplicates(self, args):
        """Check that POSCARs without duplicates are handled correctly."""
        slab, *_ = args
        with not_raises(err.AtomsTooCloseError):
            slab.check_atom_collisions()


class TestEquivalence:
    """Collection of tests for the is_equivalent method."""

    @infoless_poscar
    def test_equivalent_copy(self, args, shuffle_slab):
        """Check that a slab is equivalent to its deepcopy."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        shuffle_slab(slab_copy)  # Make sure order does not matter
        assert slab.is_equivalent(slab_copy, eps=1e-3)

    @infoless_poscar
    def test_equivalent_translated_2d(self, args, shuffle_slab):
        """Check equivalence by translating in plane by a unit vector."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        shuffle_slab(slab_copy)  # Make sure order does not matter
        for atom in slab_copy:
            atom.pos[:2] += np.random.randint(-10, 10, size=2)
        slab_copy.update_cartesian_from_fractional()
        assert slab.is_equivalent(slab_copy, eps=1e-3)

    @infoless_poscar
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
        """Check (in)equivalence of two slabs with different no. of layers."""
        slab, rpars, *_ = ag100
        thicker, _ = slab.with_extra_bulk_units(rpars, 1)
        assert not slab.is_equivalent(thicker)

    def test_more_atoms_not_equivalent(self, ag100):
        """Check (in)equivalence of two slabs with different no. of atoms."""
        slab, *_ = ag100
        two_by_one = np.diag((2, 1))
        assert not slab.is_equivalent(slab.make_supercell(two_by_one))


class TestExtraBulk:
    """Collection of tests for adding additional bulk units to slabs."""

    @fixture(name='check_doubled')
    def fixture_check_doubled(self, subtests):
        """Check correct outcome of doubling a bulk slab."""
        def _check(double, bulk, extra=''):
            with subtests.test(f'no. atoms{extra}'):
                assert double.n_atoms == 2*bulk.n_atoms
            with subtests.test(f'ab_cell unchanged{extra}'):
                assert double.ab_cell == pytest.approx(bulk.ab_cell)
            with subtests.test(f'c vector{extra}'):
                assert double.c_vector == pytest.approx(2*bulk.c_vector)
            with subtests.test(f'no. layers{extra}'):
                assert double.n_layers == 2*bulk.n_layers
            with subtests.test(f'no atom.num duplicates{extra}'):
                at_num_counts = Counter(at.num for at in double)
                assert all(c == 1 for c in at_num_counts.values())
            if not bulk.layers:
                return
            with subtests.test(f'all{extra} layers are bulk'):
                assert all(lay.is_bulk for lay in double)
        return _check

    @fixture(name='check_extra_bulk')
    def fixture_check_extra_bulk(self, subtests):
        """Check correct addition of bulk units."""
        # There are indeed quite a lot of arguments, but packing
        # them into a dedicated class seems a bit of an overkill
        # pylint: disable-next=too-many-arguments
        def _check(bulk_appended, new_atoms, n_cells,
                   slab, rpars, n_bulk_atoms):
            n_bulk_layers_before = len(slab.bulk_layers)
            expected_extra_atoms = (n_cells * n_bulk_atoms
                                    * abs(np.linalg.det(rpars.SUPERLATTICE)))
            with subtests.test('no. bulk layers unchanged'):
                assert len(slab.bulk_layers) == n_bulk_layers_before
            with subtests.test('no. atoms'):
                assert bulk_appended.n_atoms == slab.n_atoms + len(new_atoms)
            with subtests.test('no. added atoms'):
                assert len(new_atoms) == expected_extra_atoms
            with subtests.test('no atom.num duplicates'):
                # This would fail before PR #114
                at_num_counts = Counter(at.num for at in bulk_appended)
                assert all(c == 1 for c in at_num_counts.values())
        return _check

    @with_bulk_repeat
    def test_double_thickness_no_layers(self, args, check_doubled):
        """Check correct doubling of a bulk slab that has no layers."""
        slab, rpars, info = args
        TestBulkDetect.prepare_to_detect(slab, rpars,
                                         info.bulk_properties.bulk_like_below)
        slab.detect_bulk(rpars)
        bulk = slab.bulkslab
        bulk.layers = ()
        double = bulk.with_double_thickness()
        check_doubled(double, bulk)

    @with_bulk_repeat
    def test_double_thickness_twice(self, args, check_doubled):
        """Check repeated calls to with_extra_bulk_units work correctly."""
        slab, rpars, info = args
        TestBulkDetect.prepare_to_detect(slab, rpars,
                                         info.bulk_properties.bulk_like_below)
        slab.detect_bulk(rpars)
        bulk = slab.bulkslab
        double = bulk.with_double_thickness()
        quadruple = double.with_double_thickness()
        check_doubled(double, bulk, extra=' double')
        check_doubled(quadruple, double, extra=' quadruple')

    @parametrize(n_cells=(1, 2, 3, 5))
    @with_bulk_repeat
    def test_extra_bulk(self, n_cells, args, check_extra_bulk):
        """Test function for appending bulk units to a slab."""
        slab, rpars, info = args
        TestBulkDetect.prepare_to_detect(slab, rpars,
                                         info.bulk_properties.bulk_like_below)
        slab.detect_bulk(rpars)
        bulk_appended, new_atoms = slab.with_extra_bulk_units(rpars, n_cells)
        check_extra_bulk(bulk_appended, new_atoms, n_cells, slab, rpars,
                         info.bulk_properties.n_bulk_atoms)

    @parametrize(n_cells=(-2, -1.5, 0))
    def test_extra_bulk_raises_n_cells(self, ag100, n_cells):
        """Ensure complaints for invalid n_cells."""
        slab, rpars, *_ = ag100
        with pytest.raises(ValueError):
            slab.with_extra_bulk_units(rpars, n_cells)

    def test_extra_bulk_raises_no_layers(self, ag100, subtests):
        """Check correct repeated addition of bulk units."""
        slab, rpars, *_ = ag100
        with subtests.test('no bulk layers'):
            for layer in slab.layers:
                layer.is_bulk = False
            with pytest.raises(err.MissingLayersError):
                slab.with_extra_bulk_units(rpars, 1)
        with subtests.test('no layers'):
            slab.layers = ()
            with pytest.raises(err.MissingLayersError):
                slab.with_extra_bulk_units(rpars, 1)

    def test_extra_bulk_raises_no_repeat(self, ag100):
        """Check complaints if the bulk repeat has no z."""
        slab, rpars, *_ = ag100
        rpars.BULK_REPEAT = rpars.BULK_REPEAT.copy()
        rpars.BULK_REPEAT[2] = 0
        with pytest.raises(err.SlabError):
            slab.with_extra_bulk_units(rpars, 1)

    @parametrize(n_cells=(1, 2, 3, 5))
    @with_bulk_repeat
    def test_extra_bulk_thick_bulk(self, args, n_cells, with_one_thick_bulk):
        """Check correct addition of a thicker-than-repeat bulk."""
        with_one_thick_bulk(*args)
        slab, rpars, info = args
        rpars.BULK_REPEAT = info.bulk_properties.bulk_repeat
        with pytest.raises(err.SlabError) as exc:
            slab.with_extra_bulk_units(rpars, n_cells)
        assert exc.match('negative')

    def test_extra_bulk_twice(self, ag100, subtests):
        """Check correct repeated addition of bulk units."""
        slab, rpars, *_ = ag100
        once, added_once = slab.with_extra_bulk_units(rpars, 1)
        twice, added_twice = once.with_extra_bulk_units(rpars, 1)
        _, added_two = slab.with_extra_bulk_units(rpars, 2)
        n_once = len(added_once)
        n_twice = len(added_twice)
        n_two_at_once = len(added_two)
        # All of these subtests failed before PR #114
        with subtests.test('same no. atoms each time'):
            assert n_once == n_twice
        with subtests.test('same no. atoms twice and two at once'):
            assert n_two_at_once == n_once + n_twice
        with subtests.test('no duplicates'):
            n_before_removal = twice.n_atoms
            twice.remove_duplicate_atoms(rpars.SYMMETRY_EPS,
                                         rpars.SYMMETRY_EPS.z)
            assert twice.n_atoms == n_before_removal

    @parametrize(delta_c=(-0.09, 0.2))
    def test_issue_187(self, ag100, delta_c):
        """Ensure no complaints when adding bulk to a relaxed slab (#187)."""
        slab, rpars, *_ = ag100
        # Move all the non.bulk atoms down, as it would be
        # the case when non.bulk layers are allowed to relax
        non_bulk_atoms = (at for lay in slab.non_bulk_layers
                          for at in lay.atlist)
        for atom in non_bulk_atoms:
            atom.pos[2] += delta_c
        slab.update_cartesian_from_fractional(update_origin=True)
        with not_raises(err.SlabError):
            slab.with_extra_bulk_units(rpars, 1)


@if_ase
class TestFromAse:
    """Collection of tests for the from_ase class method."""

    def test_no_ase(self, monkeypatch):
        """Check complaints when no ase module is present."""
        monkeypatch.setattr(surface_slab, '_HAS_ASE', False)
        with pytest.raises(ModuleNotFoundError):
            Slab.from_ase(None)

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

    @parametrize('slab,exc', _invalid.values(), ids=_invalid)
    def test_invalid(self, slab, exc):
        """Check complaints for invalid conditions when making bulk."""
        with pytest.raises(exc):
            slab.make_bulk_slab(Rparams())

    @with_bulk_repeat
    def test_nr_of_atoms(self, args):
        """Test expected number of atoms in bulk slab for valid POSCARs."""
        slab, rpars, info = args
        bulk_slab = slab.make_bulk_slab(rpars)
        assert bulk_slab.n_atoms == info.bulk_properties.n_bulk_atoms

    @with_bulk_repeat
    def test_top_and_bottom_atoms(self, args, subtests):
        """Check that the top and bottom bulk atoms remain unchanged."""
        slab, rpars, *_ = args
        top_bulk_atom = max(slab.bulk_atoms, key=lambda at: at.pos[2])
        bot_bulk_atom = min(slab.bulk_atoms, key=lambda at: at.pos[2])
        slab.create_sublayers(rpars.SYMMETRY_EPS.z)
        top_bulk_sublayer = next(lay for lay in slab.sublayers
                                 if top_bulk_atom in lay)
        bot_bulk_sublayer = next(lay for lay in slab.sublayers
                                 if bot_bulk_atom in lay)
        bulk = slab.make_bulk_slab(rpars)
        top, bot = bulk.top_atom, bulk.bottom_atom
        with subtests.test('c centering'):
            assert top.pos[2] + bot.pos[2] == pytest.approx(1)
        with subtests.test('top atom unchanged'):
            assert top.num in {at.num for at in top_bulk_sublayer}
        with subtests.test('bottom atom unchanged'):
            assert bot.num in {at.num for at in bot_bulk_sublayer}

    @with_bulk_repeat
    def test_ucell(self, args):
        """Test expected number of atoms in bulk slab for valid POSCARs."""
        slab, rpars, info = args
        bulk_slab = slab.make_bulk_slab(rpars)
        atol = float(0.2*rpars.SYMMETRY_EPS)
        assert np.allclose(bulk_slab.ucell, info.bulk_properties.bulk_ucell,
                           atol=atol)

    def test_warning_a_larger_b(self, ag100, caplog):
        """Test expected number of atoms in bulk slab for valid POSCARs."""
        slab, rpars, *_ = ag100
        rpars.superlattice_defined = True  # Needed for check to happen
        slab.ab_cell.T[0] *= 2
        slab.update_cartesian_from_fractional()
        slab.make_bulk_slab(rpars)
        assert 'does not follow standard convention' in caplog.text


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

    @infoless_poscar
    def test_one_operation(self, args, check_identical):
        """Check correct result of reverting one unit-cell operation."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        slab.rotate_unit_cell(6)
        slab.revert_unit_cell()
        check_identical(slab, slab_copy)

    @fixture(name='with_few_operations')
    def fixture_with_few_operations(self):
        """Apply a bunch of unit-cell operations to slab."""
        rot_orders_for_cell = {
            'oblique': (),
            'rectangular': (2,),
            'hexagonal': (2, 3, 6),
            'rhombic': (2,),
            'square': (2, 4)
            }
        def _rotate_once(slab, orders):
            order = next(orders, None)
            if not order:
                return
            slab.rotate_unit_cell(order)

        def _apply(slab):
            cell_shape, _ = leedbase.checkLattice(slab.ab_cell.T)
            rot_orders = iter(rot_orders_for_cell[cell_shape])
            _rotate_once(slab, rot_orders)
            slab.transform_unit_cell_2d(((1, 1), (1, -1)))
            _rotate_once(slab, rot_orders)
            slab.translate_atoms_2d((0.25, 0.32))
            _rotate_once(slab, rot_orders)
            slab.translate_atoms_c(0.29)
            slab.collapse_cartesian_coordinates()
        return _apply

    @infoless_poscar
    def test_few_operations(self, args, with_few_operations, check_identical):
        """Same as above, but reverting a few operations."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        with_few_operations(slab)
        slab.revert_unit_cell()
        check_identical(slab, slab_copy)

    @infoless_poscar
    def test_nothing_to_undo_collapses(self, args, check_identical):
        """Check that reverting with no operations collapses coordinates."""
        slab, *_ = args
        # Translate some atoms out of the unit cell
        for atom in slab:
            atom.pos = atom.pos + (0.52, -0.58, 1.7)
        slab.update_cartesian_from_fractional(update_origin=True)
        slab_copy = deepcopy(slab)
        slab_copy.collapse_cartesian_coordinates()
        slab.revert_unit_cell()
        check_identical(slab, slab_copy)

    @infoless_poscar
    def test_undo_none_collapses(self, args, with_few_operations,
                                 check_identical):
        """Check that reverting no operations always collapses coordinates."""
        slab, *_ = args
        # Translate some atoms out of the unit cell
        for atom in slab:
            atom.pos = atom.pos + (0.52, -0.58, 1.7)
        slab.update_cartesian_from_fractional(update_origin=True)
        with_few_operations(slab)
        slab_copy = deepcopy(slab)
        slab_copy.collapse_cartesian_coordinates()
        slab.revert_unit_cell(slab.ucell_mod)
        check_identical(slab, slab_copy)

    @parametrize(n_undo=(1, 3, 4))
    @infoless_poscar
    def test_revert_in_two_steps(self, args, n_undo, check_identical,
                                 with_few_operations):
        """Check that reverting all operations in two steps is equivalent."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        with_few_operations(slab)
        modified = deepcopy(slab)

        # Undo a few operations, then the rest...
        slab.revert_unit_cell(restore_to=slab.ucell_mod[:-n_undo])
        slab.revert_unit_cell()
        # ...and also all of them at once
        modified.revert_unit_cell()

        check_identical(slab, slab_copy)
        check_identical(slab, modified)

    def test_raises(self, ag100):
        """Check complaints for invalid operation types."""
        slab, *_ = ag100
        slab.ucell_mod.append(('invalid', 'invalid'))
        with pytest.raises(RuntimeError):
            slab.revert_unit_cell()


class TestSlabLayers:
    """Collection of tests concerning slab (sub)layers."""

    @fixture(name='check_layers_correct')
    def fixture_check_layers_correct(self, subtests):
        """Check correct identification of layers from known info."""
        def _check(slab, cuts, info):
            lay_info = info.layer_properties
            with subtests.test('no. layers'):
                assert slab.n_layers == pytest.approx(lay_info.n_layers)
            with subtests.test('cut positions'):
                assert cuts == pytest.approx(lay_info.cuts)
            with subtests.test('atoms per layer'):
                n_atoms = [lay.n_atoms for lay in slab.layers]
                assert n_atoms == pytest.approx(lay_info.n_atoms_per_layer)
        return _check

    def test_auto_cut_below_bulk_cut(self, ag100, check_layers_correct):
        """Check correct layer identification with a useless auto-cut."""
        slab, rpars, info = ag100
        bulk_cuts = info.bulk_properties.bulk_cuts
        rpars.LAYER_CUTS = LayerCuts.from_string(
            f'dz(0.2) < {max(bulk_cuts) - .02} < dz(1.2)'
            )
        cuts = slab.create_layers(rpars, bulk_cuts=bulk_cuts)
        check_layers_correct(slab, cuts, info)

    @with_layers
    def test_bulk_layers(self, args, subtests):
        """Test Slab.bulk_layers property."""
        slab, rpars, info = args
        slab.create_layers(rpars)
        with subtests.test('no. bulk layers'):
            assert len(slab.bulk_layers) == info.layer_properties.n_bulk_layers
        with subtests.test('no. atoms in bulk layers'):
            n_bulk_atoms = sum(lay.n_atoms for lay in slab.bulk_layers)
            assert n_bulk_atoms == info.bulk_properties.n_bulk_atoms

    def test_bulkslab_layers_always_bulk(self, with_one_thick_bulk):
        """Check that all layers of a BulkSlab are bulk."""
        args = CasePOSCARSlabs().case_poscar_fe3o4_001_cod()
        with_one_thick_bulk(*args)
        slab, rpars, *_ = args
        bulk = slab.bulkslab
        rpars.LAYER_CUTS = LayerCuts.from_string('dz(1.0)')
        bulk.create_layers(rpars)  # 12 layers
        assert bulk.n_layers == 12
        assert all(lay.is_bulk for lay in bulk.layers)

    @with_layers
    def test_create_layers(self, args, check_layers_correct):
        """Check that layers are created correctly."""
        slab, rpars, info = args
        cuts = slab.create_layers(rpars)
        check_layers_correct(slab, cuts, info)

    @with_layers
    def test_create_sublayers(self, args):
        """Check that sublayers are created correctly."""
        slab, rpars, info = args
        slab.create_sublayers(rpars.SYMMETRY_EPS.z)
        assert slab.n_sublayers == info.layer_properties.n_sublayers

    def test_empty_layer_warning(self, ag100, caplog):
        """Check that layers are created correctly."""
        slab, rpars, *_ = ag100
        rpars.LAYER_CUTS = LayerCuts.from_string('dz(1.2) < 0.72 0.78')
        rpars.N_BULK_LAYERS = 1
        slab.create_layers(rpars)
        # check that we log a warning about the empty layer
        assert 'will be deleted' in caplog.text

    @with_layers
    def test_full_update_with_layers_defined(self, args):
        """Test full_update method with layers already defined."""
        slab, rpars, *_ = args
        n_layers_orignal = slab.n_layers
        n_atoms_per_layer_original = [lay.n_atoms for lay in slab.layers]
        # Shift some atoms out of the unit cell
        for atom in slab:
            atom.pos += 0.5
        slab.full_update(rpars)
        n_atoms_per_layer_after = [lay.n_atoms for lay in slab.layers]
        assert slab.n_layers == n_layers_orignal
        assert n_atoms_per_layer_after == n_atoms_per_layer_original
        # Check that all atoms are in the unit cell
        eps = 1e-8
        assert np.all(-eps <= at.pos <= 1 + eps for at in slab)

    @with_layers
    def test_full_update_without_topat_ori_z(self, args):
        """Test full_update when topat_ori_z is not available."""
        slab, rpars, info = args
        slab.topat_ori_z = None
        slab.full_update(rpars)
        assert slab.topat_ori_z is not None
        assert slab.n_layers == info.layer_properties.n_layers
        slab.create_sublayers(rpars.SYMMETRY_EPS.z)
        assert slab.n_sublayers == info.layer_properties.n_sublayers

    def test_from_bulk_slab_clears_layers(self, ag100):
        """Check that SurfaceSlab.from_slab(bulk) has no (sub)layers."""
        slab, rpars, *_ = ag100
        bulk = slab.make_bulk_slab(rpars)
        assert bulk.n_layers == 1
        new_slab = Slab.from_slab(bulk)
        assert not new_slab.layers
        assert not new_slab.sublayers

    def test_interlayer_gap_positive(self, ag100):
        """Ensure that gaps between layers are always positive."""
        slab, *_ = ag100
        assert all(d > 0 for d in slab.interlayer_gaps)

    @with_layers
    def test_interlayer_spacing_raises_without_layers(self, args):
        """Test that interlayer spacing raises without layers defined."""
        slab, *_ = args
        slab.layers = ()
        with pytest.raises(err.MissingLayersError):
            _ = slab.smallest_interlayer_gap

    @with_layers
    def test_n_atoms_per_sublayer(self, args):
        """Test the expected number of atoms per sublayer."""
        slab, rpars, info = args
        slab.create_sublayers(rpars.SYMMETRY_EPS.z)
        n_atoms_sublayers = [sublay.n_atoms for sublay in slab.sublayers]
        assert n_atoms_sublayers == info.layer_properties.n_atoms_per_sublayer

    @with_layers
    def test_smallest_interlayer_gap(self, args):
        """Test that interlayer spacing is correctly calculated."""
        slab, _, info = args
        expected = info.layer_properties.smallest_interlayer_gap
        spacing = slab.smallest_interlayer_gap
        assert spacing == pytest.approx(expected, abs=1e-5)

    def test_smallest_interlayer_gap_few_layers(self, manual_slab_3_atoms):
        """Check complaints when there's not enough layers."""
        slab = manual_slab_3_atoms
        with pytest.raises(err.TooFewLayersError):
            _ = slab.smallest_interlayer_gap

    @infoless_poscar
    def test_sublayer_sorting(self, args, subtests):
        """Check that sublayers are created with the expected sort order."""
        slab, rpars, *_ = args
        slab.create_sublayers(rpars.SYMMETRY_EPS.z)

        # Expect: (i) top to bottom, (ii) alphabetically
        # by element at same z (within eps)
        # Group by same z (within eps)
        z_above = slab.sublayers[0].cartbotz
        same_z = []
        by_z = [same_z]
        for layer in slab.sublayers:
            this_z = layer.cartbotz
            if abs(z_above - this_z) < rpars.SYMMETRY_EPS.z:
                same_z.append(layer)
                z_above = min(z_above, this_z)                                  # TODO: max when flipping .cartpos[2] -- Issue #174
            else:
                same_z = [layer]
                by_z.append(same_z)
                z_above = this_z
        # Compute min and max z positions
        z_min_max = []
        for same_z_group in by_z:
            pos = [lay.cartbotz for lay in same_z_group]
            z_min_max.append((min(pos), max(pos)))                              # TODO: redo when flipping .cartpos[2] -- Issue #174
        with subtests.test('z sorting'):
            assert all(bottom_of_above <= top_of_below                          # TODO: redo when flipping .cartpos[2] -- Issue #174
                       for (_, bottom_of_above), (top_of_below, _)
                       in pairwise(z_min_max))
        for same_z_group in by_z:
            z_approx = same_z_group[0].pos[2]
            elements = [lay.element for lay in same_z_group]
            with subtests.test(f'element-sort at pos~{z_approx:.4f}'):
                assert sorted(elements) == elements

    def test_sublayer_sorting_small_dz(self):
        """Check that small z displacements of sublayers maintain sorting."""
        slab, rpars, *_ = CasePOSCARSlabs().case_poscar_tio2_small()
        eps = rpars.SYMMETRY_EPS.z
        slab.create_sublayers(eps)
        sorting_before = [lay.element for lay in slab.sublayers]

        # Layers 15 and 19 are translation-equivalent (both with 2 O)           # TODO: Issue #174
        for atom in slab.sublayers[15]:
            atom.cartpos[2] += 0.05*eps
        for atom in slab.sublayers[19]:
            atom.cartpos[2] -= 0.05*eps
        slab.update_fractional_from_cartesian()
        slab.create_sublayers(eps)
        sorting_after = [lay.element for lay in slab.sublayers]
        assert sorting_after == sorting_before


class TestSlabRaises:
    """Collection of tests for diverse exception-raising conditions."""

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

    _props = {
        'ab_cell': err.InvalidUnitCellError,
        'c_vector': err.InvalidUnitCellError,
        'fewest_atoms_sublayer': err.MissingSublayersError,
        'interlayer_gaps': err.MissingLayersError,
        'smallest_interlayer_gap': err.MissingLayersError,
        }

    @parametrize('attr_name,exc', _props.items(), ids=_props)
    def test_property_raises(self, attr_name, exc):
        """Check complaints when accessing a property of an empty slab."""
        with pytest.raises(exc):
            _ = getattr(Slab(), attr_name)


@todo
class TestSlabSites:                                                            # TODO: probably better before reworking initSites
    """Collection of tests for generation of sites."""


class TestSorting:
    """Collection of tests for slab sorting."""

    @infoless_poscar
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

    @infoless_poscar
    def test_sort_original(self, args, shuffle_slab):
        """Check correct sorting of a Slab to original atom numbers."""
        slab, *_ = args
        original_atlist = deepcopy(slab.atlist)
        shuffle_slab(slab)
        slab.sort_original()
        assert all(all(at1.pos == at2.pos) and at1.el == at2.el
                   for at1, at2 in zip(slab, original_atlist))

    @infoless_poscar
    def test_simple_sort_by_z(self, args, shuffle_slab):
        """Check correct outcome of sorting by z coordinates."""
        slab, *_ = args
        shuffle_slab(slab)
        slab.sort_by_z()
        assert all(at1.pos[2] <= at2.pos[2] for at1, at2 in pairwise(slab))

    @infoless_poscar
    @parametrize(bottom_to_top=(True, False))
    def test_z_sort(self, args, bottom_to_top):
        """Check successful sorting of atoms by out-of-plane position."""
        slab, *_ = args
        slab.sort_by_z(bottom_to_top=bottom_to_top)
        _ordered = operator.ge if bottom_to_top else operator.le                # TODO: swap when flipping .cartpos -- Issue #174
        assert all(_ordered(at1.cartpos[2], at2.cartpos[2])
                   for at1, at2 in pairwise(slab))


class TestSuperAndSubCell:
    """Collection of tests for creation of larger and smaller slab versions."""

    invalid = {
        'not integer': (np.eye(2)*0.5, NonIntegerMatrixError),
        'singular': (np.eye(2)*0, SingularMatrixError),
        'non-2x2': (np.eye(3), ValueError),
        }
    valid = {
        'identity': np.diag((1, 1)),
        '2x2': np.diag((2, 2)),
        '2x1': np.diag((2, 1)),
        'off-diagonal': np.array([[2, 1], [2, 2]]),
        }

    @infoless_poscar
    def test_make_supercell_identity_does_nothing(self, args, check_identical):
        """Check that identity matrix does not change the slab."""
        slab, *_ = args
        supercell = slab.make_supercell(np.identity(2))
        check_identical(slab, supercell)

    @fixture(name='with_supercell', scope='class')
    @parametrize('transform', valid.values(), ids=valid)
    @parametrize_with_cases('args', cases=CasePOSCARSlabs.case_infoless_poscar,
                            scope='class')
    def fixture_with_supercell(self, transform, args):
        """Return a supercell of a slab.

        Parameters
        ----------
        transform : numpy.ndarray
            (Integer) transformation matrix to create a supercell.
        args : tuple
            Items are the slab whose supercell is created, an Rparams
            object to be used, and a TestInfo.

        Returns
        -------
        supercell : SurfaceSlab
            The supercell version of slab.
        transform : numpy.ndarray
            The transformation matrix used.
        slab : SurfaceSlab
            The slab whose supercell is created
        rpars : Rparams
            Run parameters object for slab
        info : TestInfo
            Information for tests relative to slab.
        """
        slab, rpars, info = args
        supercell = slab.make_supercell(transform)
        return supercell, transform, slab, rpars, info

    def test_supercell_valid(self, with_supercell, subtests):
        """Check that supercell is created correctly."""
        supercell, transform, slab, rpars, *_ = with_supercell
        # The original implementation on master from commit
        # f3512ce30cf93f78e12ffc63224a4f285625269e caused
        # duplicate atoms upon creating supercells with
        # some off-diagonal supercell matrices
        supercell.remove_duplicate_atoms(rpars.SYMMETRY_EPS,
                                         rpars.SYMMETRY_EPS.z)
        n_cells = abs(np.linalg.det(transform))
        with subtests.test('ab cell'):
            assert np.allclose(supercell.ab_cell,
                               slab.ab_cell.dot(transform.T))
        with subtests.test('no. atoms'):
            assert supercell.n_atoms == slab.n_atoms * n_cells

    @parametrize('transform', valid.values(), ids=valid)
    @infoless_poscar
    def test_supercell_does_not_change_original_slab(self, transform, args,
                                                     check_identical):
        """Check that supercell creation does not affect input slab."""
        slab, *_ = args
        slab_copy = deepcopy(slab)
        _ = slab.make_supercell(transform)
        check_identical(slab, slab_copy)

    def test_subcell_reverses_supercell(self, with_supercell, check_identical):
        """Test that make_subcell reverses make_supercell."""
        supercell, transform, slab, rpars, *_ = with_supercell
        subcell = supercell.make_subcell(rpars, transform)
        slab.sort_original()
        subcell.sort_original()
        check_identical(slab, subcell)

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

    _z_changing_matrix = {
        'reflect': np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]]),
        'mixing a, b, c': np.array([[1/np.sqrt(2), 0, 1/np.sqrt(2)],
                                    [0, 1, 0],
                                    [-1/np.sqrt(2), 0, 1/np.sqrt(2)]]),
        }

    @parametrize('matrix', _z_changing_matrix.values(), ids=_z_changing_matrix)
    def test_matrix_transform_changes_z(self, matrix, ag100):
        """Check outcome of a matrix transformation that mixes a, b, and c."""
        slab, rpars, *_ = ag100
        original_ucell = slab.ucell.copy()
        slab.make_bulk_slab(rpars)
        assert slab.bulkslab is not None
        slab.apply_matrix_transformation(matrix)
        assert np.allclose(slab.ucell, matrix.dot(original_ucell))
        # check that bulkslab is reset to None in this case
        assert slab.bulkslab is None

    _non_z_changing_matrix = {
        'mixing': np.array([[1/np.sqrt(2), -1/np.sqrt(2), 0],
                            [1/np.sqrt(2), 1/np.sqrt(2), 0],
                            [0, 0, 1]]),
        }

    @parametrize('matrix', _non_z_changing_matrix.values(),
                 ids=_non_z_changing_matrix)
    def test_matrix_transform_propagated_to_bulk(self, matrix, ag100):
        """Check that applying a matrix transform is reflected on bulkslab."""
        slab, rpars, _ = ag100
        slab.make_bulk_slab(rpars)
        original_ucell = slab.ucell.copy()
        orignal_bulk_ucell = slab.bulkslab.ucell.copy()
        slab.apply_matrix_transformation(matrix)
        assert slab.bulkslab is not None
        assert np.allclose(slab.bulkslab.ucell, matrix.dot(orignal_bulk_ucell))
        assert np.allclose(slab.ucell, matrix.dot(original_ucell))

    def test_project_c_to_z(self, subtests):
        """Check that c is along z after projection."""
        slab, *_ = CasePOSCARSlabs().case_poscar_36carbon_atoms_cm()
        slab_copy = deepcopy(slab)
        slab.project_c_to_z()
        with subtests.test('c along z'):
            assert np.allclose(slab.c_vector[:2], 0)

        # Atoms 4 & 5 are wrapped around
        atom_pairs = zip(slab.atlist[~4:5], slab_copy.atlist[~4:5])
        with subtests.test('atom positions'):
            assert all(np.allclose(at.cartpos, at_copy.cartpos)
                       for at, at_copy in atom_pairs)

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


class TestUnitCellReduction:
    """Tests for minimization of various bits of the unit cell."""

    @parametrize_with_cases('args', cases=poscar_slabs,
                            has_tag=Tag.NON_MINIMAL_CELL)
    def test_ab_cell_minimization(self, args):
        """Check correct identification of minimal ab cell (by area)."""
        slab, rpars, info = args
        ab_slab = slab.ab_cell.T.copy()
        ab_small = slab.get_minimal_ab_cell(rpars.SYMMETRY_EPS,
                                            rpars.SYMMETRY_EPS.z)
        area_ratio = np.linalg.det(ab_slab)/np.linalg.det(ab_small)
        assert area_ratio == pytest.approx(info.poscar.n_cells)

    _minimal = {'cases': poscar_slabs,
                'filter': exclude_tags(Tag.NON_MINIMAL_CELL,
                                       Tag.NO_INFO)}

    @parametrize_with_cases('args', **_minimal)
    def test_ab_cell_already_minimal(self, args):
        """Check correct identification of slabs with already minimal ab."""
        slab, rpars, *_ = args
        with pytest.raises(err.AlreadyMinimalError):
            slab.get_minimal_ab_cell(rpars.SYMMETRY_EPS, rpars.SYMMETRY_EPS.z)


class TestVacuumGap:
    """Collection of tests for checking the extent and position of vacuum."""

    @staticmethod
    def _check_raises(slab, exc):
        """Ensure that slab.check_vacuum_gap raises exc."""
        with pytest.raises(exc):
            slab.check_vacuum_gap()

    @parametrize_with_cases(
        'args', cases=CasePOSCARSlabs,
        filter=exclude_tags(Tag.VACUUM_GAP_MIDDLE,
                            Tag.VACUUM_GAP_SMALL,
                            Tag.VACUUM_GAP_ZERO,
                            Tag.NO_INFO)
        )
    def test_vacuum_ok(self, args, caplog):
        """Check no complaints if there is no problem with vacuum."""
        slab, *_ = args
        has_atoms_at_edges = any(at.pos[2] <= 1e-4  # slab is collapsed
                                 for at in slab)
        with not_raises(err.VacuumError), caplog.at_level(logging.INFO):
            slab.check_vacuum_gap()
        if has_atoms_at_edges:
            assert 'Atoms were shifted' in caplog.text

    @parametrize_with_cases('args', cases=poscar_slabs,
                            has_tag=Tag.VACUUM_GAP_MIDDLE)
    def test_vacuum_gap_mid_slab(self, args):
        """Check complaints when the vacuum is in the middle of a slab."""
        slab, *_ = args
        self._check_raises(slab, err.WrongVacuumPositionError)

    @if_ase
    def test_no_vacuum_gap(self):
        """Check complaints when atoms are very close to both c=0 and c=1."""
        slab = Slab.from_ase(
            cases_ase.ase.build.fcc110('Ni', size=(1, 1, 6), vacuum=0.0005)
            )
        with pytest.raises(err.NoVacuumError):
            slab.check_vacuum_gap()

    @parametrize_with_cases('args', cases=poscar_slabs,
                            has_tag=Tag.VACUUM_GAP_SMALL)
    def test_vacuum_gap_small(self, args):
        """Check complaints with a small vacuum gap."""
        slab, *_ = args
        self._check_raises(slab, err.NotEnoughVacuumError)


@parametrize_with_cases('args', cases=CasePOSCARSlabs,
                        has_tag=Tag.SURFACE_ATOMS)
def test_get_surface_atoms(args):
    """Check correct identification of atoms visible from vacuum."""
    slab, rpars, info = args
    surface_atoms = slab.getSurfaceAtoms(rpars)
    surface_atoms_nums = {at.num for at in surface_atoms}
    assert surface_atoms_nums == set(info.surface_atoms.surface_atom_nums)


@parametrize_with_cases('args', cases=CasePOSCARSlabs,
                        has_tag=Tag.NEAREST_NEIGHBOURS)
def test_nearest_neighbors(args):
    """Test function get_nearest_neighbors."""
    slab, _, info = args
    nearest_neighbors = slab.get_nearest_neighbours()
    expected = info.nearest_neighbors.nearest_neighbor_distances
    assert len(nearest_neighbors) == len(expected)
    for atom in slab:
        assert nearest_neighbors[atom] == pytest.approx(expected[atom.num])


def test_translation_symmetry_different_species():
    """Check that an MgO slab is not equivalent upon Mg->O translation."""
    slab, rpars, *_ = CasePOSCARSlabs().case_poscar_mgo()
    mg_to_oxygen = slab.ab_cell.T[0] / 2
    assert not slab.is_translation_symmetric(mg_to_oxygen, rpars.SYMMETRY_EPS)


def test_ucell_ori_after_clear_symmetry_and_ucell_history(ag100):
    """Ensure modifications to the unit cell don't go into ucell_ori."""
    slab, *_ = ag100
    assert slab.celltype == 'unknown'
    ucell_ori = slab.ucell.copy()
    slab.clear_symmetry_and_ucell_history()
    slab.ucell *= 3
    assert np.allclose(ucell_ori, slab.ucell_ori)
