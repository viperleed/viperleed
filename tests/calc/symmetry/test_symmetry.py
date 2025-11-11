"""Tests for module viperleed.calc.symmetry.

Contains tests for symmetry-detection routines.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-03-26'
__license__ = 'GPLv3+'

from contextlib import contextmanager
import copy
import logging

import numpy as np
import pytest
from pytest import approx

from pytest_cases import fixture, parametrize, parametrize_with_cases
from pytest_cases.filters import id_has_suffix

from viperleed.calc import symmetry
from viperleed.calc.lib.math_utils import angle as angle_radians

from ...helpers import duplicate_all
from ..poscar_slabs import make_poscar_ids
from ..tags import CaseTag
from . import simple_slabs
from .conftest import get_cases


def angle(vec1, vec2):
    """Return the (rounded) angle between vec1 and vec2 in degrees."""
    return round(np.degrees(angle_radians(vec1, vec2)))


def hermann(group):
    """Return the "Hermann-Maugin" part of group."""
    group = str(group)
    return group.split('[', maxsplit=1)[0]


def _reconstruct_case_id(case):
    """Return a full pytest-cases-style id for case."""
    base_id = case.id
    test_info = case.params.get('info', None)
    if test_info and 'poscar' in base_id:
        # Add POSCAR name
        return base_id + '-' + make_poscar_ids()(test_info)

    # See if it's a double_bulk
    lazy_bulk = case.params.get('bulk', None)
    if 'double_bulk' in base_id and lazy_bulk:
        return base_id + '-' + make_poscar_ids('thick_bulk')(lazy_bulk)

    # See if it's a CaseSimpleSlabs
    if case.func.__name__ in simple_slabs.CaseSimpleSlabs.__dict__:
        extra = (f'{k}={v}' for k, v in case.params.items())
        return '-'.join((base_id, *extra))

    return base_id


@contextmanager
def may_fail(case, known_failures, strict=False):
    """XFAIL if case is known to be problematic."""
    reason = known_failures.get(case.id, '')
    if not reason:
        reason = known_failures.get(_reconstruct_case_id(case), '')
    try:
        yield None
    except AssertionError:
        if reason:
            pytest.xfail(reason)
        raise
    # No failure if we get here
    if reason and strict:
        raise AssertionError(
            f'XPASS: This test should FAIL with reason={reason!r}'
            )


class TestPlaneGroupFinding:
    """Collection of tests for finding plane groups."""

    def test_any_plane_group_found(self, with_plane_group):
        """Ensure that symmetry detection is successful."""
        slab, *_ = with_plane_group()
        assert slab.foundplanegroup != 'unknown'

    known_incorrect_groups = {
        'hex_cmm_10': 'Known incorrect plane group p2',
        'hex_cmm_01': 'Known incorrect plane group p2',
        'poscar_diamond': 'Known incorrect plane group pm instead of rcm',
        }

    def test_correct_plane_group(self, with_plane_group, first_case):
        """Check the correct identification of the plane group."""
        slab, *_, info = with_plane_group()
        if not info.symmetry.hermann:
            pytest.skip('No symmetry information available')
        with may_fail(first_case, self.known_incorrect_groups, strict=True):
            assert hermann(slab.planegroup) == info.symmetry.hermann

    _known_incorrect_rotations = {
        'hex_cmm_10': 'BUG: no rotation applied',
        'hex_cmm_01': 'BUG: no rotation applied',
        'hex_cm_10': 'BUG: no rotation applied',
        'hex_cm_01': 'BUG: no rotation applied',
        'hex_cm_21': 'BUG: no rotation applied',
        'hex_cm_12': 'BUG: no rotation applied',
        }

    @parametrize_with_cases('args',
                            cases=get_cases('all'),
                            has_tag=CaseTag.NEED_ROTATION)
    # pylint: disable-next=too-many-arguments  # All fixtures
    def test_cell_rotated(self, args, caplog, re_match, subtests, first_case):
        """Check rotation of slabs that need one to get the group right."""
        slab, param, info, *_ = args
        a_before, _ = slab.ab_cell.T.copy()
        area_before = np.linalg.det(slab.ab_cell)
        symmetry.findSymmetry(slab, param)
        a_after, _ = slab.ab_cell.T
        area_after = np.linalg.det(slab.ab_cell)
        with subtests.test('successful logging'):
            with may_fail(first_case, self._known_incorrect_rotations):
                assert re_match(r'.*unit.*cell.*change.*higher.*symmetry.*',
                                caplog.text)
        with subtests.test('correct rotation'):
            if first_case.id == 'hex_cm_12':
                pytest.xfail('Inconsistent rotations between detection and '
                             'reduction of symmetry. Angles differ by 180')
            with may_fail(first_case, self._known_incorrect_rotations):
                assert angle(a_before, a_after) == info.symmetry.rotation
        with subtests.test('area unchanged'):
            assert area_after == approx(area_before)

    _known_incorrect_groups_fix = {
        'pmg-ucell=rectangular': 'Known incorrect plane group p2',
        'pmg-ucell=square': 'Known incorrect plane group p2',
        'pgg-ucell=rectangular': 'Known incorrect plane group rcmm',
        'pgg-ucell=square': 'Known incorrect plane group rcmm',
        'hex_cm_10': 'Known incorrect plane group p1',
        'hex_cm_01': 'Known incorrect plane group p1',
        'hex_cm_21': 'Known incorrect plane group p1',
        'hex_cm_12': 'Known incorrect plane group p1',
        'hex_cmm_10': 'Known incorrect plane group p2',
        'hex_cmm_01': 'Known incorrect plane group p2',
        }

    @fixture
    @parametrize_with_cases('args', cases=get_cases('simple'))
    def with_plane_group_fix_origin(self, args):
        """Find plane group, assuming Cartesian origin is high symmetry."""
        slab, param, info, *_ = args
        param.SYMMETRY_FIND_ORI = False
        symmetry.findSymmetry(slab, param, forceFindOri=False)
        return slab, param, info

    def test_correct_plane_group_fix_origin(self, with_plane_group_fix_origin,
                                            first_case):
        """Check the correct identification of the plane group."""
        slab, *_, info = with_plane_group_fix_origin
        with may_fail(first_case, self._known_incorrect_groups_fix):
            assert hermann(slab.planegroup) == info.symmetry.hermann


XYZ_DISPLACEMENTS = (
    ([0.2, 0, 0],),
    ([0, 0.2, 0],),
    ([0, 0, 0.2],)
    )


class TestSymmetryConstraints:
    """Collection of tests for symmetry-related atomic constraints."""

    @pytest.mark.xfail(reason='BUG: symmetry.enforceSymmetry L989',
                       strict=True)
    def test_raises_without_planegroup(self, slab_p1):
        """Check complaints when constraining without symmetry knowledge."""
        slab, param, *_ = slab_p1
        with pytest.raises(RuntimeError) as exc:
            symmetry.enforceSymmetry(slab, param)
        assert exc.match(r'.*without.*known plane\s+group')

    @pytest.mark.xfail(reason='BUG: UnboundLocalError toprotsym L1028',
                       strict=True)
    def test_raises_without_valid_subgroup(self, slab_p6m):
        """Check complaint when constraining before a valid subgroup is set."""
        slab, param, *_ = slab_p6m
        symmetry.findSymmetry(slab, param)
        with pytest.raises(ValueError) as exc:
            symmetry.enforceSymmetry(slab, param, 'cm[1 -1]')
        assert exc.match(r'.*inconsistent slab plane\s+group')

    @pytest.mark.xfail(reason='BUG: No complaint for invalid subgroup',
                       strict=True)
    def test_raises_with_invalid_subgroup(self, slab_p6m):
        """Check complaint when constraining before a valid subgroup is set."""
        slab, param, *_ = slab_p6m
        symmetry.findSymmetry(slab, param)
        with pytest.raises(ValueError) as exc:
            symmetry.enforceSymmetry(slab, param, 'p4m')
        assert exc.match('.*Not a subgroup')

    _known_invalid_constraints = {
        'pm_10': 'Known to sometimes fail with a random shift',
        'cm_1m1': 'Known to sometimes fail with a random shift',
        'rcm-ucell=rectangular': 'Known to sometimes fail with a random shift',
        'rcm-ucell=square': 'Known to sometimes fail with a random shift',
        'hex_cm_11': 'Known to sometimes fail with a random shift',
        'hex_cm_1m1': 'Known to regularly fail with a random shift',
        'hex_cm_10': 'Known to regularly fail with a random shift',
        'hex_cm_01': 'Known to regularly fail with a random shift',
        'hex_cm_21': 'Known to sometimes fail with a random shift',
        'hex_cm_12': 'Known to regularly fail with a random shift',
        'hex_cmm_10': 'Known incorrect plane group p2',
        'hex_cmm_01': 'Known incorrect plane group p2',
        'square_pm_10': 'Known to sometimes fail with a random shift',
        'square_cm_1m1': 'Known to sometimes fail with a random shift',
        }

    def test_correct_free_directions(self, with_symmetry_constraints,
                                     first_case):
        """Check the correct identification of atom constraints."""
        slab, _, info = with_symmetry_constraints()
        if not info.symmetry.link_groups:
            pytest.skip('Not enough symmetry information available')
        with may_fail(first_case, self._known_invalid_constraints):
            for atom in slab:
                if atom.num in info.symmetry.on_planes:
                    assert isinstance(atom.freedir, np.ndarray)
                elif atom.num in info.symmetry.on_axes:
                    assert not atom.freedir
                else:
                    assert atom.freedir == 1

    _known_invalid_linking = {
        'rcm-ucell=rectangular': 'Known to regularly fail with a random shift',
        'rcm-ucell=square': 'Known to regularly fail with a random shift',
        'rcmm-ucell=rectangular': 'Known to regularly fail with a random shift',
        'rcmm-ucell=square': 'Known to regularly fail with a random shift',
        'hex_cm_1m1': 'Known to often fail with a random shift',
        'hex_cm_11': 'Known to sometimes fail with a random shift',
        'hex_cm_10': 'Known to regularly fail with a random shift',
        'hex_cm_01': 'Known to regularly fail with a random shift',
        'hex_cm_21': 'Known to regularly fail with a random shift',
        'hex_cm_12': 'Known to regularly fail with a random shift',
        'hex_cmm_10': 'Known incorrect plane group p2',
        'hex_cmm_01': 'Known incorrect plane group p2',
        'square_cm_1m1': 'Known to sometimes fail with a random shift',
        'square_cm_11': 'Known to sometimes fail with a random shift',
        'square_pm_10': 'Known to sometimes fail with a random shift',
        'square_pm_01': 'Known to sometimes fail with a random shift',
        'square_pg_10': 'Known to often fail with a random shift',
        'pm_10': 'Known to sometimes fail with a random shift',
        'pg_10': 'Known to often fail with a random shift',
        'cm_1m1': 'Known to often fail with a random shift',
        'poscar_fe3o4_001_cod': 'Known to rarely fail with a random shift',
        'p6': 'Known to rarely fail with a random shift',
        }

    def test_correct_atom_linking(self, with_symmetry_constraints, first_case):
        """Check that atoms were linked as expected."""
        slab, _, info = with_symmetry_constraints()
        with may_fail(first_case, self._known_invalid_linking):
            for atom_n, linked in info.symmetry.link_groups.items():
                atom = slab.atlist.get(atom_n)
                assert set(at.num for at in atom.linklist) == linked

    _known_invalid_constrained_group = {
        'rcm-ucell=rectangular': 'Known to sometimes fail with a random shift',
        'rcm-ucell=square': 'Known to sometimes fail with a random shift',
        'pm_10': 'Known to sometimes fail with a random shift',
        'cm_1m1': 'Known to sometimes fail with a random shift',
        'hex_cm_10': 'Known to sometimes fail with a random shift',
        'hex_cm_01': 'Known to sometimes fail with a random shift',
        'hex_cm_1m1': 'Known to sometimes fail with a random shift',
        'hex_cm_11': 'Known to sometimes fail with a random shift',
        'hex_cm_21': 'Known to sometimes fail with a random shift',
        'hex_cm_12': 'Known to sometimes fail with a random shift',
        'hex_cmm_10': 'Known incorrect plane group p2',
        'hex_cmm_01': 'Known incorrect plane group p2',
        'square_pm_10': 'Known to often fail with a random shift',
        'square_pm_01': 'Known to often fail with a random shift',
        'square_cm_11': 'Known to often fail with a random shift',
        'square_cm_1m1': 'Known to often fail with a random shift',
        'poscar_diamond': 'Known incorrect plane group pm instead of rcm',
        'poscar_sb_si_111': 'Known to rarely fail with a random shift',
        'poscar_sto110_4x1': 'Known to sometimes fail with a random shift',
        }

    def test_correct_group_after_constraint(self, with_symmetry_constraints,
                                            first_case):
        """Test that applying constraints maintains the group correct."""
        slab, param, info = with_symmetry_constraints()
        if not info.symmetry.hermann:
            pytest.skip('No symmetry information available')
        symmetry.findSymmetry(slab, param)
        with may_fail(first_case, self._known_invalid_constrained_group):
            assert hermann(slab.planegroup) == info.symmetry.hermann

    _known_do_not_preserve = {
        'rcm': 'Fails in all directions due to missing centre link',
        'rcmm': 'Fails in all directions due to missing centre link',
        }

    @parametrize(displacement=XYZ_DISPLACEMENTS)
    def test_preserve_displaced_symmetry(self, displacement, displace_atoms,
                                         with_symmetry_constraints,
                                         first_case):
        """Ensure that displacements conserve symmetry."""
        slab, param, info, *_ = with_symmetry_constraints(random_shifts=False)
        displaced_slab = copy.deepcopy(slab)
        displace_atoms(displaced_slab, param, info, displacement)
        with may_fail(first_case, self._known_do_not_preserve, strict=True):
            assert displaced_slab.foundplanegroup == slab.foundplanegroup


class TestSlabSymmetrization:
    """Displace atoms at random, and check the effect of symmetrization."""

    @fixture
    @parametrize_with_cases('args', cases=get_cases('all'),
                            filter=~id_has_suffix('p1'))
    def rattled_slab(self, args):
        """Return a slab with atoms slightly misplaced."""
        slab, param, info, *rest = args
        symmetry.findSymmetry(slab, param)
        symmetry.enforceSymmetry(slab, param)
        rattled_slab, param, info = duplicate_all(slab, param, info)
        eps = 0.05
        param.SYMMETRY_EPS = param.SYMMETRY_EPS.from_value(2 * 2**0.5 * eps)
        displacements = np.random.uniform(-eps, eps, (slab.n_atoms, 2))
        for atom, delta in zip(rattled_slab, displacements):
            atom.cartpos[:2] += delta
        rattled_slab.collapse_cartesian_coordinates()
        return (rattled_slab, slab, param, info, *rest)

    _known_incorrect_rattled = {
        'pmm-ucell=square': 'Often reduced to pm',
        'pgg-ucell=rectangular': 'Often reduced to pg',
        'pgg-ucell=square': 'Often reduced to pg',
        'rcm-ucell=square': 'Sometimes reduced to pm',
        'rcm-ucell=rectangular': 'Sometimes reduced to pm',
        'rcmm-ucell=square': 'Often reduced to, e.g., pmg',
        'rcmm-ucell=rectangular': 'Often reduced to, e.g., rcm',
        'pmm-ucell=rectangular': 'Sometimes reduced to pm',
        'pmg-ucell=rectangular': 'Sometimes identified as rcm',
        'cmm': 'Sometimes reduced to cm',
        'p2': 'Sometimes reduced to p1',
        'cm_1m1': 'Sometimes reduced to p1',
        'p4': 'Sometimes reduced to p1',
        'p4m': 'Sometimes reduced to cm',
        'p4g': 'Sometimes reduced to cm',
        'p3': 'Sometimes reduced to p1',
        'p31m': 'Often reduced to cm',
        'p3m1': 'Often reduced to cm',
        'p6': 'Sometimes reduced to p2/p1',
        'p6m': 'Sometimes reduced to cm',
        'hex_cm_1m1': 'Often reduced to p1',
        'hex_cm_11': 'Often reduced to p1',
        'hex_cm_10': 'Often reduced to p1',
        'hex_cm_01': 'Often reduced to p1',
        'hex_cm_21': 'Often reduced to p1',
        'hex_cm_12': 'Often reduced to p1',
        'hex_cmm_11': 'Often reduced to cm',
        'hex_cmm_10': 'Known invalid plane group p2. May be correct cmm here',
        'hex_cmm_01': 'Known invalid plane group p2. May be correct cmm here',
        'hex_p2': 'Often reduced to p1',
        'square_cm_11': 'Sometimes reduced to pm',
        'square_cm_1m1': 'Often reduced to pm',
        'square_pmm': 'Sometimes reduced to pm',
        'square_cmm': 'Sometimes reduced to cm',
        'poscar_36carbon_atoms_p6m': 'Often reduced to cmm',
        'poscar_ag100': 'Often reduced from p4m to cm',
        'poscar-Al2O3_NiAl(111)_cHole_20061025' : (
            'Sometimes reduced to p1 from p3'
            ),
        'poscar-Cu2O(111)_1x1_surplus_oxygen': (
            'Sometimes reduced to cm from p3m1'
            ),
        'poscar_diamond': 'Known invalid group pm. May be correct rcm here',
        'poscar-Fe3O4_SCV': 'Sometimes reduced to cm from cmm',
        'poscar_fe3o4_001_cod': 'Sometimes reduced to cm/p1 from cmm',
        'poscar_mgo': 'Sometimes reduced to cmm from p4m',
        'poscar-ru': 'Sometimes reduced to cm from p3m1',
        'poscar_sb_si_111': 'Sometimes symmetry increased to rcm from pm',
        'poscar-SiC_H': 'Often reduced to cm from p3m1',
        'poscar-TiO2': 'Sometimes identified as pmg instead of pmm',
        'poscar_tio2_small': 'Sometimes reduced to pm from pmm',
        'infoless_poscar-36C_p6m': 'Sometimes reduced to cm(m)',
        'infoless_poscar-Ag(100)': 'Often reduced from p4m to cm',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025' : (
            'Sometimes reduced to p1 from p3'
            ),
        'infoless_poscar-Cu2O_111': 'Sometimes reduced to cm from p3m1',
        'infoless_poscar-Cu2O(111)_1x1_surplus_oxygen': (
            'Sometimes reduced to cm from p3m1'
            ),
        'infoless_poscar-diamond': 'Invalid pm. May be correct rcm here',
        'infoless_poscar-Fe3O4_SCV': 'Sometimes reduced to cm from cmm',
        'infoless_poscar-Fe3O4_(001)_cod1010369': (
            'Sometimes reduced to cm from cmm'
            ),
        'infoless_poscar-graphene': 'Sometimes reduced to pmg from pmm',
        'infoless_poscar-In2O3_(111)': 'Sometimes reduced to p1 from p3',
        'infoless_poscar-Ir(100)-(2x1)-O': 'Sometimes misidentified as pm',
        'infoless_poscar-MgO_cod_9006456': 'Sometimes reduced to cmm from p4m',
        'infoless_poscar-Ru(0001)-rt3Te': 'Sometimes reduced to cm from p3m1',
        'infoless_poscar-Sb_Si(111)_rect': (
            'Sometimes symmetry increased to rcm from pm'
            ),
        'infoless_poscar-SiC_H': 'Often reduced to cm from p3m1',
        'infoless_poscar-TiO2_supercell': 'Sometimes reduced to pm from pmm',
        'infoless_poscar-TiO2_small': 'Sometimes reduced to pm/p1 from pmm',
        'bulk_repeat_poscar-Cu2O_111': 'Sometimes reduced to cm from p3m1',
        'double_bulk-fe3o4': 'Often reduced to pm/p1 (from pmm)',
        'double_bulk-ru': 'Sometimes reduced to cm from p3m1',
        'fe3o4_bulk': 'Known invalid group pm. May be correct pmm here.',
        }

    def test_rattled_group(self, rattled_slab, first_case):
        """Check that the rattled slab has the right plane group."""
        rattled, original, param, info, *_ = rattled_slab
        symmetry.findSymmetry(rattled, param)
        with may_fail(first_case, self._known_incorrect_rattled):
            if info.symmetry.hermann:
                assert hermann(original.planegroup) == info.symmetry.hermann
            assert original.planegroup == rattled.planegroup


class TestBulkSymmetry:
    """Collection of bulk-symmetry finding tests."""

    @parametrize_with_cases('args', cases=get_cases('poscar'),
                            filter=id_has_suffix('_bulk'))
    def test_bulk_symmetry(self, args, first_case):
        """Assert the correct identification of bulk screws and glides."""
        bulk_slab, param, info, *_ = args
        case_id = _reconstruct_case_id(first_case)
        if case_id == 'double_bulk-fe3o4':
            pytest.xfail('Known BUG (in Slab or findBulkSymmetry?) '
                         'fixed in refactor_slab + better_symmetry')
        symmetry.findBulkSymmetry(bulk_slab, param)
        assert set(bulk_slab.bulk_screws) == info.bulk.screw_orders
        assert len(bulk_slab.bulk_glides) == info.bulk.n_glide_planes

    @parametrize_with_cases('args', cases=get_cases('poscar'),
                            has_tag=CaseTag.THICK_BULK)
    def test_thick_bulk_minimized(self, args, caplog, first_case):
        """Assert the correct identification of bulk screws and glides."""
        bulk_slab, param, *_ = args
        case_id = _reconstruct_case_id(first_case)
        if case_id == 'double_bulk-fe3o4':
            pytest.xfail('Known BUG (in Slab or findBulkSymmetry?) '
                         'fixed in refactor_slab + better_symmetry')
        with caplog.at_level(logging.DEBUG):
            symmetry.findBulkSymmetry(bulk_slab, param)
        assert 'could be reduced' in caplog.text


class TestSymmetryReduction:
    """Collection of tests for applying symmetry reductions."""

    @pytest.mark.xfail(raises=KeyError, reason='Not checked correctly')
    def test_raises_without_group_info(self, slab_p1):
        """Assert that missing plane group information raises exceptions."""
        slab, param, *_ = slab_p1
        with pytest.raises(RuntimeError) as exc:
            symmetry.setSymmetry(slab, param, 'p1')
        assert exc.match(r'.*unknown')

    @pytest.mark.xfail(reason='Does not check validity of group')
    def test_invalid_group(self, with_plane_group):
        """Assert that not providing a direction raises exceptions."""
        slab, param, *_ = with_plane_group()
        with pytest.raises(ValueError) as exc:
            symmetry.setSymmetry(slab, param, 'kt')
        assert exc.match(r'kt is not.*acceptable')

    def test_missing_direction(self, with_plane_group):
        """Assert that not providing a direction raises exceptions."""
        slab, param, *_ = with_plane_group()
        hermann_ = hermann(slab.foundplanegroup)
        if hermann_ in 'pm|pmg':
            need_direction = 'cm'
        elif hermann_ == 'p6m':
            need_direction = 'cmm'
        else:
            need_direction = 'pm'
        # with pytest.raises(MissingGroupDirectionError):
        symmetry.setSymmetry(slab, param, need_direction)

    @pytest.mark.xfail(reason='Does not check validity of subgroup')
    def test_invalid_subgroup(self, with_plane_group):
        """Make sure SlabSymmetrizer finds all atoms it should."""
        slab, param, *_ = with_plane_group()
        hermann_ = hermann(slab.foundplanegroup)
        invalid_group = 'p4' if hermann_ == 'p6m' else 'p6m'
        with pytest.raises(ValueError) as exc:
            symmetry.setSymmetry(slab, param, invalid_group)
        assert exc.match('.*not a valid symmetry reduction.*')

    @pytest.mark.xfail(reason='BUG: symmetry L724--789, missing [] at arrays',
                       raises=(ValueError, RuntimeWarning))
    def test_successful_reduction_pmg(self, slab_pmg, subtests):
        """Test that pmg can be reduced without giving a direction."""
        # Normally pm and pg require a direction. However we can
        # also not give them for the special case of pmg, as it
        # is clear which is which.
        slab, param, *_ = slab_pmg
        symmetry.findSymmetry(slab, param)
        oriplane_pmg = copy.deepcopy(slab.orisymplane)
        symmetry.setSymmetry(slab, param, 'pm')
        with subtests.test('pm'):
            assert hermann(slab.planegroup) == 'pm'
            assert all(slab.orisymplane.perp == oriplane_pmg.par)
        with subtests.test('pg'):
            symmetry.setSymmetry(slab, param, 'pg')
            assert hermann(slab.planegroup) == 'pg'
            assert all(slab.orisymplane.par == oriplane_pmg.par)

    _known_failed_reductions = {'cm_01', 'cm_21', 'cmm_01'}

    def test_successful_reduction_hex_rotation(self, slab_p6m, subtests):
        """Check unit-cell rotations when reducing p6m."""
        slab, param, *_ = slab_p6m
        symmetry.findSymmetry(slab, param)
        a_before, b_before = slab.ab_cell.T.copy()
        groups = ('cm[1 0]', 'cm[0 1]', 'cm[2 1]', 'cm[1 2]',
                  'cmm[1 0]', 'cmm[0 1]')
        # Need to get info about rotation angle. It's in the
        # TestInfo returned by CaseSimpleHexagonalSlabs methods
        hex_ = simple_slabs.CaseSimpleHexagonalSlabs()
        for group in groups:
            symmetry.setSymmetry(slab, param, group)
            a_after, b_after = slab.ab_cell.T
            # Convert group name into a method of hex_
            name = group.replace('[', '_').replace(']', '').replace(' ', '')
            method = getattr(hex_, f'case_hex_{name}')
            *_, info = method()
            with subtests.test(group):
                assert angle(a_before, b_before) == angle(a_after, b_after)
                if name in self._known_failed_reductions:
                    pytest.xfail('Known inconsistency between detection '
                                 'and symmetry reduction. Angles differ by '
                                 '180 degrees')
                assert angle(a_before, a_after) == info.symmetry.rotation
