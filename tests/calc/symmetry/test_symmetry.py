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

from ...helpers import check_draws
from ...helpers import duplicate_all
from ..poscar_slabs import make_poscar_ids
from ..tags import CaseTag
from . import simple_slabs
from .conftest import get_cases


SEEDED_SHIFT_FAILS = 'Fails with some of the seeded origin shifts'
CENTRING_NOT_LINKED = 'Known bug: centring-related atoms are not linked'


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
def may_fail(case, known_failures, strict=False, draws=None):
    """XFAIL if case is known to be problematic.

    Known failures of random-input tests refer to the default seeds:
    pass the `draws` of the test to just check if they are fresh.
    """
    reason = known_failures.get(case.id, '')
    if not reason:
        reason = known_failures.get(_reconstruct_case_id(case), '')
    if draws is not None and draws.fresh:
        reason = ''
        strict = False
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

    def test_any_plane_group_found(self, with_plane_group, draws):
        """Ensure that symmetry detection is successful."""
        def _check(draw):
            slab, *_ = with_plane_group(draw)
            assert slab.foundplanegroup != 'unknown'

        check_draws(draws, _check)

    known_incorrect_groups = {
        'hex_cmm_10': 'Known incorrect plane group p2',
        'hex_cmm_01': 'Known incorrect plane group p2',
        'poscar_diamond': 'Known incorrect plane group pm instead of rcm',
        }

    def test_correct_plane_group(self, with_plane_group, first_case, draws):
        """Check the correct identification of the plane group."""
        def _check(draw):
            slab, *_, info = with_plane_group(draw)
            if not info.symmetry.hermann:
                pytest.skip('No symmetry information available')
            assert hermann(slab.planegroup) == info.symmetry.hermann

        with may_fail(first_case,
                      self.known_incorrect_groups,
                      strict=True,
                      draws=draws):
            check_draws(draws, _check)

    _known_incorrect_rotations = {
        'hex_cmm_10': 'BUG: no rotation applied',
        'hex_cmm_01': 'BUG: no rotation applied',
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
                # The warning may be any of the logged messages
                assert any(
                    re_match(r'.*unit.*cell.*change.*higher.*symmetry', msg)
                    for msg in caplog.messages
                    )
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
        'cm_1m1': SEEDED_SHIFT_FAILS,
        'pm_10': SEEDED_SHIFT_FAILS,
        'rcm-ucell=rectangular': SEEDED_SHIFT_FAILS,
        'rcm-ucell=square': SEEDED_SHIFT_FAILS,
        'hex_cm_10': SEEDED_SHIFT_FAILS,
        'hex_cm_01': SEEDED_SHIFT_FAILS,
        'hex_cm_21': SEEDED_SHIFT_FAILS,
        'hex_cm_12': SEEDED_SHIFT_FAILS,
        'hex_cm_1m1': SEEDED_SHIFT_FAILS,
        'hex_cmm_10': 'Known incorrect plane group p2',
        'hex_cmm_01': 'Known incorrect plane group p2',
        'square_cm_1m1': SEEDED_SHIFT_FAILS,
        'square_pm_10': SEEDED_SHIFT_FAILS,
        }

    def test_correct_free_directions(self, with_symmetry_constraints,
                                     first_case, draws):
        """Check the correct identification of atom constraints."""
        def _check(draw):
            slab, _, info = with_symmetry_constraints(draw)
            if not info.symmetry.link_groups:
                pytest.skip('Not enough symmetry information available')
            for atom in slab:
                if atom.num in info.symmetry.on_planes:
                    assert isinstance(atom.freedir, np.ndarray)
                elif atom.num in info.symmetry.on_axes:
                    assert not atom.freedir
                else:
                    assert atom.freedir == 1

        with may_fail(first_case,
                      self._known_invalid_constraints,
                      strict=True,
                      draws=draws):
            check_draws(draws, _check)

    _known_invalid_linking = {
        'rcm-ucell=rectangular': CENTRING_NOT_LINKED,
        'rcm-ucell=square': CENTRING_NOT_LINKED,
        'rcmm-ucell=rectangular': CENTRING_NOT_LINKED,
        'rcmm-ucell=square': CENTRING_NOT_LINKED,
        'cm_1m1': SEEDED_SHIFT_FAILS,
        'hex_cm_10': SEEDED_SHIFT_FAILS,
        'hex_cm_01': SEEDED_SHIFT_FAILS,
        'hex_cm_21': SEEDED_SHIFT_FAILS,
        'hex_cm_12': SEEDED_SHIFT_FAILS,
        'hex_cm_1m1': SEEDED_SHIFT_FAILS,
        'hex_cmm_10': 'Known incorrect plane group p2',
        'hex_cmm_01': 'Known incorrect plane group p2',
        'square_cm_11': SEEDED_SHIFT_FAILS,
        'square_cm_1m1': SEEDED_SHIFT_FAILS,
        'square_pg_10': SEEDED_SHIFT_FAILS,
        'square_pm_10': SEEDED_SHIFT_FAILS,
        'pm_10': SEEDED_SHIFT_FAILS,
        'pg_10': SEEDED_SHIFT_FAILS,
        }

    def test_correct_atom_linking(self, with_symmetry_constraints,
                                  first_case, draws):
        """Check that atoms were linked as expected."""
        def _check(draw):
            slab, _, info = with_symmetry_constraints(draw)
            for atom_n, linked in info.symmetry.link_groups.items():
                atom = slab.atlist.get(atom_n)
                assert set(at.num for at in atom.linklist) == linked

        with may_fail(first_case,
                      self._known_invalid_linking,
                      strict=True,
                      draws=draws):
            check_draws(draws, _check)

    _known_invalid_constrained_group = {
        'pm_10': SEEDED_SHIFT_FAILS,
        'hex_cm_01': SEEDED_SHIFT_FAILS,
        'square_cm_1m1': SEEDED_SHIFT_FAILS,
        'poscar_sto110_4x1': SEEDED_SHIFT_FAILS,
        'hex_cmm_10': 'Known incorrect plane group p2',
        'hex_cmm_01': 'Known incorrect plane group p2',
        'poscar_diamond': 'Known incorrect plane group pm instead of rcm',
        }

    def test_correct_group_after_constraint(self, with_symmetry_constraints,
                                            first_case, draws):
        """Test that applying constraints maintains the group correct."""
        def _check(draw):
            slab, param, info = with_symmetry_constraints(draw)
            if not info.symmetry.hermann:
                pytest.skip('No symmetry information available')
            symmetry.findSymmetry(slab, param)
            assert hermann(slab.planegroup) == info.symmetry.hermann

        with may_fail(first_case,
                      self._known_invalid_constrained_group,
                      strict=True,
                      draws=draws):
            check_draws(draws, _check)

    _known_do_not_preserve = {
        'rcm': 'Fails in all directions due to missing centre link',
        'rcmm': 'Fails in all directions due to missing centre link',
        }

    @parametrize(displacement=XYZ_DISPLACEMENTS)
    def test_preserve_displaced_symmetry(self, displacement, displace_atoms,
                                         with_symmetry_constraints,
                                         first_case, draws):
        """Ensure that displacements conserve symmetry."""
        # No random shift is applied, so all draws would be identical
        slab, param, info, *_ = with_symmetry_constraints(
            draws[0],
            random_shifts=False,
            )
        displaced_slab = copy.deepcopy(slab)
        displace_atoms(displaced_slab, param, info, displacement)
        with may_fail(first_case, self._known_do_not_preserve, strict=True):
            assert displaced_slab.foundplanegroup == slab.foundplanegroup


class TestSlabSymmetrization:
    """Displace atoms at random, and check the effect of symmetrization."""

    @fixture
    @parametrize_with_cases('args', cases=get_cases('all'),
                            filter=~id_has_suffix('p1'))
    def rattled_slab(self, args, draws):
        """Return the slab, and a rattled copy of it for each draw.

        Returns
        -------
        slab, rattled, param, info, *rest
            `rattled` maps each draw index to a tuple
            (rattled slab, Rparams, info) for that draw.
        """
        slab, param, info, *rest = args
        symmetry.findSymmetry(slab, param)
        symmetry.enforceSymmetry(slab, param)
        eps = 0.05
        rattled = {}
        for draw in draws:
            rattled_slab, draw_param, draw_info = duplicate_all(slab,
                                                                param,
                                                                info)
            draw_param.SYMMETRY_EPS = draw_param.SYMMETRY_EPS.from_value(
                2 * 2**0.5 * eps
                )
            displacements = draw.rng.uniform(-eps, eps, (slab.n_atoms, 2))
            for atom, delta in zip(rattled_slab, displacements):
                atom.cartpos[:2] += delta
            rattled_slab.collapse_cartesian_coordinates()
            rattled[draw.index] = (rattled_slab, draw_param, draw_info)
        return (slab, rattled, param, info, *rest)

    _known_incorrect_rattled = {
        'p31m': 'Seeded rattling gives cm instead of p31m',
        'p4m': 'Seeded rattling gives cm instead of p4m',
        'p6m': 'Seeded rattling gives cm instead of p6m',
        'pgg-ucell=rectangular': 'Seeded rattling gives pg instead of pgg',
        'pmg-ucell=square': 'Seeded rattling gives rcm instead of pmg',
        'pmm-ucell=square': 'Seeded rattling gives p1 or pm instead of pmm',
        'rcm-ucell=rectangular': 'Seeded rattling gives pm instead of rcm',
        'rcmm-ucell=square': 'Seeded rattling gives pg instead of rcmm',
        'hex_cm_01': 'Seeded rattling gives p1 instead of cm',
        'hex_cm_10': 'Seeded rattling gives p1 instead of cm',
        'hex_cm_11': 'Seeded rattling gives p1 instead of cm',
        'hex_cm_1m1': 'Seeded rattling gives p1 instead of cm',
        'hex_cm_21': 'Seeded rattling gives p1 instead of cm',
        'hex_cmm_10': 'Known invalid plane group p2. May be correct cmm here',
        'hex_cmm_01': 'Known invalid plane group p2. May be correct cmm here',
        'hex_p2': 'Seeded rattling gives p1 instead of p2',
        'square_cm_11': 'Seeded rattling gives pm instead of cm',
        'square_cm_1m1': 'Seeded rattling gives pm instead of cm',
        'square_cmm': 'Seeded rattling gives cm instead of cmm',
        'square_pmm': 'Seeded rattling gives pm instead of pmm',
        'poscar-Al2O3_NiAl(111)_cHole_20061025': (
            'Seeded rattling gives p1 instead of p3'
            ),
        'poscar-Cu2O(111)_1x1_surplus_oxygen': (
            'Seeded rattling gives cm or p1 instead of p3m1'
            ),
        'poscar_diamond': 'Known invalid group pm. May be correct rcm here',
        'poscar_mgo': 'Seeded rattling gives rcm instead of p4m',
        'poscar-SiC_H': 'Seeded rattling gives cm instead of p3m1',
        'poscar_ru0001_rt3_te': (
            'Seeded rattling gives cm or p3m1 instead of p31m'
            ),
        'poscar_sb_si_111': 'Seeded rattling gives rcm instead of pm',
        'poscar_tio2_small': 'Seeded rattling gives pm instead of pmm',
        'infoless_poscar-36C_p6m': 'Seeded rattling gives cmm instead of p6m',
        'infoless_poscar-Al2O3_NiAl(111)_cHole_20061025': (
            'Seeded rattling gives p1 instead of p3'
            ),
        'infoless_poscar-Cu2O(111)_1x1_surplus_oxygen': (
            'Seeded rattling gives cm instead of p3m1'
            ),
        'infoless_poscar-diamond': 'Invalid pm. May be correct rcm here',
        'infoless_poscar-Fe3O4_(001)_cod1010369': (
            'Seeded rattling gives cm instead of cmm'
            ),
        'infoless_poscar-graphene': 'Seeded rattling gives pmg instead of pmm',
        'infoless_poscar-MgO_cod_9006456': (
            'Seeded rattling gives cmm instead of p4m'
            ),
        'infoless_poscar-Pt25Rh75(100)-p(3x1)-O': (
            'Seeded rattling gives pm instead of pmm'
            ),
        'infoless_poscar-Ru(0001)-rt3Te': (
            'Seeded rattling gives cm, p1 or p3m1 instead of p31m'
            ),
        'infoless_poscar-Sb_Si(111)_rect': (
            'Seeded rattling gives rcm instead of pm'
            ),
        'infoless_poscar-SiC_H': (
            'Seeded rattling gives cm or p3 instead of p3m1'
            ),
        'bulk_repeat_poscar-Cu2O_111': (
            'Seeded rattling gives cm or p3 instead of p3m1'
            ),
        'double_bulk-fe3o4': 'Seeded rattling gives pm instead of pmm',
        'double_bulk-ru': 'Seeded rattling gives cm instead of p3m1',
        'fe3o4_bulk': 'Seeded rattling gives pm instead of pmm',
        }

    def test_rattled_group(self, rattled_slab, first_case, draws):
        """Check that the rattled slab has the right plane group."""
        original, rattled_all, _, info, *_ = rattled_slab

        def _check(draw):
            rattled, param, _ = rattled_all[draw.index]
            symmetry.findSymmetry(rattled, param)
            if info.symmetry.hermann:
                assert hermann(original.planegroup) == info.symmetry.hermann
            assert original.planegroup == rattled.planegroup

        with may_fail(first_case,
                      self._known_incorrect_rattled,
                      strict=True,
                      draws=draws):
            check_draws(draws, _check)


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

    @pytest.mark.xfail(reason='Does not check validity of group',
                       strict=True)
    def test_invalid_group(self, with_plane_group, draws):
        """Assert that not providing a direction raises exceptions."""
        def _check(draw):
            slab, param, *_ = with_plane_group(draw)
            with pytest.raises(ValueError) as exc:
                symmetry.setSymmetry(slab, param, 'kt')
            assert exc.match(r'kt is not.*acceptable')

        check_draws(draws, _check)

    def test_missing_direction(self, with_plane_group, draws):
        """Assert that not providing a direction raises exceptions."""
        def _check(draw):
            slab, param, *_ = with_plane_group(draw)
            hermann_ = hermann(slab.foundplanegroup)
            if hermann_ in 'pm|pmg':
                need_direction = 'cm'
            elif hermann_ == 'p6m':
                need_direction = 'cmm'
            else:
                need_direction = 'pm'
            # with pytest.raises(MissingGroupDirectionError):
            symmetry.setSymmetry(slab, param, need_direction)

        check_draws(draws, _check)

    @pytest.mark.xfail(reason='Does not check validity of subgroup',
                       strict=True)
    def test_invalid_subgroup(self, with_plane_group, draws):
        """Make sure SlabSymmetrizer finds all atoms it should."""
        def _check(draw):
            slab, param, *_ = with_plane_group(draw)
            hermann_ = hermann(slab.foundplanegroup)
            invalid_group = 'p4' if hermann_ == 'p6m' else 'p6m'
            with pytest.raises(ValueError) as exc:
                symmetry.setSymmetry(slab, param, invalid_group)
            assert exc.match('.*not a valid symmetry reduction.*')

        check_draws(draws, _check)

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
