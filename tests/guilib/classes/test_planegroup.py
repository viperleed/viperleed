"""Test functionality of PlaneGroup class.

Created: 2020-01-11
Author: Michele Riva
"""

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.guilib.classes import planegroup
from viperleed.guilib.classes.planegroup import PlaneGroup


_equal = {
    'valid string': (PlaneGroup('pmg[10]'), 'pmg[1 0 ] ', True),  # Fails!
    'invalid string': (PlaneGroup('pgg'), 'p-gg', False),
    'other': (PlaneGroup('cmm'), 1, False),
    'different group': (PlaneGroup('cmm'), PlaneGroup('p6'), False),
    }

@parametrize('group1,group2,expect', _equal.values(), ids=_equal)
def test_equal(group1, group2, expect):
    """Check expected outcome of equality test."""
    _eq = group1 == group2
    assert _eq is expect


def test_highest_symmetry(subtests):
    _high = PlaneGroup.highest_symmetry_for_shape
    with subtests.test('oblique'):
        assert _high('Oblique') == 'p2'
    with subtests.test('rect'):
        assert _high('Rectangular').same_operations('pmm')
    with subtests.test('square'):
        assert _high('Square').same_operations('p4m')
    with subtests.test('hex'):
        assert _high('Hexagonal') == 'p6m'


_init = {
    'p2': ('p2', 'p2', ('E', 'C2')),
    'cmm[10]': ('cMm[10]', 'cmm[1 0]', ('E', 'C2', 'M10', 'M12')),
    'from p4m': (PlaneGroup('p4m'), 'p4m',
                 ('E', 'C2', 'C4', 'Cm4', 'Mx', 'My', 'M45', 'Mm45')),
    }

@parametrize('init,expect_group,expect_ops', _init.values(), ids=_init)
def test_init(init, expect_group, expect_ops, subtests):
    """Check correct creation of a PlaneGroup instance."""
    group = PlaneGroup(init)
    with subtests.test('group'):
        assert group.group == expect_group
    ops = group.operations(include_3d=True)  # None defined anyway
    with subtests.test('no. operations'):
        assert len(ops) == len(expect_ops)
    with subtests.test('matrices'):
        expect_matrices = {getattr(planegroup, op) for op in expect_ops}
        assert set(ops) == expect_matrices


def test_is_valid_group(subtests):
    """Check correct identification of group validity for a cell shape."""
    with subtests.test('p4m, square'):
        assert PlaneGroup.is_valid_group('p4m', cell_shape='Square')
    with subtests.test('p3m1, square'):
        assert not PlaneGroup.is_valid_group('p31m', 'Rectangular')
    with subtests.test('invalid, hex'):
        assert not PlaneGroup.is_valid_group('invalid', 'Hexagonal')


class TestGroupsCompatibleWith:
    """Tests for listing groups for a given shape and operations."""

    def test_shape_only(self):
        """Check groups compatible with a cell shape."""
        assert 'p4' in PlaneGroup.groups_compatible_with('Square')

    def test_include_operations(self):
        operations = PlaneGroup('cmm').operations()
        groups = PlaneGroup.groups_compatible_with('Square', operations)
        assert 'cm[1 -1]' in groups

    def test_p3m1_rt3(self):
        """Check correct groups for a rt3xrt3 on a p3m1 bulk."""
        rt3 = (2, 1), (-1, 1)
        p3m1 = PlaneGroup('p3m1')
        ops = [op.round().astype(int) for op in p3m1.transform(rt3)]
        assert 'p31m' in PlaneGroup.groups_compatible_with('Hexagonal', ops)


class TestProperties:
    """Collection of tests for Plane @property."""

    def test_subgroups(self):
        """Check correct subgroup identification."""
        group = PlaneGroup('p6')
        assert 'p3' in group.subgroups


class TestRaises:
    """Collection of tests for checking PlaneGroup complaints."""

    _init = {
        'not a group string': ('invalid', ValueError),
        'not a string': (1, TypeError),
        'funny characters': ('.in,valid', ValueError),
        'cannot have direction': ('p3m1[1 0]', ValueError),
        'must have direction': ('pm', planegroup.MissingDirectionError),
        }

    @parametrize('group,exc', _init.values(), ids=_init)
    def test_invalid_group(self, group, exc):
        """Check complaints when initializing with an invalid object."""
        with pytest.raises(exc):
            PlaneGroup(group)

    _screws_glides = {
        'string': ('invalid', '', ValueError, 'Invalid input'),
        'matrices wrong shape': ([np.eye(3)]*3, '', ValueError, 'shape'),
        'non-int matrices': ([np.eye(2)+0.1]*4, '', ValueError, 'integer'),
        'non-det-1 matrices': ([np.eye(2)+1]*2, '', ValueError, 'determinant'),
        'wrong type': (1, '', TypeError, 'Must be'),
        'wrong screw rotation': ('r(5)', '', ValueError, 'rotation order'),
        'wrong glide direction': ('m([  1,7])', 'Square',
                                  ValueError, 'direction'),
        'glide string no shape': ('m([1,0], [-2, 1])', None,
                                  TypeError, 'cell shape'),
        }

    @parametrize('screws_glides,shape,exc,exc_re', _screws_glides.values(),
                 ids=_screws_glides)
    def test_set_screws_glides_invalid(self, screws_glides, shape,
                                       exc, exc_re):
        """Check complaints for invalid inputs to set_screws_glides()."""
        group = PlaneGroup(group='p1')
        with pytest.raises(exc) as exc_info:
            group.set_screws_glides(screws_glides, cell_shape=shape)
        assert exc_info.match(exc_re)

    _transform = {
        'transform shape': (np.ones((2, 3)), None),
        'inverse shape': (np.eye(2), np.ones((2, 3))),
        'inverse not inverse': (np.eye(2), 2*np.eye(2)),
        }

    @parametrize('transform,inverse', _transform.values(), ids=_transform)
    def test_transform_invalid(self, transform, inverse):
        """Check complaints for invalid inputs to .transform()."""
        group = PlaneGroup(group='p1')
        with pytest.raises(ValueError):
            group.transform(transform, inverse)

    _groups_compatible = {
        'operations TypeError': ('Square', 1, TypeError),
        'cell shape': ('invalid', None, ValueError),
        'operations shape': ('Hexagonal', [np.eye(3)] * 3, ValueError),
        }

    @parametrize('shape,ops,exc', _groups_compatible.values(),
                 ids=_groups_compatible)
    def test_groups_compatible_invalid(self, shape, ops, exc):
        """Check complaints for invalid inputs to .groups_compatible_with()."""
        with pytest.raises(exc):
            PlaneGroup.groups_compatible_with(shape, operations=ops)


class TestSameOperations:
    """Tests to verify whether groups have equivalent symmetry operations."""

    _same = {
        'same group': (PlaneGroup('p1'), 'p1', True),
        'different groups': (PlaneGroup('p2'), 'p1', False),
        'different equivalent groups': (PlaneGroup('p4m'), 'p4g', True),
        'other object': (PlaneGroup('p2'), None, False),
        }

    @parametrize('group1,group2,result', _same.values(), ids=_same)
    def test_exclude_3d(self, group1, group2, result):
        """Check correct identification of identical operations."""
        assert group1.same_operations(group2) is result

    def test_include_3d(self):
        """Check that screws/glides are considered among the operations."""
        group1 = PlaneGroup('p3')
        group2 = PlaneGroup('p3')
        group1.set_screws_glides('r(2,6) m([1,2])', 'Hexagonal')
        assert not group1.same_operations(group2, include_3d=True)


class TestScrewsGlides:
    """Collection of tests for setting/checking the 3D screws/glides."""

    @parametrize(none=(None, 'none', '', tuple(), False))
    def test_none(self, none):
        """Check correct clearing of screws/glides."""
        group = PlaneGroup(group='p1')
        group.set_screws_glides(new_screws_glides=none)
        assert group.screws_glides == ()

    def test_set_from_matrices(self):
        """Check correct assignment of screws/glides from matrices."""
        group = PlaneGroup(group='p1')
        matrices = planegroup.E, planegroup.C2, planegroup.Cm6, planegroup.E
        group.set_screws_glides(new_screws_glides=matrices)
        assert len(group.screws_glides) == 3

    _string = {
        'only screws 1': ('r(2, 3)', ('C2', 'C3', 'Cm3')),
        'only screws 2': ('r(2), R(6)', ('C2', 'C3', 'Cm3', 'C6', 'Cm6')),
        'only glides 1': ('m( [  1, 2], [-1, -2], [1,0])', ('M12', 'Mx')),
        'only glides 2': ('M([1,2]) m([1,-1], [0,1])', ('M12', 'M1m1', 'My')),
        }

    @parametrize('screws_glides,expect', _string.values(), ids=_string)
    def test_set_from_string(self, screws_glides, expect):
        """Check correct assignment of screws/glides from string."""
        group = PlaneGroup('p1')
        group.set_screws_glides(screws_glides, 'Square')
        expect_ops = {getattr(planegroup, op) for op in expect}
        assert len(group.screws_glides) == len(expect)
        assert set(group.screws_glides) == expect_ops
