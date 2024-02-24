"""Test functionality of PlaneGroup class.

Created: 2020-01-11
Author: Michele Riva
"""

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.guilib.classes.planegroup import PlaneGroup


def test_planegroup_p2():
    g = PlaneGroup('p2')
    assert g.group == 'p2'

    E = [[1, 0], [0, 1]]
    C2 = [[-1, 0], [0, -1]]
    ops = g.operations()
    assert len(ops) == 2
    assert np.array_equal(ops[0], E)
    assert np.array_equal(ops[1], C2)


class TestPlaneGroupRaises:
    """Collection of tests for checking PlaneGroup complaints."""

    _init = {
        'not a group string': ('invalid', ValueError),
        'not a string': (1, TypeError),
        }

    @parametrize('group,exc', _init.values(), ids=_init)
    def test_invalid_group(self, group, exc):
        """Check complaints when initializing with an invalid object."""
        with pytest.raises(exc):
            PlaneGroup(group)

    _screws_glides = {
        'string': ('invalid', '', ValueError, 'Invalid input'),
        'matrices wrong shape': ([np.eye(3)] * 3, '', ValueError, 'shape'),
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
