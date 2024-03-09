"""Tests for the lattice2d module of viperleed.guilib.classes.

Created: 2020-01-11
Author: Michele Riva
"""

import pytest
from pytest_cases import parametrize
import numpy as np

from viperleed.guilib.classes.lattice2d import Lattice2D


_HEX_ACUTE = (1, 0), (0.5, 3**0.5/2)
_HEX_OBTUSE = (1, 0), (-0.5, 3**0.5/2)
_OBLIQUE_ACUTE = (1, 0), (1.2, 3.5)
_OBLIQUE_OBTUSE = (1, 0), (-1.2, 3.5)
_RECT = (1, 0), (0, 3.5)
_RHOMBIC_ACUTE = (2, 1), (2, -1)
_RHOMBIC_OBTUSE = (-2, 1), (2, 1)
_SQUARE = (1, 0), (0, 1)


class TestProperties:
    """Collection of tests for accessing properties."""

    _cell = {
        _HEX_ACUTE: 'Hexagonal', _HEX_OBTUSE: 'Hexagonal',
        _OBLIQUE_ACUTE: 'Oblique', _OBLIQUE_OBTUSE: 'Oblique',
        _RECT: 'Rectangular',
        _RHOMBIC_ACUTE: 'Rhombic', _RHOMBIC_OBTUSE: 'Rhombic',
        _SQUARE: 'Square',
        }

    @parametrize('basis,shape', _cell.items(), ids=(str(c) for c in _cell))
    def test_cell_shape(self, basis, shape):
        """Check correct identification of the cell shape."""
        assert Lattice2D(basis).cell_shape == shape

    _group = {
        None: 'p1',  # Default
        'pgg': 'pgg'
        }

    @parametrize('group,expect', _group.items(), ids=_group)
    def test_group_valid(self, group, expect):
        """Check assignment of a valid plane group."""
        lattice = Lattice2D(_RECT)
        lattice.group = group
        assert lattice.group == expect

    _params = {
        (_HEX_ACUTE, 'real'): (1, 1, 120),
        (_HEX_OBTUSE, 'real'): (1, 1, 120),
        (_HEX_ACUTE, 'reciprocal'): (1, 1, 60),
        (_HEX_OBTUSE, 'reciprocal'): (1, 1, 60),
        (_OBLIQUE_ACUTE, 'real'): (1, 3.7, 71.075),
        (_OBLIQUE_ACUTE, 'reciprocal'): (1, 3.7, 71.075),
        (_OBLIQUE_OBTUSE, 'real'): (1, 3.7, 108.925),
        (_OBLIQUE_OBTUSE, 'reciprocal'): (1, 3.7, 108.925),
        (_RECT, 'real'): (1, 3.5, 90),
        (_RECT, 'reciprocal'): (1, 3.5, 90),
        (_RHOMBIC_ACUTE, 'real'): (5**0.5, 5**0.5, 126.87),
        (_RHOMBIC_OBTUSE, 'real'): (5**0.5, 5**0.5, 126.87),
        (_RHOMBIC_ACUTE, 'reciprocal'): (5**0.5, 5**0.5, 53.13),
        (_RHOMBIC_OBTUSE, 'reciprocal'): (5**0.5, 5**0.5, 53.13),
        (_SQUARE, 'real'): (1, 1, 90),
        (_SQUARE, 'reciprocal'): (1, 1, 90),
        }

    @parametrize('args,expect', _params.items(), ids=(str(p) for p in _params))
    def test_lattice_parameters(self, args, expect):
        """Check correct values of vector lengths and angle."""
        basis, space = args
        lat = Lattice2D(basis, space=space)
        assert lat.lattice_parameters == pytest.approx(expect, abs=1e-3)
        assert np.linalg.det(basis) == np.linalg.det(lat.basis)


class TestRaises:
    """Collection of tests for complaints by Lattice2D objects."""

    _init = {
        'space': (_SQUARE, {'space': 'invalid_space'}, ValueError),
        'limit': (_SQUARE, {'limit': 'not_a_number'}, TypeError),
        'basis': ('not_a_list', {}, TypeError),
        'basis shape': (np.eye(3), {}, ValueError),
        'basis singular': (np.zeros((2, 2)), {}, ValueError),
        'group for shape': (_HEX_OBTUSE, {'group': 'p4'}, ValueError),
        'not a group': (_HEX_OBTUSE, {'group': 'invalid_group'}, ValueError),
        'group type': (_HEX_OBTUSE, {'group': 123}, ValueError),
        }

    @parametrize('basis,kwargs,exc', _init.values(), ids=_init)
    def test_init_invalid(self, basis, kwargs, exc):
        """Check complaints for making a Lattice2D with invalid arguments."""
        with pytest.raises(exc):
            Lattice2D(basis, **kwargs)

    _attr = {
        'n_beams': Lattice2D(_SQUARE, space='real'),
        'n_points': Lattice2D(_SQUARE, space='reciprocal')
        }

    @parametrize('attr,lattice', _attr.items(), ids=_attr)
    def test_attribute_error(self, lattice, attr):
        """Check complaints when accessing a property for the wrong space."""
        with pytest.raises(AttributeError):
            _ = getattr(lattice, attr)

