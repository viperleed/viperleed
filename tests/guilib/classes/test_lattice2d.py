"""Tests for the lattice2d module of viperleed.guilib.classes.

Created: 2020-01-11
Author: Michele Riva
"""

import pytest
from pytest_cases import parametrize
import numpy as np

from viperleed.guilib.classes.lattice2d import Lattice2D


_HEX_OBTUSE = (1, 0), (-0.5, 3**0.5/2)
_SQUARE = (1, 0), (0, 1)


class TestRaises:
    """Collection of tests for complaints by Lattice2D objects."""

    _init = {
        'space': (_SQUARE, {'space': 'invalid_space'}, ValueError),
        'limit': (_SQUARE, {'limit': 'not_a_number'}, TypeError),
        'basis': ('invalid_basis', {}, TypeError),
        'basis shape': (np.eye(3), {}, ValueError),
        'basis singular': (np.zeros((2, 2)), {}, ValueError),
        'group': (_HEX_OBTUSE, {'group': 'p4'}, ValueError),
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

