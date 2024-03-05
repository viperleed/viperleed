"""Tests for the lattice2d module of viperleed.guilib.classes.

Created: 2020-01-11
Author: Michele Riva
"""

import sys
import os

import pytest
import numpy as np

from viperleed.guilib import base
from viperleed.guilib.classes.lattice2d import Lattice2D


def test_lattice_init():
    b = [[1, 0], [0, 2]]
    g = 'pgg'
    latt = Lattice2D(b, group=g)
    assert np.array_equal(latt.basis, b)
    assert latt.group.group == g
    assert latt.cell_shape == 'Rectangular'
    assert latt.n_beams == 3

def test_lattice_wrong_basis_input():
    for b in ('', dict(), set(), 3, np.array([]), (1, 0), 2.5):
        with pytest.raises((TypeError, ValueError)) as exc:
            Lattice2D(b)
            if base.check_type(b, 'arraylike'):
                assert exc == ValueError
            else:
                assert exc == TypeError


def test_lattice_wrong_space_input():
    for space in ('', dict(), set(), 3, np.array([]), (1, 0), 2.5):
        with pytest.raises((TypeError, ValueError)) as exc:
            Lattice2D([[1, 0], [0, 1]], space=space)
            if base.check_type(space, 'str'):
                assert exc == ValueError
            else:
                assert exc == TypeError


def test_lattice_wrong_limit_input():
    for space in ('', dict(), set(), np.array([]), (1, 0), 2.5):
        with pytest.raises((TypeError, ValueError)) as exc:
            Lattice2D([[1, 0], [0, 1]], space=space)
            if base.check_type(space, 'str'):
                assert exc == ValueError
            else:
                assert exc == TypeError
