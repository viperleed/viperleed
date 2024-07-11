"""
Test functionality of Lattice class without loading the whole TLEEDMAP

Created: 2020-01-11
Author: Michele Riva
"""

import sys
import os 

tests_path = os.path.realpath(os.path.dirname(__file__))
base_path = os.path.realpath(os.path.join(tests_path, '..'))
for path in [tests_path, base_path]:
    if path not in sys.path:
        sys.path.append(path)


import pytest
import numpy as np

import viperleed.guilib.base as base


def test_lattice_init():
    b = [[1, 0], [0, 2]]
    g = 'pgg'
    latt = base.Lattice(b, group=g)
    assert np.array_equal(latt.basis, b)
    assert latt.group.group == g
    assert latt.cell_shape == 'Rectangular'
    assert latt.nbeams() == 3

def test_lattice_wrong_basis_input():
    for b in ('', dict(), set(), 3, np.array([]), (1, 0), 2.5):
        with pytest.raises((TypeError, ValueError)) as exc:
            base.Lattice(b)
            if base.checkType(b, 'arraylike'):
                assert exc == ValueError
            else:
                assert exc == TypeError


def test_lattice_wrong_space_input():
    for space in ('', dict(), set(), 3, np.array([]), (1, 0), 2.5):
        with pytest.raises((TypeError, ValueError)) as exc:
            base.Lattice([[1, 0], [0, 1]], space=space)
            if base.checkType(space, 'str'):
                assert exc == ValueError
            else:
                assert exc == TypeError


def test_lattice_wrong_limit_input():
    for space in ('', dict(), set(), np.array([]), (1, 0), 2.5):
        with pytest.raises((TypeError, ValueError)) as exc:
            base.Lattice([[1, 0], [0, 1]], space=space)
            if base.checkType(space, 'str'):
                assert exc == ValueError
            else:
                assert exc == TypeError