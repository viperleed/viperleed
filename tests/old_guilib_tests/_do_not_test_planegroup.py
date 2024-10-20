"""
Test functionality of PlaneGroup class without loading the whole TLEEDMAP

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

def test_planegroup_p2():
    g = base.PlaneGroup('p2')
    assert g.group == 'p2'

    E = [[1, 0], [0, 1]]
    C2 = [[-1, 0], [0, -1]]
    ops = g.operations()
    assert len(ops) == 2
    assert np.array_equal(ops[0], E)
    assert np.array_equal(ops[1], C2)


def test_panegroup_non_string_input_typeerror():
    with pytest.raises(TypeError):
        base.PlaneGroup(9)


def test_planegroup_invalid_group_input_valueerror():
    with pytest.raises(ValueError):
        base.PlaneGroup('a')
