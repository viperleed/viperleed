"""
Test functionality of the RealSpace class without loading the whole TLEEDMAP

Created: 2020-01-11
Author: Michele Riva
"""

import sys
import os
from copy import copy

tests_path = os.path.realpath(os.path.dirname(__file__))
base_path = os.path.realpath(os.path.join(tests_path, '..'))
for path in [tests_path, base_path]:
    if path not in sys.path:
        sys.path.append(path)

import pytest
import numpy as np

import classes

testDict = {'eMax': 15,
            'surfBasis': np.array([[0, 6.76365840], [-7.81, 0]]),
            'SUPERLATTICE': np.array([[1, 2], [-2, 0]]),
            'surfGroup': 'pgg',
            'bulkGroup': 'p6m'}


def test_realspace_init():
    rs = classes.RealSpace(testDict)
    assert np.array_equal(rs.superlattices[0], testDict['SUPERLATTICE'])
    assert np.array_equal(rs.surf.basis, testDict['surfBasis'])


def test_realspace_wrong_input_types():
    for key in testDict.keys():
        params = copy(testDict)
        params[key] = dict()
        with pytest.raises(TypeError):
            classes.RealSpace(params)
    
