"""Tests for module realspace of viperleed.gui.leedsim.classes."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-11'
__license__ = 'GPLv3+'

import numpy as np
import pytest

from viperleed.gui.leedsim.classes.realspace import RealSpace

testDict = {'eMax': 15,
            'surfBasis': np.array([[0, 6.76365840], [-7.81, 0]]),
            'SUPERLATTICE': np.array([[1, 2], [-2, 0]]),
            'surfGroup': 'pgg',
            'bulkGroup': 'p6m'}


def test_realspace_init():
    rs = RealSpace(testDict)
    with pytest.raises(AttributeError):                                         # TODO: RealSpace is not yet ready for the multi-domains, so it does not have a .superlattices attribue!
        assert np.array_equal(rs.superlattices[0], testDict['SUPERLATTICE'])
    assert np.array_equal(rs.superlattice, testDict['SUPERLATTICE'])
    assert np.array_equal(rs.surf.basis, testDict['surfBasis'])


def test_realspace_wrong_input_types():
    for key in testDict.keys():
        params = testDict.copy()
        params[key] = dict()
        with pytest.raises((ValueError, TypeError, RuntimeError)):              # TODO: RuntimeError is an awful exception
            RealSpace(params)
    
