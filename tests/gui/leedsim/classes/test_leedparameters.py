"""Tests for module leedparameters of viperleed.gui.leedsim.classes."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-10'
__license__ = 'GPLv3+'

import numpy as np
from pytest_cases import parametrize

from viperleed.gui.leedsim.classes.leedparameters import LEEDParametersList

_BASIS = ((2.0, 0.0), (-0.520944533, 2.954423259))
_MATRICES = {
    'tuple of floats': _BASIS,
    'tuple of numpy floats': tuple(tuple(np.float64(v) for v in row)
                                   for row in _BASIS),
    'array': np.array(_BASIS),
    }


@parametrize('surf_basis', _MATRICES.values(), ids=_MATRICES)
def test_matrix_input(surf_basis):
    """Check that lattice matrices are read from various types.

    The str() of numpy scalars (e.g., 'np.float64(2.0)' since numpy
    2.0) cannot be parsed back. Hence, sequences of numpy scalars must
    not be converted to strings with str().
    """
    params = LEEDParametersList({'SUPERLATTICE': ((1, 0), (0, 1)),
                                 'surfBasis': surf_basis,
                                 'surfGroup': 'p1', 'bulkGroup': 'p1',
                                 'eMax': 700})
    assert np.allclose(params[0]['surfBasis'], _BASIS)
    assert np.array_equal(params[0]['SUPERLATTICE'], ((1, 0), (0, 1)))
