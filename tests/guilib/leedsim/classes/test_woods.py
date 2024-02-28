"""Test functionality of module woods of viperleed.guilib.leedsim.classes.

Created: 2024-02-20
Author: Michele Riva (@michele-riva)
"""

import numpy as np
import pytest
from pytest_cases import parametrize

from viperleed.guilib.leedsim.classes import woods


_commensurate = {
    'integer floats': ([[1.0, 2.0], [3.0, 4.0]], True),
    'integers': ([[1, 2], [3, 4]], True),
    'singular': ([[1, 2], [2, 4]], False),
    'non int floats': ([[1.0, 2.0], [3.0, 4.01]], False),
    }

@parametrize('matrix,expect', _commensurate.values(), ids=_commensurate)
def test_is_commensurate(matrix, expect):
    """Check correct result of is_commensurate."""
    commensurate = woods.is_commensurate(matrix)
    assert commensurate is expect

