import pytest
import shutil, tempfile
import sys
import os
from pathlib import Path
from zipfile import ZipFile
from copy import deepcopy
import numpy as np

from viperleed import lib
import viperleed.lib


from viperleed.lib.classes.r_error import get_n_zero_crossings,get_zero_crossing

@pytest.mark.parametrize('x_arr,n_crossings', ((np.array([-1, 1, -1]),2 ),
                                               (np.array([0, 0, 0]), 1),
                                               ([-1, 0, -1], 1),
                                               ([-1, 1, 10, 2, -50, 0.1, 0], 3))
)
def test_get_n_zero_crossings(x_arr, n_crossings):
    assert get_n_zero_crossings(x_arr) == n_crossings

@pytest.mark.parametrize('x_arr,y_arr,crossing_x', (([0, 1, 2], [-1, 0, 1], 1),
                                                    ([0, 1, 2, 3, 4], [-1, -3, -2, 2, 5], 2.5))
                        )
def test_get_zero_crossing(x_arr, y_arr, crossing_x):
    assert get_zero_crossing(x_arr, y_arr) == crossing_x