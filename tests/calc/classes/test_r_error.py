"""Tests for viperleed.calc.classes.RError."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-03-09'
__license__ = 'GPLv3+'

import numpy as np
import pytest

from viperleed.calc.classes import r_error


class TestZeroCrossings:
    """Tests for methods detecting zeros."""

    _x_only = (  # x_values, nr_expected_crossings
        (np.array([-1, 1, -1]),       2),
        (np.array([0, 0, 0]),         1),
        ([-1, 0, -1],                 1),
        ([-1, 1, 10, 2, -50, 0.1, 0], 3),
        )

    @pytest.mark.parametrize('x_vals,n_crossings', _x_only)
    def test_nr_zero_crossings(self, x_vals, n_crossings):
        """Check that the number of detected crossings is correct."""
        assert r_error.get_n_zero_crossings(x_vals) == n_crossings

    _x_and_y = (  # (x_values, y_values), cross_x_position
        (([0, 1, 2],       [-1, 0, 1]),         1),
        (([0, 1, 2, 3, 4], [-1, -3, -2, 2, 5]), 2.5),
        (([-2, -1, 0, 1, 2], [-5, -4, -3, 1, 3]), 0.75),
        (([-2, -1, 0, 1, 2], [5, 4, 3, -1, -3]), 0.75),
        )

    @pytest.mark.parametrize('xy_vals,crossing_x', _x_and_y)
    def test_pos_zero_crossing(self, xy_vals, crossing_x):
        """Check that the x position of the detected crossing is correct."""
        assert r_error.get_zero_crossing(*xy_vals) == crossing_x
