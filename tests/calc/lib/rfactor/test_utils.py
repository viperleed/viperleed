"""Tests for viperleed.calc.lib.rfactor.utils."""

import numpy as np
import pytest

from viperleed.calc.lib.rfactor.utils import average_beam_array


def test_average_beam_array_permutation_only():
    beam_array = np.array([
        [1.0, 2.0, 3.0],
        [10.0, 20.0, 30.0],
    ])
    beam_corr = (0, 2, 1)

    averaged = average_beam_array(beam_array, beam_corr)

    expected = np.array([
        [1.0, 3.0, 2.0],
        [10.0, 30.0, 20.0],
    ])
    assert averaged == pytest.approx(expected)


def test_average_beam_array_averages_duplicate_targets():
    beam_array = np.array([
        [1.0, 3.0, 10.0],
        [2.0, 5.0, 20.0],
    ])
    beam_corr = (0, 0, 1)

    averaged = average_beam_array(beam_array, beam_corr)

    expected = np.array([
        [2.0, 10.0],
        [3.5, 20.0],
    ])
    assert averaged == pytest.approx(expected)


def test_average_beam_array_filters_minus_one():
    beam_array = np.array([
        [1.0, 3.0, 10.0],
        [2.0, 5.0, 20.0],
    ])
    beam_corr = (0, -1, 1)

    averaged = average_beam_array(beam_array, beam_corr)

    expected = np.array([
        [1.0, 10.0],
        [2.0, 20.0],
    ])
    assert averaged == pytest.approx(expected)


def test_average_beam_array_duplicate_and_minus_one_together():
    beam_array = np.array([
        [1.0, 3.0, 10.0, 100.0],
        [2.0, 5.0, 20.0, 200.0],
    ])
    beam_corr = (0, 0, -1, 1)

    averaged = average_beam_array(beam_array, beam_corr)

    expected = np.array([
        [2.0, 100.0],
        [3.5, 200.0],
    ])
    assert averaged == pytest.approx(expected)
