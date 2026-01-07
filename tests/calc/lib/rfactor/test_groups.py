"""Tests for viperleed.calc.lib.rfactor.groups."""

import numpy as np
import pytest

from viperleed.calc.lib.rfactor.groups import group_rfactors


class TestGroupedRfactors:
    def test_per_beam_true(self):
        numerators = np.array([1, 2, 3])
        denominators = np.array([4, 5, 6])
        result = group_rfactors(numerators, denominators, groups='beam')
        assert result == pytest.approx(np.array([0.25, 0.4, 0.5]))

    def test_per_beam_false(self):
        numerators = np.array([1, 2, 3])
        denominators = np.array([4, 5, 6])
        result = group_rfactors(numerators, denominators, groups=None)
        assert result == pytest.approx(np.array([0.4]))

    def test_per_beam_index_array(self):
        numerators = np.array([1, 2, 3, 4])
        denominators = np.array([4, 5, 6, 8])
        group_indices = np.array([0, 0, 1, 1])
        result = group_rfactors(
            numerators,
            denominators,
            groups=group_indices,
        )
        expected = np.array([1 / 3, 0.5])
        assert result == pytest.approx(expected)

    def test_incompatible_shapes(self):
        numerators = np.array([1, 2])
        denominators = np.array([4, 5, 6])
        with pytest.raises(ValueError):
            group_rfactors(
                numerators,
                denominators,
                groups=None,
            )
