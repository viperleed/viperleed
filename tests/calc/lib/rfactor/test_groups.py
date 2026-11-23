"""Tests for viperleed.calc.lib.rfactor.groups."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)')
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

import numpy as np
import pytest

from viperleed.calc.lib.rfactor.groups import group_rfactors


def test_group_rfactors_overall():
    """Check that groups=None gives the ratio of the summed values."""
    numerators = np.array([1.0, 2.0, 3.0])
    denominators = np.array([4.0, 5.0, 6.0])

    result = group_rfactors(numerators, denominators, groups=None)

    assert result == pytest.approx((1.0 + 2.0 + 3.0) / (4.0 + 5.0 + 6.0))


def test_group_rfactors_per_beam():
    """Check that groups='beam' gives one ratio per beam."""
    numerators = np.array([1.0, 2.0, 3.0])
    denominators = np.array([4.0, 5.0, 6.0])

    result = group_rfactors(numerators, denominators, groups='beam')

    assert result == pytest.approx(np.array([0.25, 0.4, 0.5]))


def test_group_rfactors_integer_groups():
    """Check that integer groups give one ratio of sums per group."""
    numerators = np.array([1.0, 2.0, 3.0, 4.0])
    denominators = np.array([4.0, 5.0, 6.0, 8.0])
    groups = np.array([0, 0, 1, 1])

    result = group_rfactors(numerators, denominators, groups=groups)

    expected = np.array([(1.0 + 2.0) / (4.0 + 5.0), (3.0 + 4.0) / (6.0 + 8.0)])
    assert result == pytest.approx(expected)


def test_group_rfactors_integer_groups_with_num_groups_and_gap():
    """Check that an unused group index gives NaN."""
    numerators = np.array([1.0, 2.0, 3.0])
    denominators = np.array([2.0, 4.0, 6.0])
    groups = np.array([0, 2, 2])

    result = group_rfactors(numerators,
                            denominators,
                            groups=groups,
                            num_groups=3)

    expected = np.array([
        1.0 / 2.0,
        np.nan,
        (2.0 + 3.0) / (4.0 + 6.0),
    ])
    assert result[0] == pytest.approx(expected[0])
    assert np.isnan(result[1])
    assert result[2] == pytest.approx(expected[2])


def test_group_rfactors_single_group_matches_overall():
    """Check that a single group reproduces the overall value."""
    numerators = np.array([1.0, 2.0, 3.0, 4.0])
    denominators = np.array([2.0, 3.0, 4.0, 5.0])
    groups = np.zeros(len(numerators), dtype=int)

    overall = group_rfactors(numerators, denominators, groups=None)
    grouped = group_rfactors(numerators,
                             denominators,
                             groups=groups,
                             num_groups=1)

    assert grouped.shape == (1,)
    assert grouped[0] == pytest.approx(overall)


def test_group_rfactors_unique_groups_match_per_beam():
    """Check that one group per beam reproduces the per-beam values."""
    numerators = np.array([1.0, 2.0, 3.0])
    denominators = np.array([2.0, 4.0, 5.0])
    groups = np.arange(len(numerators))

    grouped = group_rfactors(numerators,
                             denominators,
                             groups=groups,
                             num_groups=len(groups))
    per_beam = group_rfactors(numerators, denominators, groups='beam')

    assert grouped == pytest.approx(per_beam)


def test_group_rfactors_incompatible_shapes():
    """Check that numerators and denominators need the same shape."""
    numerators = np.array([1.0, 2.0])
    denominators = np.array([4.0, 5.0, 6.0])

    with pytest.raises(ValueError, match='same shape'):
        group_rfactors(numerators, denominators, groups=None)
