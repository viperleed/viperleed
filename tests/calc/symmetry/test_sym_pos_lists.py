"""Tests for the removal of duplicate positions in getSymPosLists."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-27'
__license__ = 'GPLv3+'

from types import SimpleNamespace

import numpy as np
import pytest
from pytest_cases import fixture
from scipy import spatial

from viperleed.calc import symmetry

_SMALL_INDEX = 256  # Python caches ints up to this
_KEPT_IN_CLUSTER = 2  # Origin and the first point of the cluster


def _reference_unique(points, eps):
    """Return points kept by a brute-force greedy de-duplication.

    Each point that is still kept removes all the other points that
    are within eps of it, in list order.
    """
    keep = [True] * len(points)
    for i, point in enumerate(points):
        if not keep[i]:
            continue
        for j, other in enumerate(points):
            if j != i and np.linalg.norm(point - other) <= eps:
                keep[j] = False
    return [p for p, kept in zip(points, keep) if kept]


def _candidates(points):
    """Return the candidate positions that getSymPosLists considers.

    These are the origin, the points, and the midpoints between
    any two of them (for a non-hexagonal cell).
    """
    midpoints = (
        (p + q) / 2
        for i, p in enumerate(points)
        for q in points[i + 1:]
        )
    return [np.zeros(2), *points, *midpoints]


def _reversed_tree(original):
    """Return a subclass of original that finds neighbours backwards."""

    class _ReversedTree(original):  # pylint: disable=too-few-public-methods
        """A tree returning the neighbours in reverse order."""

        def query_ball_point(self, *args, **kwargs):
            """Return the indices of neighbours, highest first."""
            found = super().query_ball_point(*args, **kwargs)
            return sorted(found, reverse=True)

    return _ReversedTree


@fixture(name='neighbours_reversed')
def fixture_neighbours_reversed(mocker):
    """Make the scipy trees return neighbours in descending order.

    scipy does not guarantee that the first index returned for a
    point is the point itself. Both tree classes are replaced, so
    that the test does not depend on the one used in symmetry.
    """
    for name in ('KDTree', 'cKDTree'):
        mocker.patch.object(
            spatial,
            name,
            _reversed_tree(getattr(spatial, name)),
            )


@fixture(name='slab_and_rpars')
def fixture_slab_and_rpars():
    """Return minimal slab and parameters for a square cell."""
    slab = SimpleNamespace(ab_cell=np.eye(2), celltype='square')
    rpars = SimpleNamespace(SYMMETRY_EPS=1e-3)
    return slab, rpars


@pytest.mark.usefixtures('neighbours_reversed')
class TestUniquePositions:
    """Tests for the de-duplication of candidate positions."""

    def test_cluster_keeps_first_point(self, slab_and_rpars):
        """Check that the first of several near-equal points is kept."""
        points = [np.array([0.3, 0.4]), np.array([0.3 + 1e-5, 0.4])]
        # The candidates are the origin, the two points, and their
        # midpoint. The last three are within SYMMETRY_EPS of each
        # other.
        candidates, hex_candidates = symmetry.getSymPosLists(
            *slab_and_rpars,
            points,
            )
        assert not hex_candidates
        assert len(candidates) == _KEPT_IN_CLUSTER
        assert np.allclose(candidates[0], (0, 0))
        assert np.allclose(candidates[1], points[0], atol=1e-9)

    def test_chain_keeps_both_ends(self, slab_and_rpars):
        """Check that a chain of near-equal points keeps both ends."""
        first = np.array([0.3, 0.4])
        last = np.array([first[0] + 16e-4, first[1]])
        # The three points are in a chain: the first is within
        # SYMMETRY_EPS of the second, the second of the third, but
        # the first and the third are not close to each other.
        middle = np.array([first[0] + 8e-4, first[1]])
        points = [first, middle, last]
        candidates, _ = symmetry.getSymPosLists(*slab_and_rpars, points)
        assert np.allclose(candidates, [(0, 0), first, last])

    def test_matches_reference_many_candidates(self):
        """Compare with a brute-force reference for many candidates."""
        eps = 0.1
        slab = SimpleNamespace(ab_cell=np.eye(2) * 3, celltype='square')
        rpars = SimpleNamespace(SYMMETRY_EPS=eps)
        points = list(np.random.default_rng(1).uniform(0, 1, (30, 2)))
        candidates, _ = symmetry.getSymPosLists(slab, rpars, points)
        expected = _reference_unique(_candidates(points), eps)
        assert len(_candidates(points)) > _SMALL_INDEX
        assert len(candidates) == len(expected)
        assert np.allclose(candidates, expected)
