"""Tests for module beam_correspondence of viperleed.calc.classes."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

import math

import numpy as np
import pytest
from pytest_cases import parametrize
from quicktions import Fraction

from viperleed.calc.classes.beam_correspondence import BeamCorrespondence


_HALF = Fraction(1, 2)


_WEIGHTS = {
    'equal': ({(1, 0): 1, (0, 1): 1}, {(1, 0): 0.5, (0, 1): 0.5}),
    'one': ({(1, 0): 3}, {(1, 0): 1.0}),
    'unequal': ({(1, 0): 1, (0, 1): 3}, {(1, 0): 0.25, (0, 1): 0.75}),
    'with zero': ({(1, 0): 0, (0, 1): 0.2}, {(0, 1): 1.0}),
    }


_INVALID = {
    'negative': {(1, 0): {(1, 0): 1, (0, 1): -0.5}},
    'all zero': {(1, 0): {(1, 0): 0, (0, 1): 0}},
    'empty': {(1, 0): {}},
    'infinite': {(1, 0): {(1, 0): math.inf, (0, 1): 1}},
    'not a number': {(1, 0): {(1, 0): math.nan, (0, 1): 1}},
    }


# Keys that differ as Python objects but are the same beam
_COLLIDING = {
    'experimental': {(0.5, 0): {(1, 0): 1}, ('1/2', '0'): {(0, 1): 1}},
    'calculated': {(1, 0): {(0.5, 0): 1, ('1/2', '0'): 1}},
    'calculated, zero weight': {(1, 0): {(0.5, 0): 1, ('1/2', '0'): 0}},
    }


class TestConstruction:
    """Tests for creating a BeamCorrespondence."""

    def test_keys_are_fractions(self):
        """Check that beam indices are stored as exact Fractions."""
        corr = BeamCorrespondence({(1, 0): {(_HALF, -1): 1}})
        (exp_key, calc_weights), = corr.assignments.items()
        calc_key, = calc_weights
        assert exp_key == (Fraction(1), Fraction(0))
        assert calc_key == (Fraction(1, 2), Fraction(-1))
        assert all(isinstance(v, Fraction) for v in exp_key + calc_key)

    def test_equal_regardless_of_input_types(self):
        """Check equality of correspondences built in different ways."""
        from_ints = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 1}})
        from_fractions = BeamCorrespondence({
            (Fraction(1), Fraction(0)): {(Fraction(0), Fraction(1)): 2.5,
                                         (Fraction(1), Fraction(0)): 2.5},
            })
        assert from_ints == from_fractions
        assert from_ints != BeamCorrespondence({(1, 0): {(1, 0): 1}})

    @parametrize('weights,expect', _WEIGHTS.values(), ids=_WEIGHTS)
    def test_weights_normalised(self, weights, expect):
        """Check that weights are normalised per beam assignment."""
        corr = BeamCorrespondence({(1, 0): weights})
        normalised, = corr.assignments.values()
        assert normalised == pytest.approx(
            {(Fraction(h), Fraction(k)): w for (h, k), w in expect.items()}
            )
        assert sum(normalised.values()) == pytest.approx(1)

    @parametrize(assignments=_INVALID.values(), ids=_INVALID)
    def test_invalid_weights(self, assignments):
        """Check complaints for negative or all-zero weights."""
        with pytest.raises(ValueError, match='weight'):
            BeamCorrespondence(assignments)

    @parametrize(assignments=_COLLIDING.values(), ids=_COLLIDING)
    def test_colliding_keys(self, assignments):
        """Check complaints for two keys that are the same beam."""
        with pytest.raises(ValueError, match=r'\(1/2 0\).*more than once'):
            BeamCorrespondence(assignments)

    def test_assignments_is_a_copy(self):
        """Check that .assignments returns a copy."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1}})
        assignments = corr.assignments
        assignments[(1, 0)][(0, 1)] = 1
        assignments[(1, 1)] = {}
        assert corr == BeamCorrespondence({(1, 0): {(1, 0): 1}})


_EQUAL = {
    'unweighted, any scale': (
        {(1, 0): {(1, 0): 1, (0, 1): 1, (1, 1): 1}},
        {(1, 0): {(1, 0): 0.7, (0, 1): 0.7, (1, 1): 0.7}},
        ),
    'unweighted, nine tenths': (
        {(1, 0): {(h, 0): 1 for h in range(9)}},
        {(1, 0): {(h, 0): 0.1 for h in range(9)}},
        ),
    'unweighted, within tolerance': (
        {(1, 0): {(1, 0): 0.3, (0, 1): 0.3, (1, 1): 0.3}},
        {(1, 0): {(1, 0): 1, (0, 1): 1 + 1e-8, (1, 1): 1}},
        ),
    'weighted, normalised': ({(1, 0): {(1, 0): 1, (0, 1): 3}},
                             {(1, 0): {(1, 0): 0.25, (0, 1): 0.75}}),
    'weighted, any scale': ({(1, 0): {(1, 0): 1, (0, 1): 3}},
                            {(1, 0): {(1, 0): 0.1, (0, 1): 0.3}}),
    'weighted, within tolerance': ({(1, 0): {(1, 0): 1, (0, 1): 3}},
                                   {(1, 0): {(1, 0): 1, (0, 1): 3 + 1e-8}}),
    }


_NOT_EQUAL = {
    'different weights': ({(1, 0): {(1, 0): 1, (0, 1): 3}},
                          {(1, 0): {(1, 0): 3, (0, 1): 1}}),
    'weighted, beyond tolerance': ({(1, 0): {(1, 0): 1, (0, 1): 3}},
                                   {(1, 0): {(1, 0): 1, (0, 1): 3.0001}}),
    'weighted and unweighted': ({(1, 0): {(1, 0): 1, (0, 1): 1}},
                                {(1, 0): {(1, 0): 1, (0, 1): 2}}),
    'different calculated beam': ({(1, 0): {(1, 0): 1, (0, 1): 1}},
                                  {(1, 0): {(1, 0): 1, (1, 1): 1}}),
    'different experimental beam': ({(1, 0): {(1, 0): 1}},
                                    {(0, 1): {(1, 0): 1}}),
    'extra assignment': ({(1, 0): {(1, 0): 1}},
                         {(1, 0): {(1, 0): 1}, (1, 1): {(1, 1): 1}}),
    }


# Weights whose normalised values change by
# rounding when they are normalised again
_ROUND_TRIP = {
    'small weight': {(0, 0): 0.762280082457942,
                     (1, 0): 0.0021060533511106927,
                     (2, 0): 0.4453871940548014},
    'similar weights': {(0, 0): 0.1859062658947177,
                        (1, 0): 0.9925434121760651,
                        (2, 0): 0.8599465287952899},
    }


class TestEquality:
    """Tests for comparing BeamCorrespondence objects."""

    @parametrize('first,second', _EQUAL.values(), ids=_EQUAL)
    def test_equal(self, first, second):
        """Check correspondences that compare equal."""
        assert BeamCorrespondence(first) == BeamCorrespondence(second)

    @parametrize('first,second', _NOT_EQUAL.values(), ids=_NOT_EQUAL)
    def test_not_equal(self, first, second):
        """Check correspondences that compare unequal."""
        assert BeamCorrespondence(first) != BeamCorrespondence(second)

    @parametrize(weights=_ROUND_TRIP.values(), ids=_ROUND_TRIP)
    def test_round_trip(self, weights):
        """Check that a correspondence equals one made from its copy."""
        corr = BeamCorrespondence({(1, 0): weights})
        assert corr.is_weighted
        assert BeamCorrespondence(corr.assignments) == corr

    def test_not_equal_other_type(self):
        """Check that other objects never compare equal."""
        assignments = {(1, 0): {(1, 0): 1}}
        corr = BeamCorrespondence(assignments)
        assert corr != assignments
        assert (corr == assignments) is False


_WITH_ZERO = {(1, 0): {(1, 0): 1, (0, 1): 0}, (0, 1): {(0, 1): 1}}


class TestZeroWeights:
    """Tests for calculated beams with zero weight."""

    def test_dropped(self):
        """Check that zero-weight calculated beams are left out."""
        without_zero = {(1, 0): {(1, 0): 1}, (0, 1): {(0, 1): 1}}
        assert (BeamCorrespondence(_WITH_ZERO)
                == BeamCorrespondence(without_zero))

    def test_not_weighted(self):
        """Check that a zero weight does not make a weighted one."""
        assert not BeamCorrespondence(_WITH_ZERO).is_weighted

    def test_not_shared(self):
        """Check that a zero-weight beam is not a shared beam."""
        assert not BeamCorrespondence(_WITH_ZERO).has_shared_beams


_SHARED = {(1, 0): {(1, 0): 1},
           (1, 1): {(1, 1): 1},
           (0, 1): {(0, 1): 1, (_HALF, 0): 1},
           (2, 0): {(2, 0): 1, (_HALF, 0): 1}}


_IS_WEIGHTED = (
    pytest.param({(1, 0): {(1, 0): 1}, (1, 1): {(1, 1): 7}}, False,
                 id='single beams'),
    pytest.param({(1, 0): {(1, 0): 2, (0, 1): 2, (-1, 0): 2}}, False,
                 id='equal'),
    pytest.param({(1, 0): {(1, 0): 1, (0, 1): 1},
                  (1, 1): {(1, 1): 3, (-1, 1): 3}}, False,
                 id='equal per assignment'),
    pytest.param({(1, 0): {(1, 0): 1, (0, 1): 1 + 1e-8}}, False,
                 id='within tolerance'),
    pytest.param(_SHARED, False, id='shared beams'),
    pytest.param({(1, 0): {(1, 0): 1, (0, 1): 2}}, True, id='unequal'),
    pytest.param({(1, 0): {(1, 0): 1},
                  (1, 1): {(1, 1): 0.4, (-1, 1): 0.6}}, True,
                 id='one of two'),
    pytest.param({(1, 0): {(1, 0): 1, (0, 1): 1 + 1e-5}}, True,
                 id='beyond tolerance'),
    )


_HAS_SHARED_BEAMS = (
    pytest.param({(1, 0): {(1, 0): 1, (0, 1): 1},
                  (1, 1): {(1, 1): 1, (-1, 1): 1},
                  (2, 0): {(2, 0): 1}}, False,
                 id='not shared'),
    pytest.param(_SHARED, True, id='shared'),
    )


class TestProperties:
    """Tests for is_weighted and has_shared_beams."""

    def test_is_weighted_order_independent(self):
        """Check that key order does not decide if weights are equal."""
        weights = {(1, 0): 1, (0, 1): 1 + 0.9e-6, (1, 1): 1 - 0.9e-6}
        reordered = {key: weights[key] for key in ((0, 1), (1, 0), (1, 1))}
        corr = BeamCorrespondence({(1, 0): weights})
        corr_reordered = BeamCorrespondence({(1, 0): reordered})
        assert corr.is_weighted is corr_reordered.is_weighted
        assert corr == corr_reordered

    @parametrize('assignments,expect', _IS_WEIGHTED)
    def test_is_weighted(self, assignments, expect):
        """Check whether a correspondence counts as weighted."""
        assert BeamCorrespondence(assignments).is_weighted is expect

    @parametrize('assignments,expect', _HAS_SHARED_BEAMS)
    def test_has_shared_beams(self, assignments, expect):
        """Check whether a calculated beam is in many assignments."""
        corr = BeamCorrespondence(assignments)
        assert corr.has_shared_beams is expect


_NO_LEGACY = {
    'weighted': {(1, 0): {(1, 0): 1, (0, 1): 2}},
    'shared': {(1, 0): {(1, 0): 1, (0, 1): 1}, (0, 1): {(0, 1): 1}},
    }


class TestLegacyVector:
    """Tests for the per-calculated-beam vector."""

    def test_vector(self, make_beams):
        """Check the experimental-beam index of each calculated beam."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 1},
                                   (1, 1): {(-1, 1): 1},
                                   (_HALF, 0): {(_HALF, 0): 1}})
        expbeams = make_beams((1, 1), (1, 0), (0.5, 0))
        ivbeams = make_beams((0, 1), (2, 0), (1, 0), (-1, 1), (0.5, 0))
        assert corr.legacy_vector(expbeams, ivbeams) == [1, -1, 1, 0, 2]

    def test_first_and_last_index(self, make_beams):
        """Check the first and last experimental beams."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1},
                                   (2, 0): {(2, 0): 1}})
        expbeams = make_beams((1, 0), (1, 1), (2, 0))
        ivbeams = make_beams((2, 0), (1, 0))
        assert corr.legacy_vector(expbeams, ivbeams) == [2, 0]

    def test_tolerance(self, make_beams):
        """Check that beams match within the tolerance of 1e-5."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 1}})
        expbeams = make_beams((1 + 4e-6, -4e-6))
        ivbeams = make_beams((1, 0), (4e-6, 1 + 4e-6), (2e-5, 1))
        assert corr.legacy_vector(expbeams, ivbeams) == [0, 0, -1]

    def test_experimental_beam_missing(self, make_beams):
        """Check calculated beams of an absent experimental beam."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1},
                                   (1, 1): {(1, 1): 1}})
        expbeams = make_beams((1, 1))
        ivbeams = make_beams((1, 0), (1, 1))
        assert corr.legacy_vector(expbeams, ivbeams) == [-1, 0]

    @parametrize(assignments=_NO_LEGACY.values(), ids=_NO_LEGACY)
    def test_raises(self, assignments, make_beams):
        """Check that weighted or shared correspondences raise."""
        corr = BeamCorrespondence(assignments)
        expbeams = make_beams((1, 0), (0, 1))
        ivbeams = make_beams((1, 0), (0, 1))
        with pytest.raises(ValueError, match='legacy'):
            corr.legacy_vector(expbeams, ivbeams)


class TestWeightMatrix:
    """Tests for the experimental x calculated weight matrix."""

    def test_matrix(self, make_beams):
        """Check the matrix of a weighted, shared correspondence."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 3},
                                   (1, 1): {(1, 1): 1, (0, 1): 1}})
        expbeams = make_beams((1, 1), (2, 0), (1, 0))
        ivbeams = make_beams((0, 1), (1, 0), (1, 1), (3, 0))
        expect = np.array([[0.5, 0, 0.5, 0],
                           [0, 0, 0, 0],
                           [0.75, 0.25, 0, 0]])
        matrix = corr.weight_matrix(expbeams, ivbeams)
        assert matrix.shape == (3, 4)
        np.testing.assert_allclose(matrix, expect)
        np.testing.assert_allclose(matrix.sum(axis=1), [1, 0, 1])

    def test_experimental_beam_missing(self, make_beams):
        """Check that absent experimental beams are skipped."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1},
                                   (2, 0): {(2, 0): 1}})
        matrix = corr.weight_matrix(make_beams((1, 0)), make_beams((1, 0)))
        np.testing.assert_allclose(matrix, [[1]])

    def test_calculated_beam_missing(self, make_beams):
        """Check complaints when a calculated beam is missing."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 1}})
        with pytest.raises(ValueError, match=r'\(0 1\)'):
            corr.weight_matrix(make_beams((1, 0)), make_beams((1, 0)))
