"""Tests for module beam_correspondence of viperleed.calc.classes."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

import numpy as np
import pytest
from pytest_cases import parametrize
from quicktions import Fraction

from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.beam_correspondence import BeamCorrespondence


_HALF = Fraction(1, 2)


def _beams(*indices):
    """Return a list of Beam objects with the given (h, k) indices."""
    return [Beam(hk) for hk in indices]


_WEIGHTS = {
    'equal': ({(1, 0): 1, (0, 1): 1}, {(1, 0): 0.5, (0, 1): 0.5}),
    'one': ({(1, 0): 3}, {(1, 0): 1.0}),
    'unequal': ({(1, 0): 1, (0, 1): 3}, {(1, 0): 0.25, (0, 1): 0.75}),
    'with zero': ({(1, 0): 0, (0, 1): 0.2},
                  {(1, 0): 0.0, (0, 1): 1.0}),
    }


_INVALID = {
    'negative': {(1, 0): {(1, 0): 1, (0, 1): -0.5}},
    'all zero': {(1, 0): {(1, 0): 0, (0, 1): 0}},
    'empty': {(1, 0): {}},
    }


_NOT_MAPPINGS = {
    'top level': [((1, 0), {(1, 0): 1})],
    'assignment': {(1, 0): [(1, 0)]},
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

    def test_len(self):
        """Check the number of beam assignments."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1}, (1, 1): {(1, 1): 1}})
        assert len(corr) == 2


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


    @parametrize(assignments=_NOT_MAPPINGS.values(), ids=_NOT_MAPPINGS)
    def test_not_a_mapping(self, assignments):
        """Check complaints for input that is not a Mapping."""
        with pytest.raises(TypeError):
            BeamCorrespondence(assignments)

    def test_assignments_is_a_copy(self):
        """Check that .assignments returns a copy."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1}})
        assignments = corr.assignments
        assignments[(1, 0)][(0, 1)] = 1
        assignments[(1, 1)] = {}
        assert corr == BeamCorrespondence({(1, 0): {(1, 0): 1}})


_UNWEIGHTED = {
    'single beams': {(1, 0): {(1, 0): 1}, (1, 1): {(1, 1): 7}},
    'equal': {(1, 0): {(1, 0): 2, (0, 1): 2, (-1, 0): 2}},
    'equal per assignment': {(1, 0): {(1, 0): 1, (0, 1): 1},
                             (1, 1): {(1, 1): 3, (-1, 1): 3}},
    'within tolerance': {(1, 0): {(1, 0): 1, (0, 1): 1 + 1e-8}},
    }


_WEIGHTED = {
    'unequal': {(1, 0): {(1, 0): 1, (0, 1): 2}},
    'one of two': {(1, 0): {(1, 0): 1},
                   (1, 1): {(1, 1): 0.4, (-1, 1): 0.6}},
    'beyond tolerance': {(1, 0): {(1, 0): 1, (0, 1): 1 + 1e-5}},
    'zero weight': {(1, 0): {(1, 0): 1, (0, 1): 0}},
    }


class TestProperties:
    """Tests for is_weighted and has_shared_beams."""


    @parametrize(assignments=_UNWEIGHTED.values(), ids=_UNWEIGHTED)
    def test_unweighted(self, assignments):
        """Check that equal weights count as unweighted."""
        assert not BeamCorrespondence(assignments).is_weighted

    @parametrize(assignments=_WEIGHTED.values(), ids=_WEIGHTED)
    def test_weighted(self, assignments):
        """Check that unequal weights count as weighted."""
        assert BeamCorrespondence(assignments).is_weighted

    def test_not_shared(self):
        """Check a correspondence without shared calculated beams."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 1},
                                   (1, 1): {(1, 1): 1, (-1, 1): 1},
                                   (2, 0): {(2, 0): 1}})
        assert not corr.has_shared_beams

    def test_shared(self):
        """Check a calculated beam in two beam assignments."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1},
                                   (1, 1): {(1, 1): 1},
                                   (0, 1): {(0, 1): 1, (_HALF, 0): 1},
                                   (2, 0): {(2, 0): 1, (_HALF, 0): 1}})
        assert corr.has_shared_beams
        assert not corr.is_weighted


_NO_LEGACY = {
    'weighted': {(1, 0): {(1, 0): 1, (0, 1): 2}},
    'shared': {(1, 0): {(1, 0): 1, (0, 1): 1}, (0, 1): {(0, 1): 1}},
    }


class TestLegacyVector:
    """Tests for the per-calculated-beam vector."""

    def test_vector(self):
        """Check the experimental-beam index of each calculated beam."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 1},
                                   (1, 1): {(-1, 1): 1},
                                   (_HALF, 0): {(_HALF, 0): 1}})
        expbeams = _beams((1, 1), (1, 0), (0.5, 0))
        ivbeams = _beams((0, 1), (2, 0), (1, 0), (-1, 1), (0.5, 0))
        assert corr.legacy_vector(expbeams, ivbeams) == [1, -1, 1, 0, 2]

    def test_first_and_last_index(self):
        """Check the first and last experimental beams."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1},
                                   (2, 0): {(2, 0): 1}})
        expbeams = _beams((1, 0), (1, 1), (2, 0))
        ivbeams = _beams((2, 0), (1, 0))
        assert corr.legacy_vector(expbeams, ivbeams) == [2, 0]

    def test_tolerance(self):
        """Check that beams match within the tolerance of 1e-5."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 1}})
        expbeams = _beams((1 + 4e-6, -4e-6))
        ivbeams = _beams((1, 0), (4e-6, 1 + 4e-6), (2e-5, 1))
        assert corr.legacy_vector(expbeams, ivbeams) == [0, 0, -1]

    def test_experimental_beam_missing(self):
        """Check calculated beams of an absent experimental beam."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1},
                                   (1, 1): {(1, 1): 1}})
        expbeams = _beams((1, 1))
        ivbeams = _beams((1, 0), (1, 1))
        assert corr.legacy_vector(expbeams, ivbeams) == [-1, 0]


    @parametrize(assignments=_NO_LEGACY.values(), ids=_NO_LEGACY)
    def test_raises(self, assignments):
        """Check that weighted or shared correspondences raise."""
        corr = BeamCorrespondence(assignments)
        expbeams = _beams((1, 0), (0, 1))
        ivbeams = _beams((1, 0), (0, 1))
        with pytest.raises(ValueError, match='legacy'):
            corr.legacy_vector(expbeams, ivbeams)


class TestWeightMatrix:
    """Tests for the experimental x calculated weight matrix."""

    def test_matrix(self):
        """Check the matrix of a weighted, shared correspondence."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 3},
                                   (1, 1): {(1, 1): 1, (0, 1): 1}})
        expbeams = _beams((1, 1), (2, 0), (1, 0))
        ivbeams = _beams((0, 1), (1, 0), (1, 1), (3, 0))
        expect = np.array([[0.5, 0, 0.5, 0],
                           [0, 0, 0, 0],
                           [0.75, 0.25, 0, 0]])
        matrix = corr.weight_matrix(expbeams, ivbeams)
        assert matrix.shape == (3, 4)
        np.testing.assert_allclose(matrix, expect)
        np.testing.assert_allclose(matrix.sum(axis=1), [1, 0, 1])

    def test_experimental_beam_missing(self):
        """Check that absent experimental beams are skipped."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1},
                                   (2, 0): {(2, 0): 1}})
        matrix = corr.weight_matrix(_beams((1, 0)), _beams((1, 0)))
        np.testing.assert_allclose(matrix, [[1]])

    def test_calculated_beam_missing(self):
        """Check complaints when a calculated beam is missing."""
        corr = BeamCorrespondence({(1, 0): {(1, 0): 1, (0, 1): 1}})
        with pytest.raises(ValueError, match=r'\(0 1\)'):
            corr.weight_matrix(_beams((1, 0)), _beams((1, 0)))
