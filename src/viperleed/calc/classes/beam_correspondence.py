"""Module beam_correspondence of viperleed.calc.classes.

Defines the BeamCorrespondence class, which holds the assignment of
calculated beams to experimental beams used in the R factor.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from collections.abc import Mapping
import math

import numpy as np
from quicktions import Fraction


# Tolerance for matching beam indices (h, k) to Beam objects
BEAM_EPS = 1e-5

# Relative tolerance below which the assignment
# weights of one beam assignment count as equal
WEIGHT_REL_TOL = 1e-6


class BeamCorrespondence:
    """The assignment of calculated beams to experimental beams.

    Each beam assignment maps one experimental beam to the calculated
    beams whose weighted average is compared with it in the R factor.
    Beams are identified by their (h, k) indices as exact fractions,
    so that equal correspondences compare equal regardless of how
    they were created. Assignment weights are normalised so that they
    sum to one in each beam assignment.

    Attributes
    ----------
    assignments : dict
        {experimental_hk: {calculated_hk: weight}}, with each hk a
        tuple of two Fractions. A copy: changing it does not change
        this BeamCorrespondence.
    is_weighted : bool
        Whether some beam assignment has unequal weights.
    has_shared_beams : bool
        Whether some calculated beam appears in more
        than one beam assignment.
    """

    def __init__(self, assignments):
        """Initialise from beam assignments.

        Parameters
        ----------
        assignments : Mapping
            {experimental_hk: {calculated_hk: weight}}. Each hk is a
            sequence of two numbers convertible to Fraction (floats
            are converted exactly: pass Fractions or integers, e.g.,
            a Beam's `hkfrac`). Weights are non-negative numbers.
            Calculated beams with zero weight are left out; the other
            weights are normalised to sum to one for each experimental
            beam. Weights that are all equal within a relative
            tolerance of WEIGHT_REL_TOL are stored as exactly 1/n
            for n calculated beams.

        Raises
        ------
        TypeError
            If `assignments`, or one of its values, is not a Mapping.
        ValueError
            If a weight is negative, or if a beam assignment
            has no calculated beam or only zero weights.
        """
        if not isinstance(assignments, Mapping):
            msg = ('assignments must be a Mapping, not '
                   f'{type(assignments).__name__!r}')
            raise TypeError(msg)
        self._assignments = {}
        for exp_hk, calc_weights in assignments.items():
            exp_key = _to_key(exp_hk)
            self._assignments[exp_key] = _normalised(exp_key, calc_weights)

    def __eq__(self, other):
        """Return whether self and other hold the same assignments."""
        if not isinstance(other, BeamCorrespondence):
            return NotImplemented
        return self._assignments == other._assignments

    __hash__ = None  # Mutable-looking container, not hashable

    def __len__(self):
        """Return the number of beam assignments."""
        return len(self._assignments)

    def __repr__(self):
        """Return a string representation of self."""
        return f'{type(self).__name__}({self._assignments!r})'

    @property
    def assignments(self):
        """Return a copy of the beam assignments."""
        return {exp_key: dict(calc_weights)
                for exp_key, calc_weights in self._assignments.items()}

    @property
    def has_shared_beams(self):
        """Return whether a calculated beam is in many assignments."""
        seen = set()
        for calc_weights in self._assignments.values():
            if seen.intersection(calc_weights):
                return True
            seen.update(calc_weights)
        return False

    @property
    def is_weighted(self):
        """Return whether some beam assignment has unequal weights."""
        return not all(_all_equal(calc_weights.values())
                       for calc_weights in self._assignments.values())

    def legacy_vector(self, expbeams, ivbeams):
        """Return the experimental-beam index for each calculated beam.

        This is the format used by TensErLEED (KAV): one entry per
        calculated beam, holding the index of its experimental beam,
        or -1 for calculated beams without experimental beam.

        Parameters
        ----------
        expbeams : Sequence of Beam
            The experimental beams.
        ivbeams : Sequence of Beam
            The calculated beams.

        Returns
        -------
        list of int
            One entry per beam in `ivbeams`: the index in `expbeams`
            of the experimental beam it is assigned to, or -1 if the
            calculated beam is in no beam assignment, or if its
            experimental beam is not in `expbeams`.

        Raises
        ------
        ValueError
            If the correspondence is weighted or a calculated beam
            appears in more than one assignment: these cannot be
            expressed as one experimental beam per calculated beam.
        """
        if self.is_weighted:
            raise ValueError('A weighted beam correspondence has no '
                             'legacy (per-calculated-beam) form.')
        if self.has_shared_beams:
            raise ValueError('A beam correspondence with calculated beams '
                             'in more than one beam assignment has no '
                             'legacy (per-calculated-beam) form.')
        calc_to_exp_index = {}
        for exp_key, calc_weights in self._assignments.items():
            exp_index = _find_index(exp_key, expbeams)
            for calc_key in calc_weights:
                calc_to_exp_index[calc_key] = exp_index
        vector = []
        for beam in ivbeams:
            calc_key = _find_key(beam, calc_to_exp_index)
            vector.append(calc_to_exp_index.get(calc_key, -1))
        return vector

    def weight_matrix(self, expbeams, ivbeams):
        """Return the assignment weights as an array.

        Parameters
        ----------
        expbeams : Sequence of Beam
            The experimental beams.
        ivbeams : Sequence of Beam
            The calculated beams.

        Returns
        -------
        numpy.ndarray
            Shape (len(expbeams), len(ivbeams)). Element [i, j] is the
            weight of calculated beam ivbeams[j] in the assignment of
            experimental beam expbeams[i]. Each row sums to one, except
            for experimental beams in no beam assignment: their rows
            are all zero.

        Raises
        ------
        ValueError
            If a calculated beam of an assignment whose experimental
            beam is in `expbeams` is not in `ivbeams`.
        """
        matrix = np.zeros((len(expbeams), len(ivbeams)))
        for exp_key, calc_weights in self._assignments.items():
            exp_index = _find_index(exp_key, expbeams)
            if exp_index < 0:
                continue
            for calc_key, weight in calc_weights.items():
                calc_index = _find_index(calc_key, ivbeams)
                if calc_index < 0:
                    msg = (f'Calculated beam {_format_key(calc_key)}, '
                           'assigned to experimental beam '
                           f'{_format_key(exp_key)}, is not among the '
                           'calculated beams.')
                    raise ValueError(msg)
                matrix[exp_index, calc_index] = weight
        return matrix


def _all_equal(weights):
    """Return whether all `weights` are equal within WEIGHT_REL_TOL."""
    first, *others = weights
    return all(math.isclose(w, first, rel_tol=WEIGHT_REL_TOL)
               for w in others)


def _find_index(key, beams):
    """Return the index of the first beam matching `key`, or -1."""
    for i, beam in enumerate(beams):
        if beam.isEqual_hk(key, eps=BEAM_EPS):
            return i
    return -1


def _find_key(beam, keys):
    """Return the first of `keys` matching `beam`, or None."""
    for key in keys:
        if beam.isEqual_hk(key, eps=BEAM_EPS):
            return key
    return None


def _format_key(key):
    """Return a string representation of a beam key."""
    first, second = key
    return f'({first} {second})'


def _normalised(exp_key, calc_weights):
    """Return {calc_key: normalised weight} for one beam assignment."""
    if not isinstance(calc_weights, Mapping):
        msg = ('Each beam assignment must be a Mapping, not '
               f'{type(calc_weights).__name__!r}')
        raise TypeError(msg)
    weights = {}
    for calc_hk, weight in calc_weights.items():
        if weight < 0:
            msg = (f'Negative weight {weight} in the assignment of '
                   f'experimental beam {_format_key(exp_key)}.')
            raise ValueError(msg)
        if weight:  # Zero-weight beams are in no assignment
            weights[_to_key(calc_hk)] = float(weight)
    total = sum(weights.values())
    if not total > 0:  # Also no calculated beams, or NaN
        msg = (f'The assignment of experimental beam {_format_key(exp_key)} '
               'has no calculated beam with a non-zero weight.')
        raise ValueError(msg)
    if _all_equal(weights.values()):
        # Unweighted: store exactly 1/n, so that the result
        # does not depend on the scale of the input weights
        return {calc_key: 1 / len(weights) for calc_key in weights}
    return {calc_key: weight / total for calc_key, weight in weights.items()}


def _to_key(indices):
    """Return a hashable (h, k) tuple of Fractions."""
    first, second = indices
    return Fraction(first), Fraction(second)
