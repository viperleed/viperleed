"""Module direction of viperleed.files.displacements.tokens."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2024-10-14"
__license__ = "GPLv3+"

import re

import numpy as np

from .base import DisplacementsFileToken, TokenParserError

COMPARE_EPS = 1E-6
DIRECTION_PATTERN = r'^(?:(?P<dir>[xyz]+))\[(?P<vec>[\d\s\.\-eE]+)\]$'
SIMPLE_DIRECTIONS = ('x', 'y', 'z')

# Note,TODO: This parser does not yet support parsing of azimuthal, and radial
# directions, as well as fractional directions (e.g., 'a', 'b', 'c'). These
# are supported in th odl viperleed.calc parser.
# They will need to be added here prior to using this parser for the TensErLEED
# backend.
# Comment by @michele-riva:
#Among the ones we currently support, we are missing in here: ab[i j] (and its
# [i j] equivalent), azi(ab[i j]) (and its equivalent azi([i j])), azi(xy[k m]),
# r(ab[i j]) (and its equivalent r([i j])), r(xy[k m]).
# 
# Others that we should probably consider (see also #413): a, b, c, ab, bc, ac
# 
# In addition, we can introduce some syntax leniency and accept also the
# following:
# 
# - various permutations of the letters when displacing in the 2D plane (i.e.,
#   ba == ab, zy == yz, etc)
# - azi[i j] == azi(ab[i j])
# - azi(i j) == azi(ab[i j])
# - r[i j] == r(ab[i j])
# - r(i j) == r(ab[i j])


class DirectionTokenParserError(TokenParserError):
    """Exception raised for unsupported direction formats."""


class DirectionToken(DisplacementsFileToken):
    """Class to parse direction specifiers in 3D space.

    DirectionToken specifiers are used to declare the allowed directions of
    movements or constraints for geometric displacements. They can specify one,
    two, or three dimensions of movements. The direction specifier can be given
    as one of:
    - 'xyz'
        Full 3D space.
    - 'xy', 'yz', 'xz'
        A 2D plane in 3D space.
    - 'x', 'y', 'z', 'xy[a b]', 'xyz[a b c]'
        A 1D vector in 3D space. For the latter two, 'a', 'b', and 'c' are
        placeholders for the components that define the direction of the vector.

    Parameters
    ----------
    direction_str : str
        The direction string to parse.

    Attributes
    ----------
    dof : int
        The number of degrees of freedom (DOF), either 1, 2, or 3.
    vectors : tuple of np.ndarray
        A tuple of one, two or three orthonormal vectors spanning the direction
        space. Note, the vectors are in the zxy convention (LEED convention).

    Raises
    ------
    ValueError
        If the direction string is invalid or if the number of components does
        not match the number of directions.
    DirectionTokenParserError
        If the direction string contains unsupported components like azimuthal
        or radial directions.
    """

    def __init__(self, direction_str):
        _dir_str = direction_str.strip()
        if not _dir_str:
            raise DirectionTokenParserError('Empty direction token.')
        self.direction_str = _dir_str
        vecs, self.dof = self._parse_direction(_dir_str)
        self._vectors = vecs

    @property
    def vectors_xyz(self):
        """Return the vectors in xyz convention."""
        return self._vectors

    @property
    def vectors_zxy(self):
        """Return the vectors in zxy convention (LEED convention)."""
        return _to_zxy(self._vectors)

    def _parse_direction(self, direction_str):
        _check_unsupported_directions(direction_str)

        # Vector direction like 'xyz[1 2 3]'
        if '[' in direction_str:
            match = re.match(DIRECTION_PATTERN, direction_str)
            if not match:
                msg = f'Invalid direction format: {direction_str}'
                raise ValueError(msg)
            dirs = list(match['dir'])
            vecs = list(map(float, match.group('vec').split()))
            if len(dirs) != len(vecs):
                msg = ('Mismatch between directions and components '
                       f'in {direction_str}')
                raise ValueError(msg)
            vec = self._embed_vector(dirs, vecs)
            return (np.array([self._normalize(vec)]), 1)
        # Basis directions like 'x', 'xy', etc.
        if not all(c in SIMPLE_DIRECTIONS for c in direction_str):
            msg = f'Invalid direction: {direction_str}'
            raise ValueError(msg)
        vecs = [self._get_basis_vector(c) for c in direction_str]
        return (np.array(vecs), len(vecs))
   
    def _get_basis_vector(self, label):
        return {
            'x': np.array([1, 0, 0]),
            'y': np.array([0, 1, 0]),
            'z': np.array([0, 0, 1]),
        }[label]

    def _embed_vector(self, dirs, comps):
        vec = np.zeros(3)
        idx_map = {'x': 0, 'y': 1, 'z': 2}
        for d, c in zip(dirs, comps):
            vec[idx_map[d]] = c
        return vec

    def _normalize(self, vec):
        norm = np.linalg.norm(vec)
        if norm < COMPARE_EPS:
            msg = f'Zero-length vector: {vec}'
            raise ValueError(msg)
        return vec / norm

    def __eq__(self, other):
        """Compare two DirectionToken objects for equality."""
        if not isinstance(other, DirectionToken):
            return False
        return self.dof == other.dof and np.allclose(
            self.vectors_xyz, other.vectors_xyz
        )

    def __str__(self):
        """Return a string representation of the DirectionToken object."""
        return f'DirectionToken(vectors={self.vectors_xyz}, dof={self.dof})'


def _check_unsupported_directions(direction_str):
    # Azimuthal & radial directions are currently not supported
    azi_rad_labels = ['azi', 'r']
    if any(label in direction_str for label in azi_rad_labels):
        msg = (
            'Azimuthal and radial directions are currently not supported. '
            f'Invalid direction: {direction_str}'
        )
        raise DirectionTokenParserError(msg)

    # Fractional directions are currently not supported
    fractional_labels = ['a', 'b', 'c']
    if any(label in direction_str for label in fractional_labels):
        msg = (
            'Fractional directions are currently not supported. '
            f'Invalid direction: {direction_str}'
        )
        raise DirectionTokenParserError(msg)

def _to_zxy(vecs):
    """Convert from xyz to zxy axis (LEED) convention."""
    # Apply axis permutation: x→y, y→z, z→x
    return vecs[:, [2, 0, 1]]  # roll last axis backward
