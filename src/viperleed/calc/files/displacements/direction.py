"""Module direction."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-14'

import re

import numpy as np

DIRECTION_PATTERN = r'^(?:(?P<dir>[xyz]+))\[(?P<vec>[\d\s\.\-eE]+)\]$'
SIMPLE_DIRECTIONS = ('x', 'y', 'z')


class UnsupportedDirectionError(Exception):
    """Exception raised for unsupported direction formats."""

    def __init__(self, message):
        super().__init__(message)
        self.message = message

class Direction:
    """Class to parse direction specifiers in 3D space.

    Direction specifiers are used to declare the allowed directions of movements
    or constraints for geometric displacements. They can specify one, two, or
    three dimensions of movements. The direction specifier can be given as one
    of:
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
        space.
    """

    def __init__(self, direction_str):
        self.direction_str = direction_str
        self.vectors, self.dof = self._parse_direction(direction_str)

    def _parse_direction(self, direction_str):
        _check_unsupported_directions(direction_str)

        # Vector direction like 'xyz[1 2 3]'
        if '[' in direction_str:
            match = re.match(DIRECTION_PATTERN, direction_str)
            if not match:
                msg = f'Invalid direction format: {direction_str}'
                raise ValueError(msg)
            dirs = list(match.group('dir'))
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
        if norm == 0:
            msg = f'Zero-length vector: {vec}'
            raise ValueError(msg)
        return vec / norm

    def __eq__(self, other):
        """Compare two Direction objects for equality."""
        if not isinstance(other, Direction):
            return False
        return self.dof == other.dof and np.allclose(
            self.vectors, other.vectors
        )

    def __repr__(self):
        """Return a string representation of the Direction object."""
        return f'Direction(vectors={self.vectors}, dof={self.dof})'


def _check_unsupported_directions(direction_str):
    # Azimuthal & radial directions are currently not supported
    azi_rad_labels = ['azi', 'r']
    if any(label in direction_str for label in azi_rad_labels):
        msg = (
            'Azimuthal and radial directions are currently not supported. '
            f'Invalid direction: {direction_str}'
        )
        raise UnsupportedDirectionError(msg)

    # Fractional directions are currently not supported
    fractional_labels = ['a', 'b', 'c']
    if any(label in direction_str for label in fractional_labels):
        msg = (
            'Fractional directions are currently not supported. '
            f'Invalid direction: {direction_str}'
        )
        raise UnsupportedDirectionError(msg)
