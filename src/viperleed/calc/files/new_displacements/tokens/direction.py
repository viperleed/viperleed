"""Module direction of viperleed.calc.files.new_displacements.tokens."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

from abc import abstractmethod
import re

import numpy as np

from viperleed.calc.constants import DISPLACEMENTS_FILE_EPS

from .base import DisplacementsFileToken, TokenParserError

CART_DIRECTION_PATTERN = r'^(?:(?P<dir>[xyz]+))\[(?P<vec>[\d\s\.\-eE]+)\]$'
FRACTIONAL_DIRECTION_PATTERN = (
    r'^(?P<dir>ab)?\[(?P<vec>[\d\s\.\-eE]+)\]$'
)
SIMPLE_DIRECTIONS = ('x', 'y', 'z')
FRACTIONAL_DIRECTIONS = ('a', 'b')
FRACTIONAL_PLANE = 'ab'
UNSUPPORTED_C_DIRECTION = 'c'


class DirectionTokenParserError(TokenParserError):
    """Exception raised for unsupported direction formats."""


class DirectionToken(DisplacementsFileToken):
    """Abstract base class for direction tokens in DISPLACEMENTS.

    This class defines the interface for direction tokens, which are
    used to specify the allowed directions of movements or constraints
    for geometric displacements.

    Supported are Cartesian directions (CartesianDirectionToken: 'x',
    'xy', 'xyz', 'xy[i j]', 'xyz[i j k]', ...) and fractional ones,
    in terms of the in-plane unit-cell vectors
    (FractionalDirectionToken: 'a', 'b', 'ab', 'ab[i j]' and its
    equivalent '[i j]'). Use make_direction_token to get the right
    one for a string.

    Not yet supported are the azimuthal and radial directions of the
    old viperleed.calc parser: azi(ab[i j]) (and its equivalent
    azi([i j])), azi(xy[k m]), r(ab[i j]) (and its equivalent
    r([i j])), r(xy[k m]). They are needed before this parser can be
    used for the TensErLEED backend. Directions along the c unit-cell
    vector ('c', 'ac', 'bc', 'abc') and directions mixing the
    Cartesian and fractional frames (e.g., 'az') are rejected.

    Comment by @michele-riva: we can introduce some syntax leniency and
    accept also the following:

    - various permutations of the letters when displacing in the 2D
      plane (i.e. ba == ab, zy == yz, etc):
    - azi[i j] == azi(ab[i j])
    - azi(i j) == azi(ab[i j])
    - r[i j] == r(ab[i j])
    - r(i j) == r(ab[i j])
    """


class CartesianDirectionToken(DirectionToken):
    """Class to parse direction specifiers in 3D space.

    CartesianDirectionToken specifiers are used to declare cartesian directions
    of movements or constraints for geometric displacements. They can specify
    one, two, or three dimensions of movements. The direction specifier can be
    given as one of:
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
    vectors_xyz : ndarray, shape (dof, 3)
        An array of one, two or three orthonormal vectors spanning the
        direction space. Note, the vectors are in the xyz convention.
        Use the `vectors_zxy` property to get them in the zxy (LEED)
        convention.

    Raises
    ------
    ValueError
        If the direction string is invalid or if the number of components does
        not match the number of directions.
    DirectionTokenParserError
        If the direction string contains unsupported components like
        azimuthal or radial directions, or if it is a fractional
        direction (use FractionalDirectionToken or make_direction_token
        for those).
    """

    def __init__(self, direction_str):
        _dir_str = direction_str.strip()
        if not _dir_str:
            raise DirectionTokenParserError('Empty direction token.')
        self.direction_str = _dir_str
        self._vectors, self.dof = self._parse_direction(_dir_str)

    def __eq__(self, other):
        """Compare two CartesianDirectionToken objects for equality."""
        if not isinstance(other, CartesianDirectionToken):
            return NotImplemented
        return self.dof == other.dof and np.allclose(
            self.vectors_xyz, other.vectors_xyz
        )

    def __str__(self):
        """Return a string representation of the CartesianDirectionToken object."""
        return f'CartesianDirectionToken(vectors={self.vectors_xyz}, dof={self.dof})'

    @property
    def vectors_xyz(self):
        """Return the vectors in xyz convention."""
        return self._vectors

    @property
    def vectors_zxy(self):
        """Return the vectors in zxy convention (LEED convention)."""
        return _xyz_to_zxy(self._vectors)

    def _parse_direction(self, direction_str):
        _check_unsupported_directions(direction_str)
        if _is_fractional(direction_str):
            msg = (
                f'{direction_str!r} is a fractional direction. Use '
                'FractionalDirectionToken or make_direction_token.'
            )
            raise DirectionTokenParserError(msg)

        # Vector direction like 'xyz[1 2 3]'
        if '[' in direction_str:
            match = re.match(CART_DIRECTION_PATTERN, direction_str)
            if not match:
                msg = f'Invalid direction format: {direction_str}'
                raise ValueError(msg)
            dirs = match['dir']
            _check_no_duplicate_directions(dirs, direction_str)
            vecs = [float(v) for v in match['vec'].split()]
            if len(dirs) != len(vecs):
                msg = (
                    'Mismatch between directions and components '
                    f'in {direction_str}'
                )
                raise ValueError(msg)
            vec = self._embed_vector(dirs, vecs)
            return np.array([self._normalize(vec)]), 1
        # Basis directions like 'x', 'xy', etc.
        # validate the direction string
        if not all(c in SIMPLE_DIRECTIONS for c in direction_str):
            msg = f'Invalid direction: {direction_str}'
            raise ValueError(msg)
        _check_no_duplicate_directions(direction_str, direction_str)
        if len(direction_str) > 3:
            msg = f'Direction string too long: {direction_str} (max 3 allowed)'
            raise ValueError(msg)
        vecs = [self._get_basis_vector(c) for c in direction_str]
        return np.array(vecs), len(vecs)

    @classmethod
    def _from_vectors(cls, direction_str, vectors_xyz):
        """Return a token for already orthonormal `vectors_xyz`."""
        token = cls.__new__(cls)
        token.direction_str = direction_str
        token._vectors = np.asarray(vectors_xyz, dtype=float)
        token.dof = len(token._vectors)
        return token

    @staticmethod
    def _get_basis_vector(label):
        return {
            'x': np.array([1, 0, 0]),
            'y': np.array([0, 1, 0]),
            'z': np.array([0, 0, 1]),
        }[label]

    @staticmethod
    def _embed_vector(dirs, comps):
        vec = np.zeros(3)
        idx_map = {'x': 0, 'y': 1, 'z': 2}
        for direction, component in zip(dirs, comps):
            vec[idx_map[direction]] = component
        return vec

    @staticmethod
    def _normalize(vec):
        norm = np.linalg.norm(vec)
        if norm < DISPLACEMENTS_FILE_EPS:
            msg = f'Zero-length vector: {vec}'
            raise ValueError(msg)
        return vec / norm


class FractionalDirectionToken(DirectionToken):
    """Class to parse directions in fractional (unit-cell) coordinates.

    Fractional directions are given in terms of the in-plane unit-cell
    vectors a and b. The direction specifier can be one of:
    - 'ab'
        The plane spanned by a and b.
    - 'a', 'b', 'ab[i j]', '[i j]'
        A 1D vector along i*a + j*b. 'a' is the same as 'ab[1 0]' and
        '[1 0]', 'b' the same as 'ab[0 1]' and '[0 1]'. The components
        'i' and 'j' may be floating-point numbers.

    Cartesian vectors are only available once the unit cell is known,
    via `to_cartesian`.

    Parameters
    ----------
    direction_str : str
        The direction string to parse.

    Attributes
    ----------
    dof : int
        The number of degrees of freedom (DOF), either 1 or 2.
    components : ndarray, shape (dof, 2)
        The fractional components of the vectors spanning the direction
        space. For a 1D direction, the components are normalized to
        unit length (the length of a direction is irrelevant).

    Raises
    ------
    DirectionTokenParserError
        If the direction string is empty, invalid, not fractional, or
        if the number of components is not two.
    """

    def __init__(self, direction_str):
        _dir_str = direction_str.strip()
        if not _dir_str:
            raise DirectionTokenParserError('Empty direction token.')
        self.direction_str = _dir_str
        self.components, self.dof = self._parse_direction(_dir_str)

    def __eq__(self, other):
        """Compare two FractionalDirectionToken objects for equality."""
        if not isinstance(other, FractionalDirectionToken):
            return NotImplemented
        return self.dof == other.dof and np.allclose(
            self.components, other.components
        )

    def __str__(self):
        """Return a string representation of this token."""
        return (
            f'FractionalDirectionToken(components={self.components}, '
            f'dof={self.dof})'
        )

    def to_cartesian(self, ab_cell):
        """Return the Cartesian direction token for a given unit cell.

        Parameters
        ----------
        ab_cell : Sequence
            Shape (2, 2). The in-plane unit cell, with the a and b
            vectors as columns (like `Slab.ab_cell`).

        Returns
        -------
        CartesianDirectionToken
            The same direction in Cartesian coordinates. Its
            `vectors_xyz` are orthonormal vectors spanning the
            direction space: for a 1D direction 'ab[i j]', the
            normalized i*a + j*b; for the 'ab' plane, the first
            vector is along a.
        """
        in_plane = self.components @ np.asarray(ab_cell, dtype=float).T
        vectors = np.zeros((self.dof, 3))
        vectors[:, :2] = _orthonormalize(in_plane)
        # pylint: disable-next=protected-access  # Same module
        return CartesianDirectionToken._from_vectors(
            self.direction_str,
            vectors,
        )

    def _parse_direction(self, direction_str):
        _check_unsupported_directions(direction_str)
        if not _is_fractional(direction_str):
            msg = f'Not a fractional direction: {direction_str}'
            raise DirectionTokenParserError(msg)
        if direction_str in FRACTIONAL_DIRECTIONS:  # 'a' or 'b'
            index = FRACTIONAL_DIRECTIONS.index(direction_str)
            return np.eye(2)[[index]], 1
        if direction_str == FRACTIONAL_PLANE:
            return np.eye(2), 2
        match = re.match(FRACTIONAL_DIRECTION_PATTERN, direction_str)
        if not match:
            msg = f'Invalid direction format: {direction_str}'
            raise DirectionTokenParserError(msg)
        try:
            comps = [float(v) for v in match['vec'].split()]
        except ValueError as exc:
            msg = f'Invalid direction components: {direction_str}'
            raise DirectionTokenParserError(msg) from exc
        if len(comps) != len(FRACTIONAL_DIRECTIONS):
            msg = (
                'Fractional directions need exactly two components. '
                f'Invalid direction: {direction_str}'
            )
            raise DirectionTokenParserError(msg)
        comps = np.array(comps)
        norm = np.linalg.norm(comps)
        if norm < DISPLACEMENTS_FILE_EPS:
            msg = f'Zero-length vector: {direction_str}'
            raise DirectionTokenParserError(msg)
        return np.array([comps / norm]), 1


def make_direction_token(direction_str):
    """Return the direction token that `direction_str` describes.

    Parameters
    ----------
    direction_str : str
        The direction string to parse. Cartesian directions ('x', 'xy',
        'xyz[1 2 3]', ...) give a CartesianDirectionToken, fractional
        ones ('a', 'b', 'ab', 'ab[1 2]', '[1 2]') a
        FractionalDirectionToken.

    Returns
    -------
    DirectionToken
        The parsed direction token.

    Raises
    ------
    DirectionTokenParserError
        If `direction_str` is empty or is a direction that is not
        supported (azimuthal, radial, along c, or mixing the Cartesian
        and fractional frames).
    ValueError
        If `direction_str` is a malformed Cartesian direction.
    """
    _dir_str = direction_str.strip()
    if not _dir_str:
        raise DirectionTokenParserError('Empty direction token.')
    _check_unsupported_directions(_dir_str)
    if _is_fractional(_dir_str):
        return FractionalDirectionToken(_dir_str)
    return CartesianDirectionToken(_dir_str)


def _check_unsupported_directions(direction_str):
    """Raise if `direction_str` is a direction we do not support."""
    # TODO: remove the azi/r check once we support azimuthal and
    # radial directions, and dispatch the appropriate token instead.
    if direction_str.startswith(('azi', 'r')):
        msg = (
            'Azimuthal and radial directions are not yet supported. '
            f'Invalid direction: {direction_str}'
        )
        raise DirectionTokenParserError(msg)
    frame = _frame_letters(direction_str)
    if UNSUPPORTED_C_DIRECTION in frame:
        msg = (
            'Directions involving the c unit-cell vector are not '
            f'supported. Invalid direction: {direction_str}'
        )
        raise DirectionTokenParserError(msg)
    is_fractional = any(c in FRACTIONAL_DIRECTIONS for c in frame)
    is_cartesian = any(c in SIMPLE_DIRECTIONS for c in frame)
    if is_fractional and is_cartesian:
        msg = (
            'Directions cannot mix Cartesian (x, y, z) and fractional '
            f'(a, b) components. Invalid direction: {direction_str}'
        )
        raise DirectionTokenParserError(msg)


def _check_no_duplicate_directions(letters, direction_str):
    """Raise ValueError if `letters` contains a letter more than once."""
    if len(set(letters)) != len(letters):
        msg = f'Duplicate directions in direction string: {direction_str}'
        raise ValueError(msg)


def _frame_letters(direction_str):
    """Return the letters before the '[' of `direction_str`."""
    return direction_str.split('[', maxsplit=1)[0].strip()


def _is_fractional(direction_str):
    """Return whether `direction_str` is in fractional coordinates."""
    frame = _frame_letters(direction_str)
    if not frame:  # Bare '[i j]'
        return direction_str.startswith('[')
    return all(c in FRACTIONAL_DIRECTIONS for c in frame)


def _orthonormalize(vectors):
    """Return orthonormal vectors spanning the same space as `vectors`.

    Gram-Schmidt on the rows of `vectors`, which keeps the first one
    parallel to its input.
    """
    result = []
    for vec in vectors:
        orthogonal = vec - sum(np.dot(vec, v) * v for v in result)
        result.append(orthogonal / np.linalg.norm(orthogonal))
    return np.array(result)


def _xyz_to_zxy(vecs):
    """Convert from xyz to zxy axis (LEED) convention."""
    # Apply axis permutation: x→y, y→z, z→x
    return vecs[:, [2, 0, 1]]  # roll last axis backward
