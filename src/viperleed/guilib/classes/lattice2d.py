"""Module lattice2d of viperleed.guilib.classes.

======================================
  ViPErLEED Graphical User Interface
======================================

Created on 2024-03-05

@author: Michele Riva (@michele-riva)

Defines the Lattice2D class, representing a two-dimensional lattice in
real or reciprocal space, with a symmetry group. Was originally part of
viperleed.guilib.base.
"""

import copy
from numbers import Real
from typing import Sequence

import numpy as np

from viperleed.calc.lib.matrix import rotation_matrix
from viperleed.guilib.classes import planegroup
from viperleed.guilib.classes.planegroup import PlaneGroup
from viperleed.guilib.helpers import array_to_string

ACUTE_TO_OBTUSE = np.array(((0, -1), (1, 0)))  # Conserves handedness
_SPECIAL_DIRECTIONS = {
    # These are the eigenvectors of mirror operations with unity
    # eigenvalue. In other words, they are the directions that
    # are left invariant by the mirror operations. Since the
    # group operations are expressed in fractional coordinates,
    # the special directions are the directions of the mirror
    # planes in fractional coordinates.
    planegroup.Mx: np.array([1, 0]),
    planegroup.My: np.array([0, 1]),
    planegroup.M11: np.array([1, 1])/2**0.5,    # Also M45
    planegroup.M1m1: np.array([-1, 1])/2**0.5,  # Also Mm45
    # Notice that the following directions appear to be swapped,
    # but are indeed right once considering that they are to be
    # pre-multiplied with the inverse of an hexagonal basis. See
    # comments in the special_directions @property for derivation.
    planegroup.M01: np.array([-1, 2])/5**0.5,
    planegroup.M10: np.array([2, -1])/5**0.5,
    planegroup.M21: np.array([1, 0]),
    planegroup.M12: np.array([0, 1]),
    }


# Disable as pylint considers twice private attributes and
# property when the property is set somewhere in the code
# pylint: disable-next=too-many-instance-attributes
class Lattice2D:
    """A two-dimensional lattice, in real or reciprocal space.

    Includes basis and its shape, plane group and whether
    it is a real- or reciprocal-space lattice.

    Attributes
    ----------
    lattice : numpy.ndarray, shape == (..., 2)
        Array of lattice points. (x, y) for real-space,
        (g_x, g_y) for reciprocal space
    hk : numpy.ndarray, shape == (..., 2)
        Array of (h, k) indices that generate the points
        in `lattice`, i.e., lattice[i] = hk[i] @ basis
    """

    # Disable due to pylint bug. Does not accept
    # 'optional' after list of allowed values.
    # pylint: disable-next=missing-type-doc
    def __init__(self, basis, space='real', group=None, limit=1):
        """Initialize Lattice2D instance.

        Parameters
        ----------
        basis : Sequence or numpy.ndarray
            Shape (2, 2). Basis vectors a and b of the lattice, with
            a == basis[0], b == basis[1]. The units are assumed to be
            Angstrom for real-space lattices, and 2*pi/Angstrom for
            reciprocal-space ones.
        space : {'real', 'reciprocal'}, optional
            Whether the lattice is a real- or reciprocal-space one.
            Default is 'real'.
        group : str or PlaneGroup or None, optional
            Plane group in Hermann-Mauguin notation. If None or
            not given, use PlaneGroup('p1'). Default is None.
        limit : int, optional
            Radius (in the same units as `basis`) used to limit the
            number of lattice points generated. Only lattice points
            closer to the origin than `limit` will be produced.
            Default is 1.

        Raises
        ------
        TypeError
            If `limit` is not a numerical value, or `basis` is not
            a numeric sequence.
        ValueError
            If `group` is not one of those acceptable for the cell
            shape given by `basis`.
        ValueError
            If `basis` is singular, or does not have shape (2, 2).
        ValueError
            If `space` is not one of 'real' or 'reciprocal'.
        """
        if space not in {'real', 'reciprocal'}:
            raise ValueError(f'{type(self).__name__}: unknown space {space!r}')
        if not isinstance(limit, Real):
            raise TypeError(f'{type(self).__name__}: limit should be a scalar,'
                            f' not {type(limit).__name__}')

        self._basis = None  # Set via .basis setter below
        self._group = None  # Set via .group setter below
        # pylint: disable-next=invalid-name  # For hk
        self.lattice, self.hk = None, None   # Set via .basis setter
        self._limit = limit
        self._shape = None       # Set via .basis setter below
        self._space = space
        self._was_acute = False  # Set via .basis setter below

        # Delegate checks to setters
        self.basis = basis   # Also updates cell shape
        self.group = group

    def __str__(self):
        """Return a string version of this Lattice2D."""
        return (f'{self.cell_shape}, {self.space}-space {type(self).__name__}'
                f'({array_to_string(self.basis)}, group={self.group})')

    def __repr__(self):
        """Return a string representation of this Lattice2D."""
        cls_name = type(self).__name__
        return (f'{cls_name}({array_to_string(self.basis)}, '
                f'space={self.space!r}, group={self.group.group!r}, '
                f'limit={self._limit})')

    @property
    def basis(self):
        """Return a 2x2 numpy.ndarray of basis vectors as rows."""
        return self._basis

    @basis.setter
    def basis(self, basis):
        """Set new basis vectors (rows) for this Lattice2D."""
        if (isinstance(basis, str)
                or not isinstance(basis, (Sequence, np.ndarray))):
            raise TypeError(
                f'{type(self).__name__}: basis must be a non-string '
                f'sequence. Found {type(basis).__name__!r} instead.'
                )
        basis = np.asarray(basis)
        if basis.shape != (2, 2):
            raise ValueError(f'{type(self).__name__}: invalid basis '
                             f'shape={basis.shape} must be (2, 2).')
        if abs(np.linalg.det(basis)) < 1e-5:
            raise ValueError(f'{type(self).__name__}: basis is singular.')

        self._basis = basis
        self._was_acute = False
        self._shape = self._get_cell_shape()  # Sets _was_acute
        self.lattice, self.hk = self._generate_lattice()

        if self.group is None:
            return

        if not PlaneGroup.is_valid_group(self.group, self.cell_shape):
            # Shape does not allow the old group.
            # Can't pick one, so use 'p1'                                       # TODO: is there a better way to do this? How to treat 3D operations (also when group is one of those of the shape)?
            self.group = 'p1'

    @property
    def cell_shape(self):
        """Return the shape of the unit cell.

        Returns
        -------
        shape : {'Oblique', 'Rectangular', 'Square', 'Rhombic', 'Hexagonal'}
        """
        return self._shape

    @property
    def group(self):
        """Return the PlaneGroup of this Lattice2D."""
        return self._group

    @group.setter
    def group(self, group):
        """Assign a new group (str or PlaneGroup) for this Lattice2D."""
        if group is None:
            self._group = PlaneGroup('p1')
            return
        if not PlaneGroup.is_valid_group(group, self.cell_shape):
            raise ValueError(f'{type(self).__name__}: invalid group {group} '
                             f'for lattice shape {self.cell_shape}')
        self._group = PlaneGroup(group)

    @property
    def lattice_parameters(self):
        """Return the lattice parameters of this Lattice2D.

        Returns
        -------
        length_a: float
            The length of the first basis vector in units of Angstrom
            or 2pi/Angstrom depending on self.space.
        length_b : float
            The length of the second basis vector in units of Angstrom
            or 2pi/Angstrom depending on self.space.
        alpha : float
            Angle between the two basis vectors in degrees.
        """
        norm_a, norm_b = np.linalg.norm(self.basis, axis=-1)
        alpha = np.arccos(np.dot(*self.basis) / (norm_a * norm_b))
        return norm_a, norm_b, np.degrees(alpha)

    @property
    def n_beams(self):
        """Return the number of LEED beams. Only for reciprocal-space."""
        if self.space == 'real':
            raise AttributeError('Real-space lattice has no LEED '
                                 'beams. Use .points instead.')
        return len(self.hk)

    @property
    def n_points(self):
        """Return the number of lattice points. Only for real-space."""
        if self.space == 'reciprocal':
            raise AttributeError('Reciprocal-space lattice has no '
                                 'points. Use .n_beams instead.')
        return len(self.hk)

    @property
    def real_basis(self):
        """Return a copy of the real-space basis of this Lattice2D.

        Notice that this is **always** the real-space basis of the
        lattice, independently of whether this lattice is real or
        reciprocal.

        Returns
        -------
        numpy.ndarray
            The real-space basis of this lattice.
        """
        if self.space == 'reciprocal':
            return self.reciprocal_basis
        return self.basis.copy()

    @property
    def reciprocal_basis(self):
        """Return the reciprocal of this lattice's basis.

        Returns
        -------
        reciprocal_basis : numpy.ndarray
            The reciprocal of this lattice's basis. Notice that if
            self.space == 'real', `reciprocal_basis` the reciprocal-
            lattice basis. However, if self.space == 'reciprocal',
            `reciprocal_basis` is the reciprocal of the reciprocal,
            i.e., the real-lattice basis.
        """
        # Working out the cross products
        #   a* = 2pi/area [b x z][0:2] = 2pi/area [b12, -b11]
        #   b* = 2pi/area [z x a][0:2] = 2pi/area [-a12, a11]
        # where
        #   area = det(basis).
        # Thus
        #  | a* |                    | b12  -a12| T
        #  |    | = 2pi/det(basis) * |          |   = 2pi * basis^(-T)
        #  | b* |                    |-b11   a11|
        return 2*np.pi*np.linalg.inv(self.basis).T

    @property
    def space(self):
        """Return the space ('real' or 'reciprocal') of this Lattice2D."""
        return self._space

    @property
    def special_directions(self):
        """Return vectors along the mirror planes and None for rotations.

        Returns
        -------
        list
            Items are in the same order as the operations of
            self.group. Each item corresponding to a rotation
            is None. Each item corresponding to a mirror plane
            is a (2,)-shaped numpy.ndarray parallel to the
            mirror plane, in the same coordinate system as
            self.real_basis.
        """
        # The special directions are those parallel to the eigenvectors
        # of the mirrors with eigenvalue == 1. We don't need to compute
        # the eigenvectors every time, though. It is sufficient to use
        # those of the group operations. A group operation expressed in
        # the coordinates of the real basis B is (see also
        # PlaneGroup.transform):
        #     R_b = inv(B) @ R @ B = Q_b @ L @ inv(Q_b),
        # where R is the fractional group operation, Q_b is the matrix
        # of eigenvectors in the coordinates of the basis, and L the
        # eigenvalue matrix. This means that
        #     R = (B @ Q_b) @ L @ inv(B @ Q_b) = Q_R @ L @ inv(Q_R),
        # that is, the eigenvector matrix in the frame of the basis
        # Q_b can be computed from the one of the operation, Q_R, as
        #     Q_b = inv(B) @ Q_R.
        # Then eigenvectors (i.e., the columns of Q_) are
        #     v_b = inv(B) @ v_R.
        # We work on the real-basis coordinates, as directions are the
        # same in real and reciprocal space. Using the reciprocal-
        # space basis would require to transpose-invert the symmetry
        # operations. Also notice that we're not including the 3d
        # operations, as special_directions make sense only for a
        # 'surface' lattice, and the role of 3d operations is merely
        # that of determining how many domains may occur in LEED                # TODO: IS THIS CORRECT??
        transform = np.linalg.inv(self.real_basis)
        fractional_directions = (_SPECIAL_DIRECTIONS.get(op, None)
                                 for op in self.group.operations())
        directions = (transform.dot(frac) if frac is not None else frac
                      for frac in fractional_directions)
        return [d/np.linalg.norm(d) if d is not None else d
                for d in directions]

    @property
    def was_acute(self):
        """Return whether the real-space basis was made obtuse."""
        return self._was_acute

    def _get_cell_shape(self):
        """Determine the shape of this lattice's basis.

        For hexagonal and rhombic lattices, this method also changes
        the real-space basis in an handedness-conserving fashion so
        that the angle between the basis vectors in real space is
        obtuse.

        Returns
        -------
        shape : {'Oblique', 'Rhombic', 'Hexagonal', 'Rectangular', 'Square'}
        """
        eps = 1e-3  # Relative tolerance for equality
        norm_a, norm_b, alpha = self.lattice_parameters
        cosine = np.cos(np.radians(alpha))
        delta = norm_a/norm_b - 1  # Relative mismatch between a and b
        if abs(cosine) < eps:  # angle is 90°
            return 'Square' if abs(delta) < eps else 'Rectangular'
        if abs(delta) >= eps:
            return 'Oblique'

        # Rhombic or hexagonal. Make real-space acute into obtuse.
        real_space_is_acute = (cosine > eps if self.space == 'real'
                               else cosine < -eps)
        self._was_acute = real_space_is_acute
        if real_space_is_acute:
            self._basis = ACUTE_TO_OBTUSE.dot(self._basis)
            cosine *= -1
        if self.space == 'reciprocal' and cosine > eps:
            # Reciprocal acute, thus real is obtuse
            cosine *= -1
        if abs(cosine + 0.5) < eps:  # angle is 120deg
            return 'Hexagonal'
        return 'Rhombic'

    def _generate_lattice(self):                                                # TODO: have to rethink the usage of limit, as we need more beams for large off-normal incidence.
        """Generate a list of lattice points from self.basis.

        Returns
        -------
        lattice : numpy.ndarray
            Lattice points. Shape == (N, 2), where N is determined
            by self._limit.
        hk : numpy.ndarray
            Integer indices of the lattice points generated,
            i.e., lattice[i] = hk[i] @ basis. Shape == (N, 2).
        """
        limit = self._limit
        basis = self.basis
        space = self.space

        if space == 'real':
            limit *= 1.5

        # get limit for the loops that follow
        shortest = min(*np.linalg.norm(basis, axis=-1),
                       np.linalg.norm(basis[0] + basis[1])/2,
                       np.linalg.norm(basis[0] - basis[1])/2)
        h_max = int(np.ceil(limit/shortest))

        # create grid of indices
        indices = np.mgrid[-h_max:h_max+1, -h_max:h_max+1].reshape(2, -1).T

        # Create lattice:
        # Notice that, given a row vector of indices (h, k) and
        # a basis in matrix form B = [[a1, a2], [b1, b2]], the
        # corresponding lattice point can be obtained as the
        # row vector
        #    L = [L1, L2] = [h, k] @ B
        lattice = np.dot(indices, basis)

        # Now find all those lattice points that lie within
        # limit, and use this as a mask for the output
        mask = np.linalg.norm(lattice, axis=1) <= limit
        return lattice[mask], indices[mask]

    def get_rotated_lattice_points(self, angle):
        """Return a copy of self.lattice rotated by angle.

        Parameters
        ----------
        angle : float
            Rotation angle in degrees. Positive values
            will give a counter-clockwise rotation.

        Returns
        -------
        rotated_lattice : numpy.ndarray
            The rotated version of self.lattice. Shape == (N, 2).
        """
        rot = rotation_matrix(np.radians(angle)).T
        return np.dot(self.lattice, rot)

    def get_rotated_basis(self, angle):
        """Return a rotated copy of this lattice's basis.

        Parameters
        ----------
        angle : float
            Rotation angle in degrees. Positive values
            will give a counter-clockwise rotation.

        Returns
        -------
        rotated basis: numpy.ndarray
            Shape == (2, 2). Unit vectors are rows.
        """
        rot = rotation_matrix(np.radians(angle)).T
        return np.dot(self.basis, rot)

    def high_symm_transform(self):
        """Return a matrix transform that gives highest-symmetry basis.

        This method **does not** transform this Lattice2D.
        Call make_high_symmetry() for that.

        Returns
        -------
        transform : numpy.ndarray
            Shape (2, 2). The transformation matrix that brings the
            basis of this lattice to its highest possible symmetry
            configuration. The high-symmetry basis is obtained by
            left-multiplying self.basis with `transform`.

        Notes
        -----
        The highest-symmetry basis generates the same lattice as
            any other basis.
        The `transform` can bring an oblique lattice into square,
            rectangular, hexagonal or rhombic, or make the basis
            vectors as close to orthogonal as possible. Lattices
            that are not oblique are considered to be already in
            their highest-symmetry form.
        The `transform` returned conserves the unit cell area.
        """
        if self.cell_shape != 'Oblique':
            # Nothing to do if it's already non-oblique
            return np.eye(2, dtype=int)

        # Will always work on the real-space lattice for convenience,
        # then convert back to the reciprocal one in case the lattice
        # was reciprocal in the first place
        basis = self.real_basis

        # As a first step, transform the basis to have the shortest
        # two vectors with angle closest to 90°. This might bring it
        # to rectangular, hexagonal or rhombic. If neither, will
        # anyway transform to have the closest to rectangular.
        # ALGORITHM for reduction to closest to rectangular:
        # This is a discrete version of Gram-Schmidt's algorithm
        # to find orthogonal bases.
        # At each iteration:
        #   - order vectors by norm, the shortest first
        #   - determine the projection of the second on
        #     the first, and calculate the nearest integer
        #   - subtract from the second the projection calculated above
        #   - check whether now the second is the smallest.
        #     If yes, repeat, otherwise finished.
        # In what follows, `t_elem` is used to define a specific
        # elementary operation to be performed on the lattice basis.
        # This is left-multiplied to `t_overall` at each elementary
        # step, so that `t_overall` contains the overall transform.
        # `swap` keeps track of whether the first and second vectors
        # are swapped at the end of this passage
        swap = np.eye(2, dtype=int)
        t_overall = np.eye(2, dtype=int)
        a_norm, b_norm = np.linalg.norm(basis, axis=1)
        while True:
            # Swap vectors if needed to get the shortest first
            t_elem = ((0, 1), (1, 0)) if a_norm > b_norm else ((1, 0), (0, 1))
            swap = np.dot(t_elem, swap)
            t_overall = np.dot(t_elem, t_overall)
            basis = np.dot(t_elem, basis)

            projection = np.dot(basis[0], basis[1])/np.dot(basis[0], basis[0])
            t_elem = (1, 0), (-round(projection), 1)
            t_overall = np.dot(t_elem, t_overall)
            basis = np.dot(t_elem, basis)

            a_norm, b_norm = np.linalg.norm(basis, axis=1)
            if a_norm <= b_norm:
                break
        # Swap vectors back if they were overall swapped
        t_overall = swap.dot(t_overall)
        basis = swap.dot(basis)

        # END OF ALGORITHM. Now basis is closest to rectangular. It
        # might be still any shape (square, rectangular, hexagonal,
        # rhombic, or oblique).

        # Use a dummy lattice to check the shape of the new basis
        dummy = Lattice2D(basis)
        if dummy.was_acute:
            # Also update the transformation matrix
            t_overall = ACUTE_TO_OBTUSE.dot(t_overall)

        # If it's still oblique, try to see if it can be transformed
        # to hex or rhombic by choosing "a" not to be the shortest
        # vector of all. If possible, keep the new transformation.
        # Otherwise, stick to the one that makes it closest to
        # rectangular.
        # All the operations that follow are stored in a matrix
        # `t_second`, to be later left-multiplied to `t_overall`
        # to get the full transformation
        if dummy.cell_shape != 'Oblique':
            t_second = np.eye(2, dtype=int)
        else:
            # Lattice is still oblique, even if closest to rectangular.
            # See if it can be made rhombic or hexagonal.
            basis = swap.dot(basis)  # Ensure basis[0] is shortest
            t_second = swap

            # The only possible combinations that can lead to a
            # rhombic/hexagonal are a'=b+a or a'=b-a, depending
            # on whether the angle is obtuse or acute, respectively
            _is_acute = np.dot(*basis) > 0
            t_elem = (
                (-1 if _is_acute else 1, 1),
                (0, 1)
                )
            t_second = np.dot(t_elem, t_second)
            basis = np.dot(t_elem, basis)
            dummy2 = Lattice2D(basis)  # May turn acute to obtuse again
            if dummy2.was_acute:
                # Also update the transformation matrix
                t_second = ACUTE_TO_OBTUSE.dot(t_second)

            if dummy2.cell_shape == 'Oblique':
                # Lattice is still oblique, no transformation is
                # needed. Will keep the one closest to rectangular
                t_second = np.eye(2, dtype=int)
        t_overall = np.dot(t_second, t_overall)

        # Finally update the transformation matrix to
        # account for the correct space of the lattice
        if self.space == 'reciprocal':
            t_overall = np.linalg.inv(t_overall).T
        return t_overall

    def is_high_symmetry(self):
        """Check whether the lattice has the highest symmetry possible."""
        return np.array_equal(np.eye(2), self.high_symm_transform())

    def make_high_symmetry(self):
        """Make the basis of this lattice highest symmetry.

        Internal attributes are also updated. Notice that if the
        basis can be transformed to higher symmetry, determining
        whether its plane group has higher symmetry than before
        is not possible.

        Returns
        -------
        transform : numpy.ndarray
            Transformation matrix that has been left-multiplied to
            the basis to make it highest symmetry. Shape (2, 2).
        """
        transform = self.high_symm_transform()
        if not np.array_equal(np.eye(2), transform):
            # Was oblique, and can be made high symmetry
            self.transform(transform)
        return transform

    def transform(self, transform):
        """Apply transform to this lattice's basis and lattice points.

        Parameters
        ----------
        transform : Sequence
            Shape (2, 2). The transformation matrix to be
            left-multiplied to the basis. This means that
            `transform` should be expressed in 'fractional'
            coordinates. This is normally one of the
            `planegroup` operations that are suited
            for self.cell_shape.

        Returns
        -------
        None.
        """
        self.basis = np.dot(transform, self.basis)

    def transformed(self, transform):
        """Return a new Lattice2D with transformed basis.

        Parameters
        ----------
        transform : Sequence
            Shape (2, 2). The transformation matrix to be
            left-multiplied to the basis. This means that
            `transform` should be expressed in 'fractional'
            coordinates. This is normally one of the
            `planegroup` operations that are suited
            for self.cell_shape.

        Returns
        -------
        transformed_lattice : Lattice2D
        """
        new_lattice = copy.deepcopy(self)
        new_lattice.transform(transform)
        return new_lattice
