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

import numpy as np

from viperleed.guilib.base import check_type
from viperleed.guilib.classes.planegroup import PlaneGroup


class Lattice2D:
    """
    Lattice2D(basis, space='real', group='p1', limit=1)

    Class for a 2D lattice. Includes basis and its shape, plane group, whether
    it is a real- or reciprocal-space lattice, and an array of lattice points.

    Parameters
    ----------
    basis : 2x2 array-like
        basis vectors a and b of the lattice, with a = basis[0], b=basis[1]
        the units are assumed to be Angstrom for real space lattices, and
        2*pi/Angstrom for reciprocal-space lattices
    space : str, default='real'
        accepts 'real' or 'reciprocal' for real and reciprocal space
        lattices, respectively
    group : str, default='p1'
        plane group in Hermann-Mauguin notation. See also the PlaneGroup
        class for acceptable values.
    limit : int, default=1
        radius used to limit the number of lattice points generated. Only
        lattice points closer to the origin than limit will be produced.

    Attributes
    ----------
    basis : 2x2 numpy.ndarray
        basis, same convention as with class instantiation

    space : str, read-only
        'real' or 'reciprocal'

    cell_shape : str, read-only
        shape of the unit cell. Can be 'Oblique', 'Rectangular',
        'Square', 'Rhombic', or 'Hexagonal'

    group : PlaneGroup
        a PlaneGroup instance representing the plane group of the lattice.
        Can be set with Lattice2D.group

    lattice : numpy.ndarray, shape == (..., 2)
        array of lattice points (x, y) [for real-space] or (gx, gy) [for
        reciprocal space]

    hk : numpy.ndarray, shape == (..., 2)
        array of (h, k) indices that generate the points in lattice, i.e.,
        lattice[i] = hk[i, 0]*basis[0] + hk[i, 1]*basis[1]

    n_beams : int, read-only
        Return number of lattice points

    real_basis : 2x2 numpy.ndarray, read only
        returns a copy of the real-space basis, independently
        of whether the lattice is real or reciprocal. Use
        Lattice2D.basis for changing the basis

    reciprocal_basis : 2x2 numpy.ndarray, read-only
        returns the reciprocal of Lattice2D.basis.
        Use Lattice2D.basis for changing the basis

    lattice_parameters : 3-tuple of floats, read-only
        return lattice parameters as (length_a, length_b, angle), where
        angle is the angle between the basis vectors in degrees. Use
        Lattice2D.basis for changing the basis

    special_directions : list of (2,) numpy.ndarrays
        Returns vectors along the directions of the mirror planes, in
        the same coordinate system as Lattice2D.basis

    Public methods
    --------------
    get_rotated_lattice(angle)
        returns a copy of Lattice2D.lattice rotated by angle (degrees),
        without actually modifying Lattice2D

    get_rotated_basis(angle)
        returns a copy of Lattice2D.basis rotated by angle, without
        actually modifying Lattice2D

    high_symm_transform()
        return matrix transform that gives highest symmetry basis
        (with the same lattice). Does not transform the Lattice2D
        though. Use Lattice2D.make_high_symmetry() to make high
        symmetry, as this also returns the same transform

    make_high_symmetry()
        transform Lattice2D so the basis has the highest possible
        symmetry, returns the transformation used

    transform(transform, as_copy=True)
        transforms the Lattice2D by applying transform TO THE LEFT of
        basis, i.e., the transform matrix passed is meant to act on ROW
        VECTORS. Can either modify Lattice2D itself (as_copy=False) or
        a copy of it. In any case, it returns the transformed Lattice2D
        instance

    equivalent_points(special_direction=None, superlattice=None)
        Returns a dictionary with (h, k): star of (h, k) entries.
        If special_direction is given, only the mirror plane
        containing it is used to compute the star of each lattice
        point. Use superlattice to return indices expressed with
        respect to a different basis

    Private methods
    ---------------
    __get_cell_shape()
        finds the shape of lattice unit cell. Use Lattice2D.cell_shape,
        unless you suspect that the shape is not up to date

    __generate_lattice()
        generates the lattice points
    """

    def __init__(self, basis, space='real', group='p1', limit=1):
        if not check_type(basis, 'arraylike'):
            raise TypeError(f"Lattice2D basis should be array-like, "
                            f"not {type(basis)}")
        if not np.shape(basis) == (2, 2):
            raise ValueError("Lattice2D basis needs to have a (2, 2) shape. "
                             f"Found {np.shape(basis)} instead")
        if not check_type(space, 'str'):
            raise TypeError(f"Lattice2D space should be a string, "
                            f"not {type(space)}")
        if space not in ('real', 'reciprocal'):
            raise ValueError(f"Lattice2D space {space} unknown")
        if not check_type(limit, 'number'):
            raise TypeError(f"Lattice2D limit should be a scalar, "
                            f"not {type(limit)}")

        self._basis = np.asarray(basis)
        self._space = space
        self._shape = self.__get_cell_shape()  # __get_cell_shape

        # check if the plane group given is consistent with the cell shape
        if not PlaneGroup.is_valid_group(group, self.cell_shape):
            raise ValueError(f"Lattice2D: invalid group {group} for lattice "
                             f"shape {self.cell_shape}")
        self._group = PlaneGroup(group)
        self._limit = limit
        self.lattice, self.hk = self.__generate_lattice()

    def __repr__(self):
        txt = (f"{self.cell_shape} "
               + f"viperleed.Lattice2D({self.basis}, ".replace('\n', '')
               + f"space={self.space}, group={self.group}, "
               + f"limit={self._limit})")
        return txt

    @property
    def basis(self):
        """
        Returns the lattice basis as a 2x2 numpy.ndarray, with a = basis[0] and
        b = basis[1]
        """
        return self._basis

    @basis.setter
    def basis(self, basis):
        """
        Sets the lattice basis to basis and updates the other attributes

        Parameters
        ----------
        basis : 2x2 array-like
        """
        if not check_type(basis, 'arraylike'):
            raise TypeError("basis must be array-like. "
                            f"Found {type(basis)} instead.")
        basis = np.asarray(basis)
        if basis.shape != (2, 2):
            raise ValueError("Lattice2D basis needs to have a (2, 2) shape. "
                             f"Found {np.shape(basis)} instead.")
        if abs(np.linalg.det(basis)) < 1e-5:
            raise ValueError("Lattice2D basis cannot be a singular matrix.")

        self._basis = basis
        self._shape = self.__get_cell_shape()
        self.lattice, self.hk = self.__generate_lattice()

        compatible_groups = self.group.groups_compatible_with(self.cell_shape)
        if self.group.group not in compatible_groups:
            # Shape does not allow the old group.
            # Can't pick one, so use 'p1'                                      # TODO: is there a better way to do this? How to treat 3D operations (also when group is one of those of the shape)?
            self.group = 'p1'

    @property
    def cell_shape(self):
        """
        Returns the shape of the unit cell as a string. The value returned is
        'Oblique', 'Rectangular', 'Square', 'Rhombic', or 'Hexagonal'.
        """
        return self._shape

    @property
    def group(self):
        """Return the PlaneGroup of the lattice."""
        return self._group

    @group.setter
    def group(self, group):
        """
        Set the plane group to a PlaneGroup instance with Hermann-Mauguin
        symbol equal to group. See PlaneGroup for a list of acceptable input
        parameters.
        """
        self._group = PlaneGroup(group)

    @property
    def lattice_parameters(self):
        """
        Returns the length of the lattice vectors and the angle alpha between
        them as a tuple (norm_a, norm_b, alpha)
        """
        norm_a, norm_b = np.linalg.norm(self.basis, axis=-1)
        alpha = np.arccos(np.dot(self.basis[0], self.basis[1])/(norm_a*norm_b))
        return norm_a, norm_b, np.degrees(alpha)

    @property
    def n_beams(self):
        """
        Returns
        -------
        int
        the number of lattice points (for real-space) or LEED beams (for
        reciprocal-space)
        """
        return len(self.hk)

    @property
    def real_basis(self):
        """
        Always returns the real-space basis of the lattice as a 2x2 np.array
        """
        if self.space == 'reciprocal':
            return self.reciprocal_basis
        # Not sure why this was like this, but seems useless. Probably was
        # just to yield a copy.
        # return np.dot(self.basis, [[1, 0], [0, 1]])
        return self.basis.copy()

    @property
    def reciprocal_basis(self):
        """
        Returns the reciprocal of the lattice basis. Notice that if
        self.space == 'real', this returns the reciprocal space basis.
        If self.space == 'reciprocal', it returns the reciprocal of
        the reciprocal, i.e., the real-lattice basis.

        This is different from the return value of real_basis
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
        """
        Returns the space of the lattice as a string. The value returned is
        either 'real' or 'reciprocal'. This property cannot be changed after
        initialization.
        """
        return self._space

    @property
    def special_directions(self):
        """
        Returns vectors along the directions of the mirror planes, in the
        same coordinate system as self.real_basis. The list returned is as long
        as there are operations in self.group, with None instead of a direction
        vector in correspondence of rotations

        Returns
        -------
        list of numpy.ndarrays
        """
        # Work on the real-basis coordinates, as directions are the
        # same in real and reciprocal space. Using a reciprocal-space basis
        # would require to transpose-invert the symmetry operations.
        # Also notice that we're not including the 3d operations, as this
        # makes sense only for a 'surface' lattice, and the role of 3d
        # operations is merely that of determining how many domains may
        # occur in LEED << IS THIS CORRECT??
        transform = np.linalg.inv(self.real_basis)
        ops = self.group.transform(transform)

        # The special directions are those parallel to the eigenvectors
        # of the mirrors with eigenvalue == 1 (use 1e-4 tolerance).
        # In fact, this is the direction that is unchanged when applying
        # the mirror operation. Rotations do not give any special direction.
        directions = []
        for op in ops:
            if np.linalg.det(op) > 0:  # Rotation
                directions.append(None)
                continue
            eigval, eigvecs = np.linalg.eig(op)
            # eigvecs has the eigenvectors as columns, ordered in the same
            # way as the eigenvalues in eigvals. Notice the use of .extend()
            # rather than .append(). This is because eigvecs is a matrix,
            # and its 1-element slice is a (1, 2) array. Using .extend() adds
            # a (2,)-shaped array.
            directions.extend(eigvecs.T[np.abs(eigval - 1) < 1e-4])
        return directions

    def __get_cell_shape(self):
        """
        Determine the shape of the lattice ('Oblique', 'Rhombic', 'Hexagonal',
        'Rectangular', or 'Square'), also changing the real-space basis in an
        handedness-conserving fashion so that the angle between the basis
        vectors in real space is obtuse (for hexagonal and rhombic only)

        Returns
        -------
        string
        """

        basis = self._basis

        # Relative tolerance factor within which things are assumed to be equal
        eps = 1e-3

        # cosine of angle between vectors
        cosine = np.dot(basis[0], basis[1])/(np.linalg.norm(basis[0])
                                             * np.linalg.norm(basis[1]))

        # Mismatch between the length of the two vectors
        delta = np.linalg.norm(basis[0])/np.linalg.norm(basis[1]) - 1
        if abs(cosine) < eps:  # angle is 90°
            if np.abs(delta) < eps:
                return 'Square'
            return 'Rectangular'
        if np.abs(delta) < eps:  # rhombic or hex
            if self.space == 'real':
                if cosine > eps:  # angle is acute -> make it obtuse
                    print('Warning: the input might be corrupted,'
                          'since the real-space basis is not obtuse!')
                    transform = (0, -1), (1, 0)  # this keeps the handedness
                else:
                    transform = (1, 0), (0, 1)
                self._basis = np.dot(transform, self._basis)
                basis = self._basis
                cosine = np.dot(basis[0],
                                basis[1])/(np.linalg.norm(basis[0])
                                           * np.linalg.norm(basis[1]))
            else:
                if cosine > eps:
                    # The reciprocal lattice will be acute, since
                    # the real one is obtuse.
                    cosine *= -1

            if abs(cosine + 0.5) < eps:  # angle is 120deg
                return 'Hexagonal'
            return 'Rhombic'
        return 'Oblique'

    def __generate_lattice(self):
        """
        Generates a list of lattice points given a basis.

        Parameters
        ----------
        limit : scalar
            Determines which portion of the lattice is generated.
            For space == 'real', the lattice is generated up to a radius
            of 1.5*limit.
                One should then plot from -limit to +limit. This should
                cover any post-rotation of the lattice that the user
                might later request.
            For space == 'reciprocal', the lattice is generated up to a
            radius of limit

        basis : 2x2 array-like, default = self.basis
            Contains the unit vectors as basis[0] and basis[1]

        space : string, {'real', 'reciprocal'}, default = self.space
            Determines whether a real-space or a reciprocal-space
            lattice is generated.
            This only affects the behavior of limit (see above)

        Returns
        -------
        Tuple (lat, hk)

        lat : 1D numpy.ndarray
            Lattice points generated

        hk : 1D numpy.ndarray, hk.shape == lat.shape
            Integer indices of the lattice points generated.
            The lattice points are lat = h*basis[0] + k*basis[1].
            hk[i] is a 2-element array of indices of lattice point i
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
        hk = np.mgrid[-h_max:h_max+1, -h_max:h_max+1].reshape(2,-1).T

        # Create lattice:
        # Notice that, given a row vector of indices (h, k) and a basis in
        # matrix form B = [[a1, a2],[b1, b2]], the corresponding lattice
        # point can be obtained as the row vector
        #    L = [L1, L2] = [h, k]*B
        lattice = np.dot(hk, basis)

        # Now find all those lattice points that lie within limit, and use
        # this as a mask for the output
        mask = np.linalg.norm(lattice, axis=1) <= limit

        return lattice[mask], hk[mask]

    def get_rotated_lattice(self, angle):
        """
        Returns a copy of Lattice2D.lattice rotated by angle

        Parameters
        ----------
        angle: float
               Rotation angle in degrees. Positive values will give a
               counterclockwise rotation

        Returns
        -------
        np.array, same shape as Lattice2D.lattice
        """
        angle = np.radians(angle)

        rot = np.array([[np.cos(angle), np.sin(angle)],
                        [-np.sin(angle), np.cos(angle)]])

        return np.dot(self.lattice, rot)

    def get_rotated_basis(self, angle):
        """
        Returns a rotated copy of the lattice basis

        Parameters
        ----------
        angle: float
               Rotation angle in degrees. Positive values will give a
               counterclockwise rotation

        Returns
        -------
        rotated basis: 2x2 np.array
        """
        angle = np.radians(angle)
        rot = np.array([[np.cos(angle), np.sin(angle)],
                        [-np.sin(angle), np.cos(angle)]])
        return np.dot(self.basis, rot)

    def high_symm_transform(self):
        """
        Returns a 2x2 np.array that brings the basis lattice into the highest
        possible symmetry configuration when left-multiplied to the lattice
        basis. The highest-symmetry basis generates the same lattice as the
        any other basis.

        The transform can bring an oblique lattice into square, rectangular,
        hexagonal or rhombic, or make the basis vectors as close to orthogonal
        as possible.

        This function DOES NOT transform the lattice. Call make_high_symmetry()
        for that.
        """
        # The following line should also redefine the basis so that it is
        # obtuse for real-space rhombic and hexagonal lattices
        _shape = self.cell_shape

        # Will always work on the real-space lattice for convenience,
        # then convert back to the reciprocal one in case the lattice was
        # reciprocal in the first place
        basis = self.real_basis

        # In what follows, t_elem is used to define a specific elementary
        # operation to be performed on the lattice basis. This is
        # left-multiplied to t_overall at each elementary step, so that
        # t_overall contains the overall transformation

        if _shape != 'Oblique':
            return (1, 0), (0, 1)

        # Lattice is oblique.
        # Transform lattice to have the shortest two vectors, with angle
        # closest to 90°.
        # This might bring it to rect, hex or rhombic.
        #
        # If neither, will anyway transform to have the closest to rect.

        # ALGORITHM for reduction to closest to rect:
        # This is a discrete version of Gram-Schmidt's algorithm to find
        # orthogonal bases
        # At each iteration:
        #   - order vectors by norm, the shortest first
        #   - determine the projection of the second on the first,
        #     and calculate the nearest integer kk
        #   - subtract from the second the projection calculated above
        #   - check whether now the second is the smallest.
        #     If yes, repeat, otherwise finished.

        # swap keeps track of whether the first and second vectors are
        # swapped at the end of this passage
        swap = (1, 0), (0, 1)
        t_overall = (1, 0), (0, 1)
        while True:
            # Swap vectors if needed to get the shortest first
            if np.linalg.norm(basis[0]) > np.linalg.norm(basis[1]):
                t_elem = (0, 1), (1, 0)
            else:
                t_elem = (1, 0), (0, 1)
            swap = np.dot(t_elem, swap)
            t_overall = np.dot(t_elem, t_overall)
            basis = np.dot(t_elem, basis)
            projection = np.dot(basis[0], basis[1])/np.dot(basis[0],
                                                           basis[0])
            projection = int(np.round(projection))
            t_elem = (1, 0), (-projection, 1)
            t_overall = np.dot(t_elem, t_overall)
            basis = np.dot(t_elem, basis)
            if np.linalg.norm(basis[0]) <= np.linalg.norm(basis[1]):
                break
        # Swap vectors back if they were overall swapped
        t_overall = np.dot(swap, t_overall)
        basis = np.dot(swap, basis)

        # END OF ALGORITHM. Now the lattice is closest to rectangular.
        # It might be still any shape (square, rect, hex, rhombic, oblique)

        # Create a dummy lattice with the new basis,
        # to check which shape it has
        _shape = Lattice2D(basis).cell_shape

        # If it's still oblique, try to see if it can be transformed to hex
        # or rhombic by choosing "a" not to be the shortest vector of all.
        #
        # If possible, keep the new transformation. Otherwise, stick to the
        # one that makes it closest to rectangular
        #
        # All the operations that follow are stored in a matrix t_second,
        # to be later left-multiplied to t_overall to get the full
        # transformation
        #
        if _shape != 'ob':
            t_second = (1, 0), (0, 1)
        else:
            # lattice is still oblique, even if closest to rectangular
            #
            # Re-swapping guarantees that that the matrix has on the first
            # line the shortest possible vector,
            # and on the second line the second shortest possible vector.
            #
            # The only possible combinations that can lead to a
            # rhombic/hex are a'=b+a or a'=b-a, depending
            # on whether the angle is acute or obtuse, respectively
            #
            t_second = swap
            basis = np.dot(swap, basis)

            t_elem = [[-int(np.sign(np.dot(basis[0], basis[1]))), 1],
                      [0, 1]]
            t_second = np.dot(t_elem, t_second)
            basis = np.dot(t_elem, basis)
            dummy2 = Lattice2D(basis)

            # The following line might change acute into obtuse
            # --> check if it was the case
            # PROBABLY IT CANNOT HAPPEN ANYWAY!! IT SHOULD RATHER BE
            # CHECKED ON DUMMY ABOVE
            _shape = dummy2.cell_shape
            sign_before = np.dot(basis[0], basis[1])
            sign_after = np.dot(dummy2.basis[0], dummy2.basis[1])
            if sign_before*sign_after < 0:
                # sign did change -> basis went from acute to obtuse
                t_second = np.dot([[0, -1], [1, 0]], t_second)

            if _shape == 'ob':
                # lattice is still oblique, no transformation is needed
                # (will keep the one closest to rect)
                t_second = (1, 0), (0, 1)
        t_overall = np.dot(t_second, t_overall)

        # Finally update the transformation matrix to account for the correct
        # space of the lattice
        if self.space == 'reciprocal':
            t_overall = np.linalg.inv(t_overall).T
        return t_overall

    def is_high_symmetry(self):
        """Check whether the lattice has the highest symmetry possible."""
        return np.array_equal(((1, 0), (0, 1)), self.high_symm_transform())

    def make_high_symmetry(self):
        """
        Transform lattice basis to the one that realizes the highest possible
        symmetry and update the value of the attributes. Notice that if the
        basis can be transformed to higher symmetry, determining whether its
        plane group has higher symmetry than before is not possible.

        Returns
        -------
        np.array, shape = (2, 2)
        transformation matrix that has been left-multiplied to the basis
        """
        transform = self.high_symm_transform()

        if not np.array_equal([[1, 0], [0, 1]], transform):
            # lattice can be higher symmetry, i.e., it was oblique
            self.basis = np.dot(transform, self.basis)
            self._shape = self.__get_cell_shape()
            self.lattice, self.hk = self.__generate_lattice()
            # no need to change the group, since group was at most p2.
            # Inferring whether the group has higher symmetry is not possible

        return transform

    def transform(self, transform, as_copy=True):
        """
        Modifies the basis and lattice points according to transform, and
        returns a copy (if copy=True) or directly modifies the values of self.

        Returns
        -------
        Lattice2D
            Either a reference to self (copy=False) or a reference
            to a new Lattice2D. In both cases the lattice returned
            is updated with the transform
        """
        new_lattice = self
        if as_copy:
            new_lattice = copy.deepcopy(self)
        new_lattice.basis = np.dot(transform, self.basis)
        return new_lattice
