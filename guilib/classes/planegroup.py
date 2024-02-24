"""Module planegroup of viperleed.guilib.classes.

======================================
  ViPErLEED Graphical User Interface
======================================

Created: 2024-02-2
Author: Michele Riva (@michele-riva)

Defines the PlaneGroup class as well as a bunch of symmetry operations.
The PlaneGroup class used to be part of the guilib.base module.
"""

from ast import literal_eval
import itertools
import re
from typing import Iterable

import numpy as np

from viperleed.guilib.helpers import remove_duplicates
from viperleed.guilib.helpers import two_by_two_array_to_tuple

# Here some shorthands for symmetry operation matrices, expressed in
# FRACTIONAL coordinates. The naming convention, unless specified
# next to the operation, is:
# Rotations:
#    Cn -> 2pi/n counter-clockwise rotation
#    Cmn -> -2pi/n counter-clockwise rotation
# Mirrors for non-orthogonal bases (a/b are the basis unit vectors):
#    Mij: mirror across a line through vector v = i*a + j*b
#    Mimj: mirror across a line through vector v = i*a - j*b
# These two are good for all cells,...
E = (1, 0), (0, 1)       # Identity
C2 = (-1, 0), (0, -1)    # == -E
# ...these two are good for rectangular and square cells,...
Mx = (1, 0), (0, -1)     # Mirror across line through 1st unit vector
My = (-1, 0), (0, 1)     # Mirror across line through 2nd unit vector
# ...these are good for square cells,...
C4 = (0, -1), (1, 0)
Cm4 = (0, 1), (-1, 0)    # == -C4
M45 = (0, 1), (1, 0)     # Across line at +45deg relative to 1st vector
Mm45 = (0, -1), (-1, 0)  # Across line at -45deg relative to 1st vector
# ...these are good for rhombic and hex (both obtuse),...
M11 = (0, 1), (1, 0)     # == M45
M1m1 = (0, -1), (-1, 0)  # == -M11 == Mm45
M01 = (-1, -1), (0, 1)
M10 = (1, 0), (-1, -1)
# ...and these are good for hex only (obtuse).
C6 = (1, 1), (-1, 0)
Cm6 = (0, -1), (1, 1)
C3 = (0, 1), (-1, -1)
Cm3 = (-1, -1), (1, 0)
M21 = (1, 1), (0, -1)
M12 = (-1, 0), (1, 1)


_GROUPS_FOR_SHAPE = {
    'Oblique': ('p1', 'p2'),
    'Rectangular': (
        'p1', 'p2', 'pm', 'pm[1 0]', 'pm[0 1]', 'pg[1 0]', 'pg[0 1]',
        'rcm[1 0]', 'rcm[0 1]', 'pmm', 'pmg[1 0]', 'pmg[0 1]',
        'pgg', 'rcmm'
        ),
    'Square': (
        'p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'pg[1 0]', 'pg[0 1]',
        'cm[1 1]', 'cm[1 -1]', 'pmm', 'pmg[1 0]', 'pmg[0 1]',
        'pgg', 'cmm', 'p4', 'p4m', 'p4g'
        ),
    'Rhombic': (
        'p1', 'p2', 'cm[1 1]', 'cm[1 -1]', 'cmm'
        ),
    'Hexagonal': (
        'p1', 'p2', 'cm[0 1]', 'cm[1 0]', 'cm[1 1]', 'cm[1 -1]', 'cm[1 2]',
        'cm[2 1]', 'cmm', 'p3', 'p3m1', 'p31m', 'p6', 'p6m'
        )
    }
_GROUP_TO_OPS = {
    'p1': (E,),
    'p2': (E, C2),
    'pm[1 0]': (E, Mx), 'pm[0 1]': (E, My),
    'pg[1 0]': (E, Mx), 'pg[0 1]': (E, My),
    'cm[1 0]': (E, M10), 'cm[0 1]': (E, M01),
    'cm[1 1]': (E, M11), 'cm[1 -1]': (E, M1m1),
    'cm[1 2]': (E, M12), 'cm[2 1]': (E, M21),
    'rcm[1 0]': (E, Mx), 'rcm[0 1]': (E, My),
    'pmm': (E, Mx, My, C2),
    'pmg[1 0]': (E, Mx, My, C2), 'pmg[0 1]': (E, Mx, My, C2),
    'pgg': (E, C2, Mx, My),
    'cmm': (E, C2, M11, M1m1), 'rcmm': (E, C2, Mx, My),
    'p4': (E, C2, C4, Cm4),
    'p4m': (E, Mx, My, M45, Mm45, C2, C4, Cm4),
    'p4g': (E, Mx, My, M45, Mm45, C2, C4, Cm4),
    'p3': (E, C3, Cm3),
    'p3m1': (E, M12, M21, M1m1, C3, Cm3),
    'p31m': (E, M10, M11, M01, C3, Cm3),
    'p6': (E, C6, Cm6, C3, C2, Cm3),
    'p6m': (E, M10, M01, M12, M21, M11, M1m1, C6, Cm6, C3, C2, Cm3)
    }
_MAP_GLIDE_OPS = {  # For set_screws_glides. Direction to matrix
    (1, 0): M10,
    (0, 1): M01,
    'x': Mx,
    'y': My,
    (1, 1): M11,    # Also OK for +45deg on square cells
    (1, -1): M1m1,  # Also OK for -45deg on square cells
    (1, 2): M12,
    (2, 1): M21,
    }
_MAP_SCREW_OPS = {  # For set_screws_glides. Rot. orders to matrices
    2: (C2,),
    3: (C3, Cm3),
    4: (C4, C2, Cm4),
    6: (C6, C3, C2, Cm3, Cm6),
    }
_SUBGROUPS = {
    'p1': ('p1'),
    'p2': ('p1', 'p2'),
    'pm[1 0]': ('p1', 'pm[1 0]'),
    'pm[0 1]': ('p1', 'pm[0 1]'),
    'pg[1 0]': ('p1', 'pg[1 0]'),
    'pg[0 1]': ('p1', 'pg[0 1]'),
    'cm[1 0]': ('p1', 'cm[1 0]'),
    'cm[0 1]': ('p1', 'cm[0 1]'),
    'cm[1 1]': ('p1', 'cm[1 1]'),
    'cm[1 -1]': ('p1', 'cm[1 -1]'),
    'cm[1 2]': ('p1', 'cm[1 2]'),
    'cm[2 1]': ('p1', 'cm[2 1]'),
    'rcm[1 0]': ('p1', 'pm[1 0]', 'pg[1 0]', 'rcm[1 0]'),
    'rcm[0 1]': ('p1', 'pm[0 1]', 'pg[0 1]', 'rcm[0 1]'),
    'pmm': ('p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'pmm'),
    'pmg[1 0]': ('p1', 'p2', 'pg[1 0]', 'pm[0 1]', 'pmg[1 0]'),
    'pmg[0 1]': ('p1', 'p2', 'pg[0 1]', 'pm[1 0]', 'pmg[0 1]'),
    'pgg': ('p1', 'p2', 'pg[1 0]', 'pg[0 1]', 'pgg'),
    'cmm': ('p1', 'p2', 'cm[1 1]', 'cm[1 -1]', 'cmm'),
    'rcmm': ('p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'pg[1 0]', 'pg[0 1]', 'pgg',
             'rcm[1 0]', 'rcm[0 1]', 'pmm', 'pmg[1 0]', 'pmg[0 1]', 'rcmm'),
    'p4': ('p1', 'p2', 'p4'),
    'p4m': ('p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'cm[1 1]', 'cm[1 -1]', 'pmm',
            'cmm', 'p4', 'p4m'),
    'p4g': ('p1', 'p2', 'pg[1 0]', 'pg[0 1]', 'cm[1 1]', 'cm[1 -1]', 'pgg',
            'cmm', 'p4', 'p4g'),
    'p3': ('p1', 'p3'),
    'p3m1': ('p1', 'cm[1 -1]', 'cm[2 1]', 'cm[1 2]', 'p3', 'p3m1'),
    'p31m': ('p1', 'cm[1 0]', 'cm[0 1]', 'cm[1 1]', 'p3', 'p31m'),
    'p6': ('p1', 'p2', 'p3', 'p6'),
    'p6m': ('p1', 'p2', 'cm[1 0]', 'cm[0 1]', 'cm[1 1]', 'cm[1 -1]', 'cm[2 1]',
            'cm[1 2]', 'cmm[1 2]', 'cmm[1 0]', 'cmm[2 1]', 'cmm[0 1]',
            'cmm[1 -1]', 'cmm[1 1]', 'p3', 'p3m1', 'p31m', 'p6', 'p6m')
    }


class PlaneGroup:
    """Class representing a planar 2D group."""

    def __init__(self, group='p1'):
        """Initialize plane-group instance.

        Parameters
        ----------
        group: str or PlaneGroup, optional
            Hermann-Mauguin notation for the 2D plane group (with some
            extensions). Acceptable values:
            'p1', 'p2', 'pm[1 0]', 'pm[0 1]', 'pg[1 0]', 'pg[0 1]',
            'cm[1 0]', 'cm[0 1]', 'cm[1 1]', 'cm[1 -1]', 'cm[1 2]',
            'cm[2 1]', 'rcm[1 0]', 'rcm[0 1]', 'pmm', 'pmg[1 0]',
            'pmg[0 1]', 'pgg', 'cmm', 'rcmm', 'p4', 'p4m', 'p4g',
            'p3', 'p3m1', 'p31m', 'p6', 'p6m'. Default is 'p1'.
            See docs/_static/planegroups.pdf for more info.

        Raises
        ------
        TypeError
            If group is neither a string nor a PlaneGroup.
        ValueError
            If group (including its direction) is not an acceptable
            plane group.
        """
        if isinstance(group, PlaneGroup):
            bulk_3d = group.screws_glides
            group = group.group
        else:
            bulk_3d = tuple()
        group = self.__check_group_name(group)
        self.group = group

        # The next one will be a tuple of the 2x2 matrices representing
        # the isomorphism part of screws and glide planes perpendicular
        # to the surface
        self.__ops_3d = bulk_3d

    def __eq__(self, other):
        """Return whether self is equal to other.

        This is a relatively simple check, that only looks at whether
        the Hermann-Mauguin names of self and other are the same.

        Use self.same_operations(other, inlude_3d) to explicitly check
        the set of symmetry operations.
        """
        # Notice that Python 3 handles correctly (and in a
        # faster way) cases in which __ne__ is not reimplemented.
        # This is why __ne__ is not reimplemented here.
        if not isinstance(other, PlaneGroup):
            try:
                other = PlaneGroup(other)
            except (TypeError, ValueError):
                return NotImplemented
        return self.group == other.group

    def __repr__(self):
        """Return a string representation of PlaneGroup."""
        return f'PlaneGroup({self.group!r})'

    def __str__(self):
        """Return the Hermann-Mauguin name as a string."""
        return self.group

    @staticmethod
    def groups_compatible_with(cell_shape, operations=(), include_3d=False):
        """Return the groups compatible with cell_shape and operations.

        Given a certain cell shape, the groups compatible with
        that shape and a set of symmetry operations are those
        that are (i) groups that can exists for that cell shape,
        and (ii) also subgroups of the list of operations given.

        Parameters
        ----------
        cell_shape : {'Oblique', 'Rectangular', 'Square',
                      'Rhombic', 'Hexagonal'}
            The shape of the cell of the lattice. Used to initially
            pick the list of compatible groups.
        operations : iterable, optional
            Operation matrices. If not given or empty, the list
            selected from the shape is returned. Otherwise, only
            those groups in the list that are subgroups of the
            list of operations given are returned. If given, the
            elements of 'operations' are assumed to be represented
            in fractional coordinates. Notice that passing a
            generator will consume it. Default is an empty tuple.
        include_3d : bool, optional
            Whether also 3D operations (screw axes, glide planes)
            should be included in the check. Default is False.

        Returns
        -------
        tuple of strings
            The Hermann-Mauguin names of compatible groups
        """
        try:
            compatible_with_shape = _GROUPS_FOR_SHAPE[cell_shape]
        except KeyError as err:
            raise ValueError(
                f'PlaneGroup: invalid cell_shape {cell_shape}. Should be one '
                'of ' + ', '.join(repr(k) for k in _GROUPS_FOR_SHAPE)
                ) from err
        if not operations:
            return compatible_with_shape

        if not isinstance(operations, Iterable):
            raise TypeError('PlaneGroup: operations should be an iterable')

        # Convert any numpy array to tuples,
        # for the 'not in' comparison below
        operations = list(operations)
        for i, operation in enumerate(operations):
            if np.shape(operation) != (2, 2):
                raise ValueError('PlaneGroup: operations should be '
                                 'an iterable of (2, 2) matrices')
            if isinstance(operation, np.ndarray):
                operations[i] = two_by_two_array_to_tuple(operation)

        compatible = []
        for group in compatible_with_shape:
            group = PlaneGroup(group)
            if any(operation not in operations
                   for operation in group.operations(include_3d)):
                continue
            compatible.append(group.group)
        return tuple(compatible)

    @property
    def screws_glides(self):
        """Return 3D symmetry operations as 2x2 tuples.

        The 2x2 tuples of integers returned represent the isomorphic
        part of screw axes and glide planes perpendicular to the
        surface. Use .set_screws_glides to modify this attribute.

        Returns
        -------
        tuple of tuples
        """
        return self.__ops_3d

    def set_screws_glides(self, new_screws_glides, cell_shape=None):
        """Set the isomorphic part of 3D symmetry operations.

        Parameters
        ----------
        new_screws_glides : str or Sequence
            The specification to assign a new value to screws_glides.
            May be a sequence, containing 2x2 matrices of integers,
            or a string of the following forms (with or without spaces,
            not case-sensitive):
                'r(#, #, ...), m([#, #], [#, #], ...)'
                'm([#, #], [#, #], ...), r(#, #, ...)'
                'r(#, #, ...)'
                'm([#, #], [#, #], ...)'
                'None'
            For screws, the list in 'r()' provides the order of
            rotations.  For glides, an entry '[i,j]' means the
            glide plane leaves the in-plane direction i*a + j*b
            unmodified.
        cell_shape : {None, 'Oblique', 'Rectangular', 'Square',
                      'Hexagonal', 'Rhombic'}
            Cell shape of the lattice. This parameter is mandatory if
            new_screws_glides is a string and it defines glide planes.
            Default is None.

        Raises
        ------
        TypeError
            If new_screws_glides is neither a string nor an iterable.
        TypeError
            If new_screws_glides is a string specifying glide planes,
            but cell_shape was not given.
        ValueError
            If new_screws_glides does not match one of the acceptable
            string specifications
        ValueError
            If screw orders or glide-plane directions are not among the
            acceptable ones when specified as a string.
        ValueError
            If new_screws_glides is a non-string sequence, but it does
            not contain 2x2 integer matrices.
        """
        if not new_screws_glides:
            self.__ops_3d = tuple()
            return
        if isinstance(new_screws_glides, str):
            self._set_new_screws_glides_from_string(new_screws_glides,
                                                    cell_shape)
        elif isinstance(new_screws_glides, Iterable):
            self._set_new_screws_glides_from_matrices(new_screws_glides)
        else:
            raise TypeError(
                f'{type(self).__name__}.set_screws_glides: Invalid input. '
                f'Must be string or Iterable, not {type(input).__name__}'
                )

    def _set_new_screws_glides_from_string(self, screws_glides, cell_shape):
        """Assign a new value to .screws_glides from a string specification."""
        if screws_glides.lower() == 'none':
            self.__ops_3d = tuple()
            return
        screws_glides = re.sub(r'\s+', '', screws_glides)
        self.__ops_3d = remove_duplicates(
            itertools.chain(self._parse_screws_spec(screws_glides),
                            self._parse_glide_spec(screws_glides, cell_shape)),
            return_type=tuple
            )
        if not self.__ops_3d:
            raise ValueError(f'{type(self).__name__}.set_screws_glides: '
                             'Invalid input.')

    @classmethod
    def _parse_screws_spec(cls, screws_spec):
        """Yield screw-axes point operations from a string specification."""
        # Pattern: an 'r()' token containing a single digit, optionally
        # followed by more single digits, comma separated. We accept
        # multiple tokens in screws_spec.
        screw_re = re.compile(r'R\((\d(?:,\d)*)\)', re.I)
        try:
            screw_orders = literal_eval(
                ','.join(screw_re.findall(screws_spec)) + ','
                )
        except (ValueError, TypeError, SyntaxError,
                MemoryError, RecursionError):
            # Either no match, or some problem in the screws_spec
            return
        screws = (screw
                  for order in screw_orders
                  for screw in _MAP_SCREW_OPS[order])
        try:
            yield from screws
        except KeyError as exc:
            raise ValueError(f'{cls.__name__}.set_screws_glides: Invalid '
                             f'rotation order {exc} in the input. Only 2-, '
                             '3-, 4-, and 6-fold orders allowed.') from None

    @classmethod
    def _parse_glide_spec(cls, glide_spec, cell_shape):
        """Yield glide-plane point operations from a string specification."""
        # Pattern: an 'm()' token containing one pair of '[N,N]',
        # optionally followed by more pairs of '[N,N]', comma separated
        # Here N is any positive or negative single-digit number.
        # Notice the use of non-capturing groups for the inner '[N,N]'
        # to avoid them as separate tokens in findall.
        glide_re = re.compile(r'M\(((?:\[-?\d,-?\d\])(?:,\[-?\d,-?\d\])*)\)',
                              re.I)
        try:
            glides = literal_eval(
                ','.join(glide_re.findall(glide_spec)) + ','
                )
        except (ValueError, TypeError, SyntaxError,
                MemoryError, RecursionError):
            # Either no match, or some problem in the glide_spec
            return
        if not cell_shape:
            raise TypeError(f'{cls.__name__}.set_screws_glides: cell '
                            'shape is required when glide planes are given.')
        is_orthogonal = cell_shape in {'Square', 'Rectangular'}
        for direction in glides:
            if direction[0] < 0 or not direction[0] and direction[1] < 0:
                direction = [-d for d in direction]
            # The only directions that require special attention are
            # [1, 0] and [0, 1], since the matrices depend on the shape
            # of the cell
            if direction == [1, 0] and is_orthogonal:
                key = 'x'
            elif direction == [0, 1] and is_orthogonal:
                key = 'y'
            else:
                key = tuple(direction)
            try:
                yield _MAP_GLIDE_OPS[key]
            except KeyError:
                raise ValueError(
                    f'{cls.__name__}.set_screws_glides: invalid '
                    f'direction {direction} for glide plane. The only allowed '
                    f'directions are: {_MAP_GLIDE_OPS.keys() - {"x", "y"}}'
                    ) from None

    def _set_new_screws_glides_from_matrices(self, matrices):
        """Assign new screws_glides from a sequence of 2x2 integer matrices."""
        matrices = np.asarray(matrices)
        if len(matrices.shape) != 3 or matrices.shape[1:] != (2, 2):
            raise ValueError(
                f'{type(self).__name__}.set_screws_glides: a sequence input '
                'should be a 1D "list" of 2x2 "matrices". Found incompatible '
                f'shape {matrices.shape}.'
                )
        if np.any(abs(matrices - matrices.round()) > 1e-4):
            raise ValueError(
                f'{type(self).__name__}.set_screws_glides: a sequence input '
                'should contain only integer-valued matrices.'
                )

        matrices = matrices.round().astype(int)
        # Finally, check determinants
        determinants = (abs(np.linalg.det(op)) for op in matrices)
        if any(abs(d - 1) > 1e-4 for d in determinants):
            raise ValueError(
                f'{type(self).__name__}.set_screws_glides: a sequence '
                'input should contain uni-determinant matrices.'
                )
        self.__ops_3d = remove_duplicates(
            (two_by_two_array_to_tuple(op) for op in matrices),
            return_type=tuple
            )

    @property
    def subgroups(self):
        """Return the 2D subgroups of self as a set of strings."""
        return set(_SUBGROUPS[self.group])

    @classmethod
    def highest_symmetry_for_shape(cls, cell_shape):
        """Return the highest symmetry group compatible with cell_shape."""
        *_, high_sym = cls.groups_compatible_with(cell_shape)
        return cls(high_sym)

    def is_valid_group(self, group, cell_shape):
        """Checks if group is a valid group for a given cell_shape.

        No type checking is done on group other than string
        and PlaneGroup, under the assumption that group will
        be used to create a PlaneGroup itself, which does
        the type checking in the constructor.

        Parameters
        ----------
        group : str or PlaneGroup
            The group to be checked. When a string, assume
            it is the Hermann-Mauguin name.
        cell_shape : {'Oblique', 'Rectangular', 'Square',
                      'Hexagonal', 'Rhombic'}
            Shape of the lattice for which compatibility is
            checked.

        Returns
        -------
        valid_group : bool
        """
        valid_for_shape = self.groups_compatible_with(cell_shape)
        if isinstance(group, str):
            group = self.__check_group_name(group)
            return group in valid_for_shape
        if isinstance(group, PlaneGroup):
            return group.group in valid_for_shape
        return False

    def operations(self, include_3d=False):
        """Return point symmetry operations as 2x2 matrices.

        Parameters
        ----------
        include_3d : bool
            Whether the returned 'list' should include also
            the isomorphic part of the 3D screws and glides
            perpendicular to the surface.

        Returns
        -------
        operations : tuple
            Items are 2x2 tuples representing the isomorphic part
            (i.e., excluding translations for screws/glides) of
            the symmetry operations of this group. The matrices
            are expressed in 'fractional' coordinates.
        """
        ops = list(_GROUP_TO_OPS[self.group])
        if include_3d:
            ops.extend(self.screws_glides)
        return tuple(ops)

    def same_operations(self, other, include_3d=False):
        """Return whether self has the same operations as other.

        The comparison does not check whether the Hermann-Mauguin
        names are the same. Use self == other to test that. Thus
        PlaneGroup('pm[1, 0]').same_operations(PlaneGroup('pg[1, 0]'))
        returns True.

        Parameters
        ----------
        other : str or PlaneGroup
            The other group whose operations are to be compared to
            theo ones of self.
        include_3d : bool, optional
            Whether the comparison should also include the isomorphic
            part of the 3D screw/glide operations. Default is False.

        Returns
        -------
        bool
        """
        if not isinstance(other, PlaneGroup):
            try:
                other = PlaneGroup(other)
            except (TypeError, ValueError):
                return False
        self_ops = set(self.operations(include_3d))
        other_ops = set(other.operations(include_3d))
        return self_ops == other_ops

    def transform(self, transform, inverse=None, include_3d=False):
        """Transform the group operations to new coordinates.

        No assumption is made on the coordinate transform.
        This means that the operation matrices returned may
        have non-integer values.

        Parameters
        ----------
        transform : Sequence
            Shape == (2, 2). The coordinate transformation to be
            applied. This is interpreted as follows: assume that
            the operations are currently expressed in a coordinate
            system with basis B_0 = (a_0, b_0), with a_0 and b_0
            row vectors. The new basis is
            B_1 = (a_1, b_1) = transform @ B_0.
        inverse : Sequence, optional
            Shape == (2, 2). The inverse of the transform. This can
            be given for performance reasons. If not given or None,
            inverse = numpy.linalg.inv(transform).  If given, it is
            checked that transform @ inverse is the identity matrix.
            Default is None.
        include_3d : bbol, optional
            Whether the operations transformed should also
            include the isomorphic part of the 3D screws/glides.

        Returns
        -------
        tuple of numpy.ndarrays
        """
        if np.shape(transform) != (2, 2):
            raise ValueError('PlaneGroup.transform requires a 2x2 '
                             'array-like as coordinate transform matrix. '
                             f'Found shape {np.shape(transform)} instead.')
        if inverse is None:
            inverse = np.linalg.inv(transform)
        elif np.shape(inverse) != (2, 2):
            raise ValueError('PlaneGroup.transform requires a 2x2 array-like '
                             'as the inverse of the coordinate transform. '
                             f'Found shape {np.shape(transform)} instead.')
        elif not np.allclose(np.dot(transform, inverse), ((1, 0), (0, 1))):
            raise ValueError('PlaneGroup.transform transformation matrix and '
                             'inverse are inconsistent.')
        return tuple(np.linalg.multi_dot((transform,
                                          op,
                                          inverse))
                     for op in self.operations(include_3d))

    def __check_group_name(self, group):
        """Check that a string is an acceptable Hermann-Mauguin name.

        Before checking, spaces are fixed to allow checking directions.

        Parameters
        ----------
        group : str
            String to be fixed and checked.

        Returns
        -------
        fixed_group : str
            A Hermann-Mauguin with correct spacings, if the
            string given was an acceptable PlaneGroup name.

        Raises
        ------
        TypeError
            If group is not a string
        ValueError
            If group is not an acceptable Hermann-Mauguin name.
        """
        if not isinstance(group, str):
            raise TypeError('"group" should be a string. '
                            f'Found {type(group).__name__} instead')
        group_re = re.compile(
            r'''
            (?P<hermann>[\w]+)      # Hermann-Mauguin
            (?:[\[]\s*              # Optional direction opening bracket
            (?P<dir1> [-]*[\d]+)    # First direction
            [\s]*                   # Optional space
            (?P<dir2> [-]*[\d]+)    # Second direction
            [\]])*                  # Optional direction closing bracket
            ''',
            re.VERBOSE)

        group_match = group_re.match(group)
        if group_match is None:
            raise ValueError(f'{group} is not an acceptable plane group.')
        group = group_match.group('hermann')

        if group_match.group('dir1') is not None:
            group = (f'{group}[{group_match.group("dir1")} '
                     f'{group_match.group("dir2")}]')

        if group not in _GROUP_TO_OPS:
            raise ValueError(f'{group} is not an acceptable plane group.')

        return group
