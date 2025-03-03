"""Module bulk_slab of viperleed.calc.classes.slab.

Defines the BulkSlab class, a BaseSlab describing a 3D-periodic
crystal. This module was created as part of the refactoring of
slab.py (originally created on 2019-06-13).
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-02-21'
__license__ = 'GPLv3+'

import copy
import functools
import itertools

import numpy as np
from scipy.spatial.distance import cdist as euclid_distance

from viperleed.calc.lib import leedbase
from viperleed.calc.lib.coordinates import add_edges_and_corners
from viperleed.calc.lib.coordinates import collapse
from viperleed.calc.lib.itertools_utils import cycle
from viperleed.calc.lib.matrix import rotation_matrix_order

from .base_slab import BaseSlab
from .errors import AlreadyMinimalError
from .errors import MissingSublayersError
from .errors import SlabError


class BulkSlab(BaseSlab):
    """A class representing an infinite solid, period in 3D.

    Contains unit cell, element information and atom coordinates.
    Also has a variety of convenience functions for manipulating
    and updating the atoms.

    Attributes
    ----------
    ucell : np.array
        The unit cell as vectors a, b, c (columns)
    poscar_scaling : float
        The original scaling factor from POSCAR
    elements : tuple of str
        Element labels as read from POSCAR
    chemelem : set of str
        Chemical elements in the slab, including from `ELEMENT_MIX`
    n_per_elem : dict {str: int}
        The number of atoms per POSCAR element.
    atlist : AtomList
        List of all atoms in the slab.
    layers : tuple of Layer
        Each `layer` is a composite of sublayers, as in TensErLEED
    sublayers : tuple of SubLayer
        Each SubLayer contains atoms of equal element and Z coordinate
    sitelist : list of Sitetype
        List of distinct sites as Sitetype, storing information
        on vibration and concentration
    ucell_mod : list of tuples (str, numpy.ndarray)
        Stored modifications made to the unit cell; each is a tuple
        of (type, array), where type is 'lmul', 'rmul', 'add', or
        'c_shift'.
    topat_ori_z : float
        Stores the original position of the topmost atom in Cartesian
        coordinates
    celltype : str                                                              # TODO: would be nicer with an Enum
        Unit-cell shape as string. Values: 'oblique', 'rhombic',
        'rectangular', 'square', 'hexagonal'
    planegroup : str
        Symmetry group of the slab. May be reduced by the user
        relative to `foundplanegroup`.
    foundplanegroup : str
        Highest symmetry found. Doesn't get modified when user
        reduces symmetry manually.
    orisymplane : SymPlane
        Only stored if the `planegroup` is ambiguous as to which unit
        vector the symmetry plane at the origin is parallel to
    linklists : list of list of Atom
        List of lists of atoms which are linked by a symmetry operation
    bulk_screws : list of int
        Integer list of rotation orders present in the bulk.
    bulk_glides : list of SymPlane
        List of glide-symmetry planes present in the bulk.
    """

    def __init__(self):
        """Initialize instance."""
        super().__init__()
        del self.bulkslab
        self.bulk_screws = []
        self.bulk_glides = []

    @property
    def is_bulk(self):
        """Return whether this is a bulk slab."""
        return True

    @property
    def smallest_interlayer_gap(self):
        """Return the smallest z gap between two adjacent layers.

        Make sure to update_layer_coordinates() before.

        Returns
        -------
        min_dist : float
            The smallest of the z distances between adjacent layers.
            Distances are calculated between the bottommost atom of
            the higher layer and the topmost atom of the lower one.

        Raises
        ------
        MissingLayersError
            If no layers are available.
        """
        if self.n_layers == 1:
            return self.c_vector[2] - self.layers[0].thickness
        return super().smallest_interlayer_gap

    @classmethod
    def from_slab(cls, other):
        """Return a `cls` instance with attributes deep-copied from `other`.

        Parameters
        ----------
        other : BaseSlab
            The slab whose attributes are to be copied.

        Returns
        -------
        new_slab : SurfaceSlab
            A new slab instance with attributes copied from
            `other`. Notice that no modification will occur
            for new_slab.layers. This means that some of the
            .layers of `new_slab` may be non-bulk. The same
            holds true for atoms. The caller is responsible
            for modifying layer.is_bulk after this call.
        """
        return super().from_slab(other)

    # Disabled too-many-arguments below because 7/5 seem better than
    # packing these arguments into some data structure which would
    # make the call to this method harder to follow.
    # pylint: disable-next=too-many-arguments
    def apply_bulk_cell_reduction(self, eps, epsz=None,
                                  new_c_vec=None,
                                  new_ab_cell=None,
                                  recenter=True,
                                  z_periodic=True):
        """Reduce bulk unit cell in- and/or out-of-plane.

        Extra atoms left after reduction of the unit cell size
        are removed.

        If desired, and the unit cell is actually changed, the
        atom that is currently at the top can be kept at the top
        also later. In this case, the z position of the new unit
        cell is also centred relative to the geometric centre of
        the atom coordinates.

        Parameters
        ----------
        eps : float
            Cartesian tolerance (Angstrom) for in-plane comparison
            of atomic coordinates.
        epsz : float or None, optional
            Cartesian tolerance (Angstrom) for comparison of atomic
            coordinates in the direction orthogonal to the surface.
            If not given or None, it is taken equal to `eps`. Default
            is None.
        new_c_vec : Sequence or None, optional
            Cartesian coordinates of the new repeat vector. If not
            given or None, the current repeat vector is retained.
            Shape (3,). Default is None.
        new_ab_cell : Sequence or None, optional
            Cartesian coordinates of the new surface unit cell. The
            new unit vectors are rows, i.e., a, b = `new_ab_cell`.
            If not given or None, the current surface unit cell is
            retained. Default is None.
        recenter : bool, optional
            Whether atom coordinates along the c axis should be
            modified so that the atom that was highest when this
            method was called (or one equivalent to it upon
            translation) is also highest afterwards. Default
            is True.
        z_periodic : bool, optional
            Whether the new_c_vec assigned will make this BulkSlab
            periodic along the z direction. Default is True.

        Returns
        -------
        None.
        """
        if new_ab_cell is None and new_c_vec is None:
            return  # Nothing to do

        if epsz is None:
            epsz = eps
        update_origin = new_c_vec is not None

        # Keep track of the atom that is currently at the top:
        # after reducing the unit-cell vector along c, the fractional
        # positions are screwed. We rely on the Cartesian ones, and,
        # if requested, will 'shift' the fractional ones such that
        # the topmost atom now is still the topmost atom later
        top_atom = self.top_atom if recenter else None

        # Make sure Cartesians are up to date,
        # then reduce c direction if needed
        self.update_cartesian_from_fractional()
        if new_c_vec is not None:
            self.c_vector[:] = new_c_vec

        # Reduce in-plane dimension
        if new_ab_cell is not None:
            self.ab_cell[:] = new_ab_cell.T

        # Collapse the atom coordinates to the new cell. We collapse
        # the coordinates in a different manner if this slab is turned
        # into a periodic one: atoms close to the edges in z are
        # collapsed towards c=0. This prevents atoms closer to the
        # edges than eps/epsz from being assigned different sublayers
        # later on.
        if not z_periodic:
            self.collapse_cartesian_coordinates(update_origin=update_origin)
        else:
            # Here we do a version of collapse_cartesian_coordinates
            # that accounts for eps.
            releps = np.dot((eps, eps, epsz), np.linalg.inv(self.ucell.T))
            self.update_fractional_from_cartesian()
            self.collapse_fractional_coordinates(releps=releps)
            self.update_cartesian_from_fractional(update_origin=update_origin)

        # Get rid of duplicates resulting from the size reduction
        self.remove_duplicate_atoms(eps, epsz)

        if not recenter:
            return

        # As mentioned above, fractional positions in c
        # are now screwed. Make top atom the topmost again...
        top_atom_cfrac = top_atom.pos[2]
        for atom in self:
            atom.pos[2] = (atom.pos[2] - top_atom_cfrac + 0.9999) % 1.0
        # ...then center fractional c coordinates around cell midpoint
        atoms_c_frac = [at.pos[2] for at in self]
        midpos = (max(atoms_c_frac) + min(atoms_c_frac)) / 2
        for atom in self:
            atom.pos[2] = (atom.pos[2] - midpos + 0.5) % 1.0
        self.update_cartesian_from_fractional(update_origin=update_origin)

    def ensure_minimal_c_vector(self, rpars, z_periodic=True):
        """Reduce the c vector to its minimum value.

        This method not only finds the shortest c vector, but also
        reduces the extension of this slab to have the new c vector
        as the out-of-plane repeat. It also modifies `rpars` with
        the minimal vector found. Use `get_minimal_c_vector` if you
        only want to detect the shortest repeat, without affecting
        either the slab or `rpars`.

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS read from file. Attributes used:
            read-only: SYMMETRY_EPS.
            overwritten: BULK_REPEAT. The minimal c vector
            found is stored in the BULK_REPEAT attribute of
            `rpars`.
        z_periodic : bool, optional
            Whether the slab is to be considered periodic along
            the direction perpendicular to the surface while the
            c vector is minimized. Use False, unless the current
            c vector is already a repeat vector for this slab.
            Default is True.

        Returns
        -------
        shortest_c : numpy.ndarray
            The new shortest c vector found.

        Raises
        ------
        AlreadyMinimalError
            If no shorter c vector is found.
        SlabError
            If this slab has too large of a vacuum gap at the bottom.
        """
        self.create_sublayers(rpars.SYMMETRY_EPS.z)

        # May raise AlreadyMinimalError
        rpars.BULK_REPEAT = self.get_minimal_c_vector(rpars.SYMMETRY_EPS,
                                                      rpars.SYMMETRY_EPS.z,
                                                      z_periodic=z_periodic)
        bottom_atom = self.bottom_atom
        delta_z = self.topat_ori_z - bottom_atom.cartpos[2]                     # TODO: flip with .cartpos[2] -- Issue #174
        if not z_periodic and delta_z > rpars.SYMMETRY_EPS.z:
            # Here we could also remove_vacuum_at_bottom ourselves, but
            # it is a bit of a pain as we would need to shift things
            # back up later in order to keep the z centering.
            raise SlabError(
                f'{type(self).__name__}.ensure_minimal_c_vector: cannot reduce'
                ' c vector for a non-z-periodic slab that has a vacuum gap at '
                f'the bottom (gap={delta_z:.3f} > {rpars.SYMMETRY_EPS.z:.3f}).'
                f' Shift all atoms down along z using '
                'slab.remove_vacuum_at_bottom(rpars)'
                )

        # Now we're sure that this slab will become
        # periodic after modifying the c vector
        self.apply_bulk_cell_reduction(rpars.SYMMETRY_EPS,
                                       epsz=rpars.SYMMETRY_EPS.z,
                                       new_c_vec=rpars.BULK_REPEAT,
                                       z_periodic=True)
        self.create_sublayers(rpars.SYMMETRY_EPS.z)
        return rpars.BULK_REPEAT

    def get_bulk_3d_str(self):
        """Return info about bulk screw axes and glide planes as a string.

        Returns
        -------
        bulk_3d_str : str
            Format is 'r(2, 4), m([1,1], [ 1,-1])' if there is
            any screw axes or glide planes, otherwise 'None'.
        """
        return leedbase.bulk_3d_string(self.bulk_screws,
                                       (p.par for p in self.bulk_glides))

    def get_bulk_repeat(self, rpars, only_z_distance=False):
        """Return the bulk repeat vector (with positive z).

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS to be interpreted. This is unused for
            a bulk slab.
        only_z_distance : bool, optional
            Whether a distance in the direction perpendicular to the
            surface (i.e., not necessarily along the c axis) should
            be returned rather than a full vector. This is ignored
            if `rpars.BULK_REPEAT` is a vector. Default is False.

        Returns
        -------
        bulk_repeat_vector : numpy.ndarray or float
            Bulk repeat vector pointing from the bulk to the surface,
            or its component along z.
        """
        bulkc = self.c_vector.copy()
        if bulkc[2] < 0:
            bulkc *= -1

        if isinstance(rpars.BULK_REPEAT, np.ndarray) or not only_z_distance:
            return bulkc
        return bulkc[2]

    def get_candidate_layer_periods(self, epsz=1e-4):
        """Find at which offsets the sublayers repeat in a bulk slab.

        Parameters
        ----------
        epsz : float, optional
            Tolerance (Cartesian) on z distances. Default is 1e-4.

        Returns
        -------
        candidate_periods : list
            Each element is an integer corresponding to the
            number of sublayers that potentially constitute
            a period for the slab. They are 'potential' in
            that only consistency of chemical elements and
            number of atoms are considered.

        Raises
        ------
        MissingSublayersError
            If this slab has no sublayers defined.
        """
        if not self.sublayers:
            raise MissingSublayersError
        if self.n_sublayers < 2:
            return []

        # Start from one layer (e.g., the first one) and accumulate
        # potential periods as the offsets of all other layers with
        # same element and same number of atoms. Stop at half the
        # unit cell height, as periodicity cannot go beyond h/2.
        # The same is true for the number of sublayers.
        n_layers = self.n_sublayers
        candidate_periods = []
        ucell_h = self.c_vector[2]
        ref_lay = self.sublayers[0]
        for i, lay in enumerate(self.sublayers[1:], start=1):
            if (abs(lay.cartbotz - ref_lay.cartbotz) > ucell_h / 2 + epsz
                    or i > n_layers / 2):
                break
            if (lay.element == ref_lay.element
                    and lay.n_atoms == ref_lay.n_atoms):
                candidate_periods.append(i)

        # Now make sure that also the layers
        # in between have the same periods
        return [p for p in candidate_periods
                if self._is_candidate_period_acceptable(p, epsz)]

    def _is_candidate_period_acceptable(self, period, epsz):
        """Return whether a sublayer period is an acceptable candidate.

        Parameters
        ----------
        period : int
            Index distance between sublayers to be checked. This method
            assumes that self.sublayers[0] and self.sublayers[period]
            are consistent (same element and number of atoms).
        epsz : float
            Tolerance (Cartesian) on z distances.

        Returns
        -------
        acceptable : bool
            Whether `period` is a suitable sublayer period. This
            means that all the pairs of sublayers of the form
            (self.sublayers[0<i<period], self.sublayers[i+period])
            match one another (same element and number of atoms),
            and have the same z distance as the one between
            self.sublayers[0] and self.sublayers[period].
        """
        n_layers = self.n_sublayers
        ucell_h = self.c_vector[2]

        # Keep track of the distance between layers that
        # we know match: the first one and the one at period.
        # All pairs of layers in between must also be at the
        # same distance.
        z_dist = self.sublayers[period].cartbotz - self.sublayers[0].cartbotz
        layer_pairs = zip(self.sublayers, cycle(self.sublayers, period))
        for i, (lay, lay_plus_period) in enumerate(layer_pairs, start=1):
            if i > n_layers / 2:
                return True
            z_dist_this = lay_plus_period.cartbotz - lay.cartbotz
            z_dist_this %= ucell_h
            if (lay_plus_period.element != lay.element
                    or lay_plus_period.n_atoms != lay.n_atoms
                    or abs(z_dist - z_dist_this) > epsz):
                return False
        return True

    def get_minimal_c_vector(self, eps, epsz=None, z_periodic=True):
        """Return the smallest bulk c vector, if any.

        Notice that this method does not modify this slab. If
        you want this slab to have the shortest c vector, use
        `ensure_minimal_c_vector` instead.

        Parameters
        ----------
        eps : float
            Cartesian tolerance for in-plane distances
        epsz : float or None, optional
            Cartesian tolerance for distances in the z direction,
            i.e., perpendicular to the surface. If not given or
            None, it is taken equal to `eps`.
        z_periodic : bool, optional
            Whether the current slab should be considered periodic
            in the direction of the c vector. This is commonly True
            for all bulk slabs, unless the current c vector is not
            a repeat vector (e.g., while c is being identified using
            detect_bulk). Default is True.

        Returns
        -------
        shortest_c : numpy.ndarray
            The minimal repeat c vector found, in Cartesian
            coordinates with (x, y) in the surface plane, and
            z directed from the solid towards the surface, i.e.,
            opposite to the usual LEED convention. The in-plane
            components are always minimized.

        Raises
        ------
        AlreadyMinimalError
            If a shorter repeat vector could not be found.
        """
        if epsz is None:
            epsz = eps
        periods = self.get_candidate_layer_periods(epsz)
        if not periods:
            raise AlreadyMinimalError(
                f'{type(self).__name__}.get_minimal_c_vector: no '
                'candidate sublayer periods for reducing c vector'
                )

        # Work with a deepcopy, as we will play around with sublayers
        self_copy = copy.deepcopy(self)
        self_copy.update_cartesian_from_fractional()
        self_copy.create_sublayers(epsz)
        n_layers = self_copy.n_sublayers

        # Pick one layer (e.g., the first one) and test translations
        # for all vectors connecting an atom of this layer (e.g., the
        # first one) to all atoms of all other layers.
        ori = self_copy.sublayers[0].cartpos

        # Since the periods are sorted from smaller to larger, we need
        # only keep track of one repeat vector, stopping at the earliest
        _is_symmetric = functools.partial(self_copy.is_translation_symmetric,
                                          eps=eps, z_periodic=z_periodic)
        repeat_c = None
        for period in periods:
            other_layer = self_copy.sublayers[period % n_layers]
            # Vectors from surface to bulk (the way we want them
            # later), but flip them for testing matching, which
            # means translating the slab up towards the surface.
            # This is important when, e.g., we're looking for the
            # c vector in SurfaceSlab.detect_bulk: the upper part
            # of the slab is not bulk.
            test_vecs = (ori - atom.cartpos for atom in other_layer)
            repeat_c = next((v for v in test_vecs if _is_symmetric(-v)), None)
            if repeat_c is not None:
                break

        if repeat_c is None:
            raise AlreadyMinimalError(
                f'{type(self).__name__}.get_minimal_c_vector: none '
                'of the candidate sublayer periods is an actual period'
                )

        # Flip z, because we store it opposed in Atom.cartpos[2], but
        # want repeat_c with the same coordinate system as unit cell
        repeat_c[2] *= -1                                                       # TODO: .cartpos[2] -- Issue #174

        # Minkowski-reduce the repeat vector to have
        # it shortest and as close to z as possible
        leedbase.reduce_c_vector(repeat_c, self.ab_cell.T)
        return repeat_c

    def is_bulk_glide_symmetric(self, symplane, sublayer_period, eps):
        """Return if the slab has a 3D glide plane.

        Parameters
        ----------
        symplane : SymPlane
            The mirror part of the glide plane to test.
        sublayer_period : int
            Number of sublayers, corresponding to the translation
            part of the glide operation to be tested: translations
            are tested for all vectors connecting pairs of sublayers
            whose indices differ by `sublayer_period`. Normally one
            of the values returned by `get_candidate_layer_periods`.
        eps : float
            Tolerance (Cartesian) on atomic positions.

        Returns
        -------
        has_glide : bool
            Whether the slab has a bulk glide plane with the
            given `symplane` and `sublayer_period` translation.
        """
        matrix = symplane.point_operation(n_dim=3)
        return self._is_bulk_transform_symmetric(matrix, sublayer_period, eps)

    def is_bulk_screw_symmetric(self, order, sublayer_period, eps):
        """Return if the slab has a 3D screw axis.

        Parameters
        ----------
        order : int
            The rotation order of the screw axis to test.
        sublayer_period : int
            Number of sublayers, corresponding to the translation
            part of the screw operation to be tested: translations
            are tested for all vectors connecting pairs of sublayers
            whose indices differ by `sublayer_period`. Normally one
            of the values returned by `get_candidate_layer_periods`.
        eps : float
            Tolerance (Cartesian) on atomic positions.

        Returns
        -------
        has_screw : bool
            Whether the slab has a bulk screw with the given
            `order` and `sublayer_period` translation.
        """
        matrix = rotation_matrix_order(order, dim=3)
        return self._is_bulk_transform_symmetric(matrix, sublayer_period, eps)

    # Disabled too-many-locals below: While there are indeed quite
    # a few locals (20/15), refactoring this would require either
    # helper methods that take a lot of arguments or some special
    # class that holds them. That's because all the stuff is computed
    # the very least number of times to speed up execution. Let's
    # prefer execution time and use a lot of comments to help.
    # pylint: disable-next=too-many-locals
    def _is_bulk_transform_symmetric(self, matrix, sublayer_period, eps):
        """Return if this slab is equivalent under a 3D screw/glide.

        Parameters
        ----------
        matrix : numpy.ndarray
            The rotational or mirror part of the symmetry
            transformation to be tested. Shape (3, 3). The
            transformation `matrix` is applied from the left to
            the Cartesian coordinates (taken as column vectors).
            It should represent a point operation that does not
            change the z coordinates, i.e., it is block-diagonal
            with non-trivial elements only in [:2, :2].
        sublayer_period : int
            Number of sublayers to be considered for constructing
            test translation vectors. Should be one of the candidate
            periods returned by `self.get_candidate_layer_periods()`.
            Translational symmetry is tested for all vectors connecting
            atoms between sublayer pairs whose indices differ by
            `sublayer_period`.
        eps : float
            Tolerance (Cartesian) for position equivalence.

        Returns
        -------
        is_symmetric : bool
            True if screw/glide symmetric, else False.
        """
        # Let's use the better version of the unit cell, with
        # unit vector as rows
        ucell = self.ucell.T
        ucell_inv = np.linalg.inv(ucell)
        releps = eps / np.linalg.norm(ucell, axis=1)

        # Get translation vectors to check. Notice that it is enough
        # to pick any single atom from any layer as 'reference', and
        # test translations to all atoms in another layer. To minimize
        # the number of test translations, use pairs of low-occupancy
        # layers. The next lines assume that sublayer_period comes from
        # a call to get_candidate_layer_periods.
        lowocclayer = self.fewest_atoms_sublayer
        ori = matrix.dot(lowocclayer.cartpos)
        ref_layers = itertools.islice(self.sublayers,  # From next one
                                      lowocclayer.num + 1, None)
        # compare_layers: sublayers starting at num+period index,
        # wrapped around. Notice that num+period is the buddy of
        # lowocclayer. We take it out of the iterator while making
        # the translations. After then, compare_layers is aligned
        # with ref_layers.
        compare_layers = cycle(self.sublayers,
                               lowocclayer.num + sublayer_period)
        translations = [at.cartpos - ori for at in next(compare_layers)]

        for ref_layer, compare_layer in zip(ref_layers, compare_layers):
            # Prepare matrix-transformed versions of the Cartesian
            # coordinates of this sublayer, after collapsing to the
            # base unit cell. Translations will be applied later.
            matrix_transformed, _ = collapse(
                np.array([at.cartpos for at in ref_layer]),
                ucell, ucell_inv
                )
            matrix_transformed = matrix_transformed.dot(matrix.T)

            # Get coordinates of the sublayer to compare to, also
            # collapsed to base cell, and including extra atoms
            # for those close to edges and corners
            to_compare, frac_coords = collapse(
                np.array([at.cartpos for at in compare_layer]),
                ucell, ucell_inv
                )
            to_compare, _ = add_edges_and_corners(to_compare,
                                                  frac_coords,
                                                  releps, ucell)
            j = 0
            while j < len(translations):
                transformed_3d = matrix_transformed + translations[j]
                transformed_3d, _ = collapse(transformed_3d, ucell, ucell_inv)
                distances = euclid_distance(transformed_3d, to_compare)
                if any(distances.min(axis=1) > eps):
                    translations.pop(j)
                else:
                    j += 1
            if not translations:
                return False
        return True

    def with_double_thickness(self, new_atoms_start_idx=None):
        """Return a copy of this bulk slab which is twice as thick."""
        double_slab = copy.deepcopy(self)
        c_vec = double_slab.c_vector

        # For atoms that are added, because we use z flipped                    # TODO: .cartpos[2] -- Issue #174
        c_vec_atoms = c_vec.copy()
        c_vec_atoms[2] *= -1
        double_slab.update_cartesian_from_fractional(update_origin=True)

        if new_atoms_start_idx is None:
            # BulkSlab objects tend to have a non-continuous
            # distribution of atom numbers that usually come
            # from the parent SurfaceSlab for which this slab
            # is the bulk. We cannot use the normal n_atoms + 1,
            # as we may end up in a conflict of atom numbers.
            new_atoms_start_idx = max(at.num for at in double_slab) + 1

        # Now decide which way to go, depending on
        # whether there are layers already defined
        if double_slab.layers:
            # pylint: disable=protected-access
            double_slab._add_one_bulk_cell(double_slab.layers, c_vec,
                                           c_vec_atoms, 0, new_atoms_start_idx)
        else:
            for atom in double_slab.atlist.copy():
                atom.duplicate(num=new_atoms_start_idx)
                atom.cartpos += c_vec_atoms
                new_atoms_start_idx += 1
            c_vec[:] *= 2
        double_slab.collapse_cartesian_coordinates(update_origin=True)
        double_slab.sublayers = ()  # They are outdated
        return double_slab
