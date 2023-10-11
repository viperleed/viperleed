# -*- coding: utf-8 -*-
"""Module bulk_slab of viperleed.tleedmlib.classes.slab.

Created on 2023-02-21, originally Jun 13 2019

@author: Florian Kraushofer (@fkraushofer)
@author: Michele Riva (@michele-riva)

Defines the BulkSlab class, a BaseSlab describing a 3D-periodic crystal.
This module was created as part of the refactoring of slab.py.
"""

import copy
import itertools

import numpy as np
from scipy.spatial.distance import cdist as euclid_distance

from viperleed.tleedmlib import leedbase
from viperleed.tleedmlib.base import add_edges_and_corners, collapse
from viperleed.tleedmlib.base import rotation_matrix_order

from .base_slab import BaseSlab
from .slab_utils import _cycle


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
    layers : list of Layer
        List of Layer objects, where each `layer` is a composite
        of sublayers, as in TensErLEED
    sublayers : list of SubLayer
        List of SubLayer objects, each containing atoms of equal
        element and Z coordinate
    sitelist : list of Sitetype
        List of distinct sites as Sitetype, storing information
        on vibration and concentration
    ucell_mod : list of tuples (str, numpy.ndarray)
        Stored modifications made to the unit cell; each is a tuple
        of (type, array), where type is 'lmul', 'rmul', or 'add'
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
    symbaseslab : Slab or None                                                  # TODO: do we need one??
        Slab with the smallest in-plane unit-cell area that shows
        the full symmetry of the slab.
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
        bulkc = self.ucell.T[2].copy()
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
        """
        if self.n_sublayers < 2:
            return []

        # Start from one layer (e.g., the first one) and accumulate
        # potential periods as the offsets of all other layers with
        # same element and same number of atoms. Stop at half the
        # unit cell height, as periodicity cannot go beyond h/2                 # TODO: @fkrausofer: probably we should also stop at floor(len(sublayers) / 2), as we cannot have a sub-period longer than half the number of sublayers?
        candidate_periods = []
        ucell_h = self.ucell[2, 2]
        ref_lay = self.sublayers[0]
        for i, lay in enumerate(self.sublayers[1:], start=1):
            if abs(lay.cartbotz - ref_lay.cartbotz) > ucell_h / 2 + epsz:
                break
            if (lay.element == ref_lay.element
                    and lay.n_atoms == ref_lay.n_atoms):
                candidate_periods.append(i)

        # Now make sure that also the layers
        # in between have the same periods
        n_layers = self.n_sublayers
        i = 0
        while i < len(candidate_periods):
            period_is_ok = True
            period = candidate_periods[i]

            # Keep track of the distance between layers that
            # we know match: the first one and the one at period.
            # All pairs of layers in between must also be at the
            # same distance.
            z_dist = self.sublayers[period].cartbotz - ref_lay.cartbotz
            layer_pairs = zip(self.sublayers, _cycle(self.sublayers, period))
            for j, (lay, lay_plus_period) in enumerate(layer_pairs, start=1):
                if j > n_layers / 2:
                    break
                z_dist_this = lay_plus_period.cartbotz - lay.cartbotz
                z_dist_this %= ucell_h
                if (lay_plus_period.element != lay.element
                        or lay_plus_period.n_atoms != lay.n_atoms
                        or abs(z_dist - z_dist_this) > epsz):
                    period_is_ok = False
                    break
            if period_is_ok:  # Period passed the tests
                i += 1
            else:
                candidate_periods.pop(i)
        return candidate_periods

    def getMinC(self, rp, z_periodic=True):
        """Checks whether there is a vector c with a smaller length than
        the current one. If so, returns the minimized vector, else returns
        None."""
        eps = rp.SYMMETRY_EPS
        pcands = self.get_candidate_layer_periods(eps)
        if len(pcands) == 0:
            return None
        ts = copy.deepcopy(self)
        ts.update_cartesian_from_fractional()
        ts.create_sublayers(eps)
        baseLayer = ts.sublayers[0]
        baseInd = ts.sublayers.index(baseLayer)
        nl = ts.n_sublayers
        ori = baseLayer.cartpos  # compare displacements from here
        repeatC = None
        for per in pcands:
            ind = (baseInd + per) % nl
            for at in ts.sublayers[ind]:
                v = ori - at.cartpos
                if ts.is_translation_symmetric(v, eps, z_periodic=z_periodic):
                    repeatC = at.cartpos - ori
                    break
            if repeatC is not None:
                break
        if repeatC is None:
            return None
        # optimize C vector to be close to Z, if possible
        cFracBase = np.dot(np.linalg.inv(ts.ab_cell), repeatC[:2]) % 1.0
        newC = np.append(np.dot(ts.ab_cell, cFracBase), -repeatC[2])
        for (i, j) in [(0, -1), (-1, 0), (-1, -1)]:
            v = np.dot(ts.ab_cell, cFracBase + np.array([i, j]))
            if (np.linalg.norm(np.append(v, -repeatC[2]))
                    < np.linalg.norm(newC)):
                newC[:2] = v
        return newC

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
        compare_layers = _cycle(self.sublayers,
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
        *_, c_vec = double_slab.ucell.T

        # For atoms that are added, because we use z flipped                    # TODO: .cartpos[2]
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
            c_vec *= 2
        double_slab.collapse_cartesian_coordinates(update_origin=True)
        double_slab.sublayers.clear()  # They are outdated
        return double_slab
