# -*- coding: utf-8 -*-
"""Module base_slab of viperleed.tleedmlib.classes.slab.

Created on 2023-02-21, originally Jun 13 2019

@author: Florian Kraushofer (@fkraushofer)
@author: Michele Riva (@michele-riva)

Defines the BaseSlab class, useful to describe collections of atoms in
crystalline form. This is the abstract base class at the basis of both
BulkSlab and SurfaceSlab classes (for 3D- and 2D-periodic systems,
respectively), and contains generic functionality that is common to
both. This module contains refactored and modified functionality that
used to be contained in the original slab module by F. Kraushofer.
"""

from abc import ABC, abstractmethod
from collections import Counter
import copy
import itertools
import logging
from numbers import Real
from operator import attrgetter, itemgetter
import re

import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist as euclid_distance

from viperleed.tleedmlib import leedbase
from viperleed.tleedmlib.base import add_edges_and_corners, angle, collapse
from viperleed.tleedmlib.base import collapse_fractional, pairwise
from viperleed.tleedmlib.base import rotation_matrix, rotation_matrix_order
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.layer import Layer
from viperleed.tleedmlib.classes.sitetype import Sitetype

from .slab_errors import InvalidUnitCellError
from .slab_errors import NeedsLayersError, NeedsSublayersError, SlabError
from .slab_utils import _z_distance


_LOGGER = logging.getLogger('tleedm.slab')


# TODO: would it make sense to have cartpos[2] go the other way? Watch out for the
#       "reference" z in LEED that should always be the z of the topmost atom AFTER REFCALC
#       Could store the top_z, then have a .leed_pos attribute that returns top_z - cartpos[2]
# TODO: too-many-instance-attributes
# TODO: layer coordinates may not be up to date after we do update_origin
# TODO: a huge fraction of the time spent when dealing with slab symmetry
#       operations is actually taken up by calls to deepcopy. It would help
#       quite a bit to have a .light_copy(which_purpose) method that makes
#       appropriate copies of selected parts of the slab. Most of the times,
#       we only need ucell and some lightweight version of atlist (maybe as
#       a tuple of lightweight Atoms that do not hold references to this slab).
class BaseSlab(ABC):
    """An abstract base class representing a solid.

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
    atlist : list of Atom
        List of all atoms in the slab.
    layers : list of Layer
        List of Layer objects, where each `layer` is a composite
        of sublayers, as in TensErLEED
    sublayers : list of Layer
        List of Layer objects, each containing atoms of equal
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
    symbaseslab : Slab or None
        Slab with the smallest in-plane unit-cell area that shows
        the full symmetry of the slab.
    """

    def __init__(self):
        """Initialize instance."""
        self.ucell = np.array([])                                               # base
        self.poscar_scaling = 1.                                                # base
        self.chemelem = set()                                                   # base
        self.n_per_elem = {}                                                    # base
        self.atlist = []                                                        # base
        self.layers = []                                                        # base
        self.sublayers = []                                                     # base
        self.sitelist = []                                                      # base
        self.ucell_mod = []                                                     # base
        self.ucell_ori = np.array([])                                           # base
        self.topat_ori_z = None                                                 # base (non-bulk after we fix the cartpos[2] flip)
        self.celltype = 'unknown'                                               # base
        self.planegroup = 'unknown'                                             # base
        self.foundplanegroup = 'unknown'                                        # base
        self.orisymplane = None                                                 # base
        self.linklists = []                                                     # base?
        self.symbaseslab = None                                                 # non-bulk?
        self.bulkslab = None  # Deleted in BulkSlab.__init__

        # Remember the last value of the ELEMENT_MIX parameter that
        # was applied. Prevents repeated applications
        self._last_element_mix = None                                           # base?

    def __contains__(self, item):
        """Return whether an atom, layer, site, ... is in this slab."""
        if isinstance(item, Atom):
            return item in self.atlist
        if isinstance(item, Layer):
            return item in self.sublayers or item in self.layers
        if isinstance(item, Sitetype):
            return item in self.sitelist
        if isinstance(item, str):   # Only chemical element?
            return item in self.chemelem
        return False

    def __iter__(self):
        """Return an iterator of Atoms in self."""
        return iter(self.atlist)

    def __reversed__(self):
        """Return a reversed iterator over Atoms in self."""
        return reversed(self.atlist)

    @property
    def ab_cell(self):
        """Return the 2D portion of the unit cell."""
        try:
            return self.ucell[:2, :2]
        except IndexError:  # Uninitialized
            raise InvalidUnitCellError(
                f'{type(self).__name__} has no unit cell defined.'
                ) from None

    @property
    def angle_between_ucell_and_coord_sys(self):
        """Return angle between first unit-cell vector and coordinate system.

        Returns
        -------
        angle : float
            Angle between first slab unit cell vector and Cartesian
            coordinate system in degrees.

        Raises
        ------
        InvalidUnitCellError
            If this property is accessed before there is a unit cell
            defined.
        """
        a_vec, *_ = self.ab_cell.T
        # NB: arctan2 requires (y, x) order
        return np.degrees(np.arctan2(a_vec[1], a_vec[0]))

    @property
    def bulk_layers(self):
        """Return the layers of self that are bulk."""
        return [lay for lay in self.layers if lay.isBulk]

    @property
    def elements(self):
        """Return a tuple of elements in this slab, as originally read."""
        return tuple(self.n_per_elem.keys())

    @property
    def fewest_atoms_sublayer(self):
        """Return the sublayer with fewest atoms."""
        if not self.sublayers:
            raise NeedsSublayersError(f'{type(self).__name__} has '
                                      'no sublayers defined')
        return min(self.sublayers, key=lambda lay: len(lay.atlist))

    @property
    @abstractmethod
    def is_bulk(self):
        """Return whether this is a bulk slab."""
        return False

    @property
    def n_atoms(self):
        """Return the number of atoms in this slab."""
        return len(self.atlist)

    @property
    def n_layers(self):
        """Return the number of composite layers of this slab."""
        return len(self.layers)

    @property
    def n_sublayers(self):
        """Return the number of sublayers in this slab."""
        return len(self.sublayers)

    #                                                                           TODO: remove. Used only once. Also confusing because it's only in-plane
    @property
    def reciprocal_vectors(self):
        """Returns the reciprocal lattice vectors as an array.

        Reciprocal vectors are defined by
        $a_i \dot$ b_j = 2\pi\delta_{ij}$.
        We need to transpose here again, because of swap row <-> col
        when going between real and reciprocal space.
        TensErLEED always does this calculation explicitly and
        normalizes to area, but here the inverse already contains a
        factor of 1/det.

        Returns
        -------
        np.ndarray, shape=(2, 2)
            Array of *reciprocal* lattice vectors *as rows*.
        """
        return 2*np.pi*np.linalg.inv(self.ab_cell.T).T

    @property
    def smallest_interlayer_spacing(self):
        """Return the smallest z gap between two adjacent layers.

        Returns
        -------
        min_dist : float
            The smallest of the z distances between adjacent layers.
            Distances are calculated between the topmost atom of the
            lower layer and the bottommost one of the higher. Zero if
            there is only one layer.

        Raises
        ------
        NeedsLayersError
            If no layers are available
        """
        if not self.layers:
            raise NeedsLayersError

        if self.n_layers == 1:
            return 0.                                                           # TODO: I don't think it's right that it is zero if there's only one layer. Think about it.
        # self.update_cartesian_from_fractional()                               # TODO: I don't think this is needed. It does not update anything for layers; only makes sense if we also .update_layer_coordinates. Would make more sense if layer.cartori and .cartbotz were @property

        # Recall that z increases moving deeper into the solid
        return min(lay_below.cartori[2] - lay_above.cartbotz                    # TODO: change when flipping .cartpos[2]
                   for lay_above, lay_below in pairwise(self.layers))

    @classmethod
    def from_slab(cls, other):
        """Return a `cls` instance with attributes deep-copied from `other`."""
        if not isinstance(other, BaseSlab):
            raise TypeError(f'{cls.__name__}.from_slab: other is not a slab.')
        if type(other) is cls:
            return copy.deepcopy(other)

        instance = cls()
        memo = {id(other): instance}
        for attr, value in other.__dict__.items():
            if not hasattr(instance, attr):
                # Skip attributes that do not belong to cls instances
                continue
            setattr(instance, attr, copy.deepcopy(value, memo))
        return instance

    def addBulkLayers(self, rp, n=1):
        """Returns a copy of the slab with n bulk units appended at the
        bottom, and a list of the new atoms that were added."""
        ts = copy.deepcopy(self) # temporary slab
        newbulkats = []
        duplicated = []
        zdiff = 0.
        for _ in range(n):
            blayers = ts.bulk_layers
            if isinstance(rp.BULK_REPEAT, np.ndarray):
                bulkc = np.copy(rp.BULK_REPEAT)
                if bulkc[2] < 0:
                    # perhaps vector given from surface to bulk instead of reverse...
                    bulkc = -bulkc
            else:
                cvec = ts.ucell[:, 2]
                if zdiff == 0. and rp.BULK_REPEAT is None:
                    # assume that interlayer vector from bottom non-bulk to top
                    # bulk layer is the same as between bulk units
                    zdiff = (blayers[-1].cartbotz
                             - ts.layers[blayers[0].num-1].cartbotz)
                elif zdiff == 0. and isinstance(rp.BULK_REPEAT,
                                                (float, np.floating)):
                    zdiff = rp.BULK_REPEAT
                bulkc = cvec * zdiff / cvec[2]
            ts.update_cartesian_from_fractional()
            cfact = (ts.ucell[2, 2] + abs(bulkc[2])) / ts.ucell[2, 2]
            ts.ucell[:, 2] = ts.ucell[:, 2] * cfact
            bulkc[2] = -bulkc[2]
            original_atoms = ts.atlist[:] # all atoms before adding layers

            # split bulkc into parts parallel and orthogonal to unit cell c
            # this allows to keep the same ucell and shift correctly the new bulk layers
            c_direction = ts.ucell[:, 2] / np.dot(ts.ucell[:, 2], ts.ucell[:, 2])
            bulkc_project_to_c = np.dot(bulkc, ts.ucell[:, 2]) * c_direction
            bulkc_perp_to_c = bulkc - bulkc_project_to_c
            added_this_loop = []
            for at in original_atoms:
                if at.layer.isBulk and at not in duplicated:
                    new_atom = at.duplicate()
                    newbulkats.append(new_atom)
                    duplicated.append(at)
                    added_this_loop.append(new_atom)
                    new_atom.oriN = ts.n_atoms

                # old atoms get shifted up along ucell c
                at.cartpos += bulkc_project_to_c
            for at in added_this_loop:
                # new atoms get shifted perpendicular to ucell c
                at.cartpos -= bulkc_perp_to_c
            # TODO: could be done outside loop?
            ts.collapse_cartesian_coordinates(update_origin=True)
            ts.sort_original()
        return ts, newbulkats

    def check_a_b_in_plane(self):
        """Raise InvalidUnitCellError if a, b have out-of-plane components."""
        if any(self.ucell[2, :2]):
            _err = ('Unit cell a and b vectors must not '
                    'have an out-of-surface (Z) component!')
            _LOGGER.error(_err)
            raise InvalidUnitCellError(_err)

    def clear_symmetry_and_ucell_history(self):                                 # base or only surface?
        """Set all symmetry information back to default values.

        This method also completely erases the history of unit-cell
        modifications, and sets the current unit cell as the 'original'
        one. This means that, if the unit cell was modified prior to
        a call to this method, the original unit cell **cannot** be
        recovered by a call to revertUnitCell.

        Returns
        -------
        None.
        """
        self.ucell_mod = []
        self.ucell_ori = self.ucell.copy()
        self.celltype = 'unknown'
        self.foundplanegroup = self.planegroup = 'unknown'
        self.orisymplane = None

    def createLayers(self, rparams, bulk_cuts=[]):
        """Creates a list of Layer objects based on the N_BULK_LAYERS and
        LAYER_CUTS parameters in rparams. If layers were already defined,
        overwrite. The bulk_cuts kwarg allows specifically inserting
        automatically detected bulk layer cuts. Returns the cuts as a sorted
        list of floats."""
        # first interpret LAYER_CUTS parameter - can be a list of strings
        self.check_a_b_in_plane()

        if self.is_bulk:
            bulk_cuts = ()

        ct = []
        rgx = re.compile(r'\s*(dz|dc)\s*\(\s*(?P<cutoff>[0-9.]+)\s*\)')
        al = self.atlist[:]
        al.sort(key=lambda atom: atom.pos[2])
        for (i, s) in enumerate(rparams.LAYER_CUTS):
            if type(s) == float:
                ct.append(s)
                continue
            s = s.lower()
            if 'dz' in s or 'dc' in s:
                m = rgx.match(s)
                if not m:
                    _LOGGER.warning('Error parsing part of LAYER_CUTS: ' + s)
                    continue
                cutoff = float(m.group('cutoff'))
                lowbound = 0.
                if bulk_cuts:
                    lowbound = max(bulk_cuts)
                highbound = 1.
                val = None
                if (i > 1) and (rparams.LAYER_CUTS[i-1] in ['<', '>']):
                    try:
                        val = float(rparams.LAYER_CUTS[i-2])
                    except ValueError:
                        _LOGGER.warning('LAYER_CUTS: Error parsing left-hand '
                                        'boundary for ' + s)
                    if val is not None:
                        if rparams.LAYER_CUTS[i-1] == '<':
                            lowbound = val
                        else:
                            highbound = val
                if i < len(rparams.LAYER_CUTS) - 2 and (rparams.LAYER_CUTS[i+1]
                                                        in ['<', '>']):
                    try:
                        val = float(rparams.LAYER_CUTS[i+2])
                    except ValueError:
                        _LOGGER.warning('LAYER_CUTS: Error parsing right-hand '
                                        'boundary for ' + s)
                    if val is not None:
                        if rparams.LAYER_CUTS[i+1] == '>':
                            lowbound = val
                        else:
                            highbound = val
                if 'dc' in s:
                    cutoff *= (self.ucell[2, 2]
                               / np.linalg.norm(self.ucell[:, 2]))
                for i in range(1, len(al)):
                    if ((abs(al[i].cartpos[2]-al[i-1].cartpos[2]) > cutoff)
                            and al[i].pos[2] > lowbound
                            and al[i].pos[2] < highbound
                            and al[i-1].pos[2] > lowbound
                            and al[i-1].pos[2] < highbound):
                        ct.append(abs((al[i].pos[2]+al[i-1].pos[2])/2))
            elif s not in ['<', '>']:
                try:
                    ct.append(float(s))
                except ValueError:
                    _LOGGER.warning('LAYER_CUTS: Could not parse value: ' + s)
                    continue
        if bulk_cuts:
            ct = [v for v in ct if v > max(bulk_cuts) + 1e-6] + bulk_cuts
        ct.sort()
        self.layers = []
        tmplist = self.atlist[:]
        self.sort_by_z()
        laynum = 0
        b = True if rparams.N_BULK_LAYERS > 0 else False
        newlayer = Layer(self, 0, b)
        self.layers.append(newlayer)
        for atom in self:
            # only check for new layer if we're not in the top layer already
            if laynum < len(ct):
                if atom.pos[2] > ct[laynum]:
                    # if atom is higher than the next cutoff, make a new layer
                    laynum += 1
                    b = True if rparams.N_BULK_LAYERS > laynum else False
                    newlayer = Layer(self, laynum, b)
                    self.layers.append(newlayer)
                    check = True    # check for empty layer
                    while check:
                        if laynum >= len(ct):
                            check = False
                        elif atom.pos[2] <= ct[laynum]:
                            check = False
                        else:
                            laynum += 1
                            b = (True if rparams.N_BULK_LAYERS > laynum
                                 else False)
                            newlayer = Layer(self, laynum, b)
                            self.layers.append(newlayer)
            atom.layer = newlayer
            newlayer.atlist.append(atom)
        dl = []
        for layer in self.layers:
            if not layer.atlist:
                _LOGGER.warning('A layer containing no atoms was found. Layer '
                                'will be deleted. Check LAYER_CUTS parameter.')
                rparams.setHaltingLevel(2)
                dl.append(layer)
        for layer in dl:
            if layer.isBulk:
                self.layers[layer.num+1].isBulk = True
            self.layers.remove(layer)
            del layer
        self.layers.reverse()
        for i, layer in enumerate(self.layers):
            layer.getLayerPos()
            layer.num = i
        self.atlist = tmplist
        return ct

    def _get_sublayers_for_el(self, element, eps):
        """Yield Layer objects for atoms of `element` within `eps`."""
        sublists = [[a for a in self if a.el == element]]
        # First, split at points where two atoms are more than eps apart
        i = 0
        while i < len(sublists):
            atoms = sublists[i]
            for j, atom_pair in enumerate(pairwise(atoms)):
                if _z_distance(*atom_pair) > eps:
                    sublists[i:i+1] = atoms[:j+1], atoms[j+1:]
                    break
            i += 1

        # Now go through again and split sublayers at greatest
        # interlayer distance, if they are too thick overall
        i = 0
        while i < len(sublists):
            atoms = sublists[i]
            if len(atoms) < 2 or _z_distance(atoms[0], atoms[-1]) <= eps:
                i += 1
                continue
            distances = ((j, _z_distance(*atoms))
                         for j, atoms in enumerate(pairwise(atoms), start=1))
            maxdistindex, _ = max(distances, key=itemgetter(1))
            sublists[i:i+1] = atoms[:maxdistindex], atoms[maxdistindex:]

        # Finally create sublayers based on sublists
        for atoms in sublists:
            new_layer = Layer(self, 0)
            new_layer.atlist = atoms
            new_layer.cartbotz = atoms[0].cartpos[2]
            yield new_layer

    def create_sublayers(self, eps=0.001):
        """Assign the atoms in the slab to sublayers.

        Each sublayer contains only atoms with the same chemical
        element, and such that all atoms in the sublayer are closer
        than `eps` to one another in the z direction.

        After calling this method, `slab.atlist` is sorted along
        the z Cartesian coordinate (bottommost atom first), and the
        `slab.sublayers` attribute is up to date. Calling this method
        more than once erases previous sublayer lists.

        The sublayers created are stored in the following order:
        - Sublayers towards the top of the slab come earlier
        - Sublayers with the same z coordinate (within eps) are
          next to each other, sorted by alphabetical element order

        Parameters
        ----------
        eps : float, optional
            Limit value for z difference to assign atoms to different
            sublayers. Default is 0.001.

        Returns
        -------
        None.
        """
        self.check_a_b_in_plane()
        self.sort_by_z()
        subl = []
        for element in self.elements:                                           # TODO: should we complain if there's no elements?
            subl.extend(self._get_sublayers_for_el(element, eps))

        # subl is sorted element-first. Re-sort it as in the __doc__
        self.sublayers = []

        # Work with subl sorted from bottom to top, and pop the last
        # element each time (i.e., the topmost layer to be processed)
        subl.sort(key=attrgetter('cartbotz'), reverse=True)
        while subl:
            same_z = [subl.pop()]  # accumulate sublayers with same z
            this_z = same_z[0].cartbotz
            while subl:
                if abs(subl[-1].cartbotz - this_z) < eps:
                    same_z.append(subl.pop())
                else:
                    break
            # Finally re-sort by element
            same_z.sort(key=lambda lay: lay.atlist[0].el)                       # TODO: is this even necessary? .sort should be stable, so element order should be preserved
            self.sublayers.extend(same_z)
        for i, layer in enumerate(self.sublayers):
            layer.num = i

    def collapse_cartesian_coordinates(self, update_origin=False):
        """Ensure all atoms are inside the unit cell.

        This method updates both Cartesian and fractional coordinates,
        using the current Cartesian coordinates as starting point.

        Parameters
        ----------
        update_origin : bool, optional
            Whether the z component of the Cartesian origin should
            also be updated. Default is False.

        Returns
        -------
        None.
        """
        self.update_fractional_from_cartesian()
        self.collapse_fractional_coordinates()
        self.update_cartesian_from_fractional(update_origin=update_origin)

    def collapse_fractional_coordinates(self):                                  # TODO: would be nicer to have the public method update everything correctly, and rather make this a private one.  USED ONLY IN SLABs and in iosearch
        """Ensure all atoms are inside the unit cell.

        Atoms whose fractional positions are outside the unit cell
        are "back-folded" inside.

        Notice that calling this method DOES NOT UPDATE the Cartesian
        coordinates of the atoms. It must be preceded or followed by
        a call to `update_cartesian_from_fractional` to ensure that
        fractional and Cartesian coordinates are up to date.

        Returns
        -------
        None.

        See Also
        --------
        collapse_cartesian_coordinates
        """
        for atom in self:
            collapse_fractional(atom.pos, in_place=True)

    def fullUpdate(self, rpars):
        """readPOSCAR initializes the slab with information from POSCAR;
        fullUpdate re-initializes the atom list, then uses the information
        from the parameters file to create layers, calculate cartesian
        coordinates (absolute and per layer), and to update elements and
        sites."""
        self.collapse_fractional_coordinates()
        self.update_cartesian_from_fractional()
        if not self.layers:
            self.createLayers(rpars)
        self._update_chem_elements(rpars)
        if not self.sitelist:
            self.initSites(rpars)
        if rpars.fileLoaded['VIBROCC']:
            for at in self:
                at.initDisp()

    @abstractmethod
    def getBulkRepeat(self, rp):
        """Based on a pre-existing definition of the bulk, tries to identify
        a repeat vector for which the bulk matches the slab above. Returns that
        vector in cartesian coordinates, or None if no match is found."""

    def getMinUnitCell(self, rp, warn_convention=False):
        """Check if there is a 2D unit cell smaller than the current one.

        Parameters
        ----------
        rp : RunParams
            The current parameters. The only attributes
            used are SYMMETRY_EPS and SYMMETRY_EPS_Z.
        warn_convention : bool, optional
            If True, warnings are added to the current
            logger in case making the reduced unit cell
            stick to the conventions would result in a
            sub-optimal superlattice matrix. Default is
            False.

        Returns
        -------
        can_be_reduced : bool
            True if there is a smaller 2D unit cell. The
            unit cell is considered minimizable if there
            is a mincell with area smaller than the one
            of the current cell. A lower limit for the
            area of mincell is taken as 1 A**2.
        mincell : np.ndarray
            The minimal 2D unit cell, if it can be reduced,
            otherwise the current one. Notice that mincell is
            such that (a, b) = mincell, i.e., it is transposed
            with respect to self.ucell.
        """
        # TODO: write a testcase for the reduction of POSCAR Sb on Si(111)
        eps = rp.SYMMETRY_EPS
        epsz = rp.SYMMETRY_EPS_Z
        abst = self.ab_cell.T

        # Create a test slab: C projected to Z
        ts = copy.deepcopy(self)
        ts.project_c_to_z()
        ts.sort_by_z()
        ts.create_sublayers(epsz)

        # Use the lowest-occupancy sublayer (the one
        # with fewer atoms of the same site type)
        lowocclayer = ts.fewest_atoms_sublayer
        n_atoms = len(lowocclayer.atlist)
        if n_atoms < 2:
            # Cannot be smaller if there's only 1 atom
            return False, abst

        # Create a list of candidate translation vectors, selecting
        # only those for which the slab is translation symmetric
        plist = [at.cartpos[0:2] for at in lowocclayer.atlist]
        vlist = ((p1 - p2) for (p1, p2) in itertools.combinations(plist, 2))
        tvecs = [v for v in vlist if ts.isTranslationSymmetric(v, eps)]
        if not tvecs:
            return False, abst

        # Now try to reduce the cell: test whether we can use a pair of
        # vectors from [a, b, *tvecs] to make the cell smaller. Keep in
        # mind that with n_atoms, we cannot reduce the area by more than
        # a factor 1/n_atoms (which would give 1 atom per mincell).
        mincell = abst.copy()
        mincell_area = abs(np.linalg.det(mincell))
        smaller = False
        smallest_area = mincell_area / n_atoms
        for vec in tvecs:
            # Try first replacing the current second unit vector
            tcell = np.array([mincell[0], vec])
            tcell_area = abs(np.linalg.det(tcell))
            if (tcell_area >= smallest_area - eps**2
                    and tcell_area < mincell_area - eps**2):
                mincell = tcell
                mincell_area = tcell_area
                smaller = True
                continue

            # Try replacing the current first unit vector instead
            tcell = np.array([mincell[1], vec])
            tcell_area = abs(np.linalg.det(tcell))
            if (tcell_area >= smallest_area - eps**2
                    and tcell_area < mincell_area - eps**2):
                mincell = tcell
                mincell_area = tcell_area
                smaller = True

        if not smaller:
            return False, abst

        # Use Minkowski reduction to make mincell high symmetry
        mincell, _, _ = leedbase.reduceUnitCell(mincell)

        # Cosmetic corrections
        if abs(mincell[0, 0]) < eps and abs(mincell[1, 1]) < eps:
            # Swap a and b when matrix is off-diagonal
            mincell[[0, 1]] = mincell[[1, 0]]
        if abs(mincell[1, 0]) < eps and abs(mincell[0, 1]) < eps:
            # If matrix is diagonal, make elements positive
            mincell = abs(mincell)
        # By convention, make the shorter vector the first one
        if np.linalg.norm(mincell[0]) > np.linalg.norm(mincell[1]) + eps:
            if abs(mincell[1, 0]) < eps and abs(mincell[0, 1]) < eps:
                # if matrix is diagonal, DO NOT make it off-diagonal
                if warn_convention:
                    _LOGGER.warning(
                        'The unit cell orientation does not follow '
                        'standard convention: to keep SUPERLATTICE matrix '
                        'diagonal, the first bulk vector must be larger '
                        'than the second. Consider swapping the unit cell '
                        'vectors.'
                        )
            else:
                mincell = np.dot([[0, 1], [-1, 0]], mincell)
        # Finally, make sure it's right-handed
        if angle(mincell[0], mincell[1]) < 0:
            mincell = np.dot([[1, 0], [0, -1]], mincell)
        return True, mincell

    def get_nearest_neigbours(self):
        """Return the nearest-neighbour distance for all atoms.

        Returns
        -------
        nn_distances : dict
            Form is {`atom`: `nn_distance`}, with `nn_distance`
            the shortest distance between `atom` and all its
            neighbours, taking into account the periodicity of
            the slab in the direction parallel to the surface.
        """
        # Work with a supercell to account for atoms at edges/corners.
        # Decide based on unit vector lengths how many times we want
        # to repeat in each direction: the more the unit vectors are
        # mismatched the more repetitions.                                      # TODO: there may be a better way taking angle into account as well.
        *ab_cell, _ = self.ucell.T  # ab_cell is 2x3
        ab_norms = np.linalg.norm(ab_cell, axis=1)
        repeats = np.ceil(ab_norms.max() / ab_norms)  # '+' direction
        repeats = 2*repeats + 1           # '+', '-', and centre cell
        supercell = self.make_supercell(np.diag(repeats))

        # For NN query use KDTree from scipy.spatial
        atom_coords_supercell = [atom.cartpos for atom in supercell]
        tree = KDTree(atom_coords_supercell)

        # Coordinates of all atoms in the centre cell:
        center_cell_offset = (repeats + 1) / 2  # Indices only
        center_cell_offset = center_cell_offset.dot(ab_cell)
        atom_coords_center = [atom.cartpos + center_cell_offset
                              for atom in self]

        # Get all distances. Use k==2 as k==1 would just return
        # a bunch of zeros since all atom_coords_center are a
        # subset of atom_coords_supercell: their first nearest
        # neighbour is themselves
        distances, _ =  tree.query(atom_coords_center, k=2)
        distances = distances[:, 1]  # First column is zeros
        return dict(zip(self, distances))

    def initSites(self, rp):
        """Goes through the atom list and supplies them with appropriate
        SiteType objects, based on the SITE_DEF parameters from the supplied
        Rparams."""
        atlist = self.atlist[:]     # copy to not have any permanent changes
        atlist.sort(key=lambda atom: atom.oriN)
        sl = []
        for el in rp.SITE_DEF:
            for sitename in rp.SITE_DEF[el]:
                newsite = Sitetype(el, sitename)
                sl.append(newsite)
                for i in rp.SITE_DEF[el][sitename]:
                    try:
                        if atlist[i-1].el != el:
                            _LOGGER.warning(
                                'SITE_DEF tries to assign atom number '
                                + str(i) + ' as ' + el + ', but POSCAR has it '
                                'as '+atlist[i-1].el+'. Atom will be skipped '
                                'and left as default site type!')
                            rp.setHaltingLevel(1)
                        else:
                            atlist[i-1].site = newsite
                    except IndexError:
                        _LOGGER.error('SITE_DEF: atom number out of bounds.')
                        raise
        for el in self.elements:
            newsite = Sitetype(el, 'def')
            found = False
            for at in atlist:
                if at.el == el and at.site is None:
                    at.site = newsite
                    found = True
            if found:
                sl.append(newsite)
        for site in [s for s in sl if s.el in rp.ELEMENT_MIX]:
            site.mixedEls = rp.ELEMENT_MIX[site.el][:]
        self.sitelist = sl

    def isEquivalent(self, slab, eps=0.001):
        """Compares the slab to another slab, returns True if all atom cartpos
        match (with at least one other atom, if there are duplicates), False
        if not. Both slabs are copied and collapsed to the (0,0) cell
        before."""
        slab1 = copy.deepcopy(self)
        slab2 = copy.deepcopy(slab)
        slab1.collapse_cartesian_coordinates()
        slab2.collapse_cartesian_coordinates()
        # reorder sublayers by Z to then compare by index
        slab1.sublayers.sort(key=lambda sl: sl.cartbotz)
        slab2.sublayers.sort(key=lambda sl: sl.cartbotz)
        ab = self.ab_cell
        for (i, sl) in enumerate(slab1.sublayers):
            if (len(sl.atlist) != len(slab2.sublayers[i].atlist)
                    or abs(sl.cartbotz-slab2.sublayers[i].cartbotz) > eps
                    or sl.atlist[0].el != slab2.sublayers[i].atlist[0].el):
                return False
            for at1 in sl.atlist:
                complist = [at1.cartpos[0:2]]
                # if we're close to an edge or corner, also check translations
                for j in range(0, 2):
                    releps = eps / np.linalg.norm(ab[:, j])
                    if abs(at1.pos[j]) < releps:
                        complist.append(at1.cartpos[:2] + ab[:, j])
                    if abs(at1.pos[j]-1) < releps:
                        complist.append(at1.cartpos[:2] - ab[:, j])
                if len(complist) == 3:
                    # coner - add the diagonally opposed one
                    complist.append(complist[1] + complist[2] - complist[0])
                found = False
                for at2 in slab2.sublayers[i].atlist:
                    for p in complist:
                        if np.linalg.norm(p-at2.cartpos[0:2]) < eps:
                            found = True
                            break
                    if found:
                        break
                if not found:
                    return False
        return True

    def make_supercell(self, transform):                                        # TODO: surface only?
        """Return a copy of the slab replicated according to `transform`.

        The 'inverse' (i.e., leading to a size reduction) of this
        operation can be obtained by calling `makeSymBaseSlab` with
        the same `transform`.

        Parameters
        ----------
        transform : numpy.ndarray
            Shape (2, 2). All elements must be (close to) integers.
            The transposed unit cell is transformed by multiplication
            with `transform` on the left. Notice that this is NOT the
            usual way of doing it, e.g., as VESTA does. The equivalent
            matrix in VESTA is `transform.T`.

        Returns
        -------
        super_slab : Slab
            Copy of this slab replicated according to `transform`.

        Raises
        ------
        NonIntegerMatrixError
            If `transform` has non-integer entries.
        SingularMatrixError
            If `transform` is singular.
        ValueError
            If `transform` has unacceptable shape.
        """
        super_slab = copy.deepcopy(self)
        try:
            repeats = leedbase.get_superlattice_repetitions(transform)
        except ValueError as exc:
            raise type(exc)(
                f'{type(self).__name__}.make_supercell: {exc}'
                ) from None

        if any(r > 1 for r in repeats):  # Need to duplicate atoms
            for atom in super_slab.atlist.copy():
                for i, j in itertools.product(*repeats):
                    # pylint: disable=compare-to-zero
                    if i == j == 0:  # Skip already existing atom
                        continue
                    # Duplicate the atom, and implicitly store it
                    # in super_slab. Then change its fractional pos
                    duplicate_atom = atom.duplicate()
                    duplicate_atom.pos[:1] += (i, j)

        # We now may have added new atoms that will have fractional
        # coordinates outside the original cell. Make sure we have
        # the correct Cartesian coordinates before transforming the
        # unit cell.
        super_slab.update_cartesian_from_fractional()                           # TODO: was update_origin_True, but should not be needed as we have only added atoms by moving them in-plane
        super_slab.ab_cell[:] = super_slab.ab_cell.dot(transform.T)

        # Starting from the stored Cartesian coordinates, collapse
        # all atoms (incl. fractional coordinates) to the new cell
        super_slab.collapse_cartesian_coordinates()                             # TODO: was update_origin_True, same comment as above
        return super_slab

    def makeSymBaseSlab(self, rp, transform=None):
        """Copies self to create a symmetry base slab by collapsing to the
        cell defined by rp.SYMMETRY_CELL_TRANSFORM, then removing duplicates.
        Also assigns the duplicateOf variable for all atoms in self.atlist.
        By default, the transformation matrix will be taken from rp, but a
        different matrix can also be passed."""
        ssl = copy.deepcopy(self)
        ssl.clear_symmetry_and_ucell_history()
        ssl.update_cartesian_from_fractional()
        # reduce dimensions in xy
        transform3 = np.identity(3, dtype=float)
        if transform is not None:
            transform3[:2, :2] = transform
        else:
            transform3[:2, :2] = rp.SYMMETRY_CELL_TRANSFORM
        ssl.ucell = np.dot(ssl.ucell, np.linalg.inv(np.transpose(transform3)))
        ssl.collapse_cartesian_coordinates(update_origin=True)
        ssl.ucell_mod = []
        # if self.ucell_mod is not empty, don't drag that into the new slab.
        # remove duplicates
        ssl.create_sublayers(rp.SYMMETRY_EPS_Z)
        newatlist = []
        for subl in ssl.sublayers:
            i = 0
            while i < len(subl.atlist):
                j = i+1
                baseat = [a for a in self
                          if a.oriN == subl.atlist[i].oriN][0]
                while j < len(subl.atlist):
                    if subl.atlist[i].isSameXY(subl.atlist[j].cartpos[:2],
                                               eps=rp.SYMMETRY_EPS):
                        for a in [a for a in self
                                  if a.oriN == subl.atlist[j].oriN]:
                            a.duplicateOf = baseat
                        subl.atlist.pop(j)
                    else:
                        j += 1
                i += 1
            newatlist.extend(subl.atlist)
        ssl.atlist = newatlist
        ssl.update_element_count()   # update number of atoms per element again
        # update the layers. Don't use Slab.createLayers here to keep it
        #   consistent with the slab layers
        for i, layer in enumerate(ssl.layers):
            layer.slab = ssl
            layer.getLayerPos()
            layer.num = i
            layer.atlist = [at for at in layer.atlist if at in ssl]
        return ssl

    def project_c_to_z(self):                                                   # TODO: surface only?
        """Make the c vector of the unit cell perpendicular to the surface.

        All atom coordinates are updated to fit the new basis.

        Returns
        -------
        None.
        """
        c_vec_xy = self.ucell.T[2, :2]
        if any(c_vec_xy):  # Non-zero components
            self.update_cartesian_from_fractional()
            c_vec_xy[:] = 0
            self.collapse_cartesian_coordinates()  # Also updates fractional

    # def update_atom_numbers(self):
    def resetAtomOriN(self):
        """Gets new 'original' numbers for atoms in the slab. If a bulkslab
        is defined, also updates the numbers there to keep the two consistent.
        """
        self.sort_original()
        self.sort_by_element()
        bulkAtsRenumbered = []
        for (i, at) in enumerate(self):
            if self.bulkslab is not None:
                for bat in [a for a in self.bulkslab
                            if a.oriN == at.oriN
                            and a not in bulkAtsRenumbered]:
                    bat.oriN = i+1
                    bulkAtsRenumbered.append(bat)
            at.oriN = i+1

    def revertUnitCell(self, restoreTo=None):
        """If the unit cell in a and b was transformed earlier, restore the
        original form and coordinates. If a 'restoreTo' argument is passed,
        restore only back to the point defined by the argument."""
        if restoreTo is None:
            restoreTo = []
        if len(self.ucell_mod) > 0:
            self.update_cartesian_from_fractional()
            oplist = self.ucell_mod[len(restoreTo):]
            for op in list(reversed(oplist)):
                if op[0] == 'add':
                    for at in self:
                        at.cartpos[0:2] -= op[1]
                    self.collapse_cartesian_coordinates()
                elif op[0] == 'lmul':
                    self.ucell = np.dot(np.linalg.inv(op[1]), self.ucell)
                    self.collapse_cartesian_coordinates()
                elif op[0] == 'rmul':
                    self.ucell = np.dot(self.ucell, np.linalg.inv(op[1]))
                    self.collapse_cartesian_coordinates()
            self.ucell_mod = self.ucell_mod[:len(restoreTo)]

    def sort_by_element(self):
        """Sort `slab.atlist` by element, preserving the element order."""
        _map = {e: i for i, e in enumerate(self.elements)}
        try:
            self.atlist.sort(key=lambda atom: _map[atom.el])
        except KeyError as missing_el:
            _err = ('Unexpected point encountered in '
                    f'{type(self).__name__}.sort_by_element: Could '
                    f'not find element {missing_el} in element list')
            _LOGGER.error(_err)
            raise SlabError(
                'Perhaps you added some atoms and did not '
                f'update_element_count()? {self.elements=}'
                ) from missing_el

    def sort_by_z(self, bottom_to_top=True):
        """Sort `slab.atlist` by z coordinate."""
        self.atlist.sort(key=lambda at: at.pos[2], reverse=not bottom_to_top)

    def sort_original(self):
        """Sort `slab.atlist` by original atom order from POSCAR."""
        self.atlist.sort(key=attrgetter('oriN'))

    def update_cartesian_from_fractional(self, update_origin=False):            # TODO: should we do anything to the .bulkslab too?
        """Assign absolute Cartesian coordinates to all atoms.

        The frame of reference has (x, y) as defined by the a and b
        unit-cell vectors, z == 0 for the topmost atom and positive
        going down through the slab (from the surface into the solid).

        This method can be used, e.g., to (re)calculate all Cartesian
        coordinates when distortions have been applied to the unit
        cell (but atoms should retain their position relative to the
        modified cell).

        Parameters
        ----------
        update_origin : bool, optional
            Whether the z component of the Cartesian origin should
            be updated. If True, `topat_ori_z` is modified to be the
            current z coordinate of the topmost atom. This value is
            ignored if `topat_ori_z` was never initialized. Set this
            to True if the topmost atom is likely to have moved in z.
            Default is False.

        Returns
        -------
        None.
        """
        if update_origin or self.topat_ori_z is None:
            topat = max(self, key=lambda atom: atom.pos[2])
            self.topat_ori_z = np.dot(self.ucell, topat.pos)[2]
        for atom in self:
            atom.cartpos = np.dot(self.ucell, atom.pos)
            atom.cartpos[2] = self.topat_ori_z - atom.cartpos[2]

    def _update_chem_elements(self, rpars):                                     # TODO: @fkraushofer why aren't we also taking into account ELEMENT_RENAME here?
        """Update elements based on the ELEMENT_MIX parameter.

        Issue warnings in case of a naming conflict.

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS to be used. Only the ELEMENT_MIX
            attribute is used.

        Returns
        -------
        None.
        """
        if self._last_element_mix == rpars.ELEMENT_MIX:
            return     # don't update if up to date

        # Check for overlapping element names
        overlapping = {e for mix in rpars.ELEMENT_MIX.values() for e in mix}
        overlapping &= set(self.elements)
        for element in overlapping:
            _LOGGER.warning(
                f'Element name {element} given in ELEMENT_MIX is '
                'also an element name in POSCAR. It is recommended '
                'you rename the element in the POSCAR file.'
                )
        self.chemelem = set()
        for element in self.elements:
            new_els = rpars.ELEMENT_MIX.get(element, (element,))
            self.chemelem.update(e.capitalize() for e in new_els)
        self._last_element_mix = copy.deepcopy(rpars.ELEMENT_MIX)

    def update_element_count(self):
        """Update the number of atoms per element."""
        n_per_elem = Counter(at.el for at in self)

        # Keep original element order, if there were any elements
        elements = self.elements or n_per_elem
        self.n_per_elem = {el: n_per_elem[el]
                           for el in elements
                           if el in n_per_elem}

    def update_fractional_from_cartesian(self):                                 # TODO: should we do anything to the .bulkslab too?
        """Calculate atoms' fractional coordinates from Cartesian ones.

        This method is commonly used when `ucell` is modified, and
        the atoms should retain their absolute Cartesian positions
        while recomputing their fractional coordinates relative to
        the new unit cell.

        Returns
        -------
        None.
        """
        uci = np.linalg.inv(self.ucell)
        for atom in self:
            # Flip over the z Cartesian coordinate, as
            # we store it as "positive going down"
            cartpos = atom.cartpos.copy()
            cartpos[2] = self.topat_ori_z - cartpos[2]                          # TODO: edit when flipping .cartpos[2]
            atom.pos = np.dot(uci, cartpos)

    def update_layer_coordinates(self):
        """Update the Cartesian position of all `layers`."""
        for layer in self.layers:
            layer.getLayerPos()

    # ----------------------  TRANSFORMATIONS  ------------------------

    def apply_matrix_transformation(self, trafo_matrix):
        """Apply an orthogonal transformation to the unit cell and all atoms.

        The transformation is given as an orthogonal transformation
        matrix (O) which is applied to BOTH the unit cell and all
        Cartesian atomic coordinates. The unit cell (U, unit vectors
        as columns) is transformed to U' = O @ U. Atomic coordinates
        (v, as column vectors) are transformed to v' = O @ v. This
        transformation is essentially equivalent to a change of basis.

        This method differs from `rotate_atoms`, `mirror_atoms` and
        `rotateUnitCell` in that the former two only cause a
        rotation of the atoms, but not of the unit cell, whereas
        the latter rotates the unit cell but not the atoms. Here both
        unit cell and atoms are transformed.

        If the transformation is an out-of-plane rotation/mirror (i.e.,
        it changes the z components of unit vectors), layers, bulkslab,
        and sublayers are discarded and will need to be recalculated.
        Otherwise, the same coordinate transform is also applied to
        the `bulkslab`, if present.

        Parameters
        ----------
        trafo_matrix : Sequence
            `trafo_matrix` must be an orthogonal 3-by-3 matrix.
            Contains the transformation matrix (O) describing
            the applied transformation.

        Raises
        ------
        ValueError
            If `trafo_matrix` is not 3-by-3 or not orthogonal.

        Examples
        --------
        Apply a rotation by 90 deg around the z axis to the unit cell
        (in positive direction, i.e. clockwise when looking along z)
        >>> theta = np.pi/2
        >>> rot_mat = [[np.cos(theta), -np.sin(theta), 0],
                       [np.sin(theta),  np.cos(theta), 0],
                       [0, 0, 1]]
        >>> slab.apply_matrix_transformation(rot_mat)
        """
        trafo_matrix = np.asarray(trafo_matrix)
        if trafo_matrix.shape != (3, 3):
            raise ValueError('apply_matrix_transformation: '
                             'not a 3-by-3 matrix')
        if not np.allclose(np.linalg.inv(trafo_matrix), trafo_matrix.T):
            raise ValueError('apply_matrix_transformation: matrix is not '
                             'orthogonal. Consider using apply_scaling.')

        # Determine whether trafo_matrix will change
        # the z component of the unit vectors
        changes_z = not np.allclose(trafo_matrix[2], (0, 0, 1))

        self.ucell = trafo_matrix.dot(self.ucell)
        self.ucell[abs(self.ucell) < 1e-5] = 0.
        self.update_cartesian_from_fractional(update_origin=changes_z)

        # Update also 'layers', sublayers and bulkslab: if the
        # transformation touched 'z' we invalidate everything
        if changes_z:
            self.layers.clear()
            self.sublayers.clear()

        if self.is_bulk:
            return

        if changes_z:
            self.bulkslab = None
        elif self.bulkslab:
            self.bulkslab.apply_matrix_transformation(trafo_matrix)

    def apply_scaling(self, *scaling):
        """Rescale the unit-cell vectors.

        This can be used to stretch/compress along unit cell vectors
        in order to change lattice constants in some direction or to
        apply an isotropic scaling in all directions. To apply other
        (orthogonal) transformations (e.g., rotation, flipping), use
        `apply_matrix_transformation`.

        The same scaling is also applied to `bulkslab`, if this slab
        has one.

        Parameters
        ----------
        *scaling : Sequence
            If only one number, an isotropic scaling is applied to
            the unit cell and atom positions. If a sequence with
            three entries, the scaling will be applied along the
            unit-cell vectors in the given order.

        Returns
        -------
        scaling_matrix : numpy.ndarray
            The matrix used for scaling the unit vectors.

        Raises
        ------
        TypeError
            If `scaling` has neither 1 nor three elements, or
            any of the elements is not a number.
        ValueError
            If `scaling` would make the unit cell singular
            (i.e., reduce the length of a unit vector to zero)

        Examples
        ----------
        Stretch the unit cell by a factor of 2 along a, b and c.
        This doubles the lattice constant and increases the volume
        8-fold:
        >>> slab.apply_scaling(2)

        Compresses the unit cell by a factor of 3 along c:
        >>> slab.apply_scaling(1, 1, 1/3)
        """
        if len(scaling) not in (1, 3):
            raise TypeError(f'{type(self).__name__}.apply_scaling: '
                            'invalid number of arguments. Expected '
                            f'one or three, got {len(scaling)}.')
        if not all(isinstance(s, Real) for s in scaling):
            raise TypeError(f'{type(self).__name__}.apply_scaling: '
                            f'invalid scaling factor. Expected one '
                            'or three numbers.')
        if len(scaling) == 1:
            scaling *= 3
        if any(abs(s) < 1e-5 for s in scaling):
            raise ValueError(f'{type(self).__name__}.apply_scaling: cannot '
                             'reduce unit vector(s) to zero length')

        # Apply to unit cell (basis). Notice the inverted order,
        # because the unit cell is stored with unit vectors as
        # columns (i.e., a = ucell[:, 0])
        scaling_matrix = np.diag(scaling)
        self.ucell = self.ucell.dot(scaling_matrix)
        self.update_cartesian_from_fractional(update_origin=scaling[2] != 1)

        try:
            self.bulkslab.apply_scaling(*scaling)
        except AttributeError:
            pass
        return scaling_matrix

    def mirror_atoms(self, symplane, glide=False):
        """Apply a mirror or glide transform across `symplane`.

        Notice that the unit cell stays unchanged. Coordinates
        are collapsed to the base cell after transformation.

        Parameters
        ----------
        symplane : SymPlane
            The mirror/glide plane. If `symplane.is_glide` a
            glide operation is performed rather than a mirror.
        glide : bool, optional
            Whether a translation should also be applied, to
            realize a glide-symmetry transformation. This is
            ignored if `symplane.is_glide`. Default is False.

        Returns
        -------
        None.
        """
        mirror = symplane.point_operation(n_dim=3)
        glide_vec = np.zeros(2)
        if symplane.is_glide or glide:
            glide_vec = symplane.glide_vector
        self._transform_atoms_2d(mirror, center=symplane.pos, shift=glide_vec)

    def rotate_atoms(self, order, axis=(0., 0.)):
        """Apply an `order`-fold 2D rotation around axis.

        Notice that the unit cell stays unchanged. Coordinates
        are collapsed to the base cell after transformation.

        Parameters
        ----------
        order : int
            Order of rotation, as accepted by `rotation_matrix_order`.
            When looking down at the slab from vacuum, atoms will be
            rotated counter-clockwise by 2pi/`order` radians.
        axis : Sequence, optional
            In-plane Cartesian position of the rotation axis. Shape
            should be (2,). Default is (0, 0) (i.e., the origin).

        Returns
        -------
        None.
        """
        rotation = rotation_matrix_order(order, dim=3)
        self._transform_atoms_2d(rotation, center=axis, shift=0)

    def rotateUnitCell(self, order, append_ucell_mod=True):
        """Rotates the unit cell (around the origin), leaving atom positions
        the same. Note that this rotates in the opposite direction as
        rotate_atoms."""
        self.update_cartesian_from_fractional()
        m = rotation_matrix_order(order)
        m3 = np.identity(3, dtype=float)
        m3[:2, :2] = m
        self.ucell = np.dot(m3, self.ucell)
        if append_ucell_mod:
            self.ucell_mod.append(('lmul', m3))
        self.update_fractional_from_cartesian()

    def _transform_atoms_2d(self, transform, center, shift):
        """Apply matrix `transform` at `center` to all atoms, then `shift`.

        Fractional coordinates are collapsed (in 2D) at the end.

        Parameters
        ----------
        transform : Sequence
            Shape (3, 3). 2D transformation matrix. Will be applied
            to the left to the Cartesian coordinates of each atom,
            taken as a column vector.
        center : Sequence
            Shape (2,). Atoms are translated to `center` before applying
            `transform`, then translated back.
        shift : float or Sequence
            Extra translation to be applied to all atoms at the end
            of the centred transform `transform`. If a sequence, should
            have shape (2,).

        Returns
        -------
        None.
        """
        # All the transformations are done on the in-plane
        # coordinates only. Make sure we have the out-of-plane
        # right as well.
        if self.topat_ori_z is None:
            self.update_cartesian_from_fractional()

        # Let's work with fractional coordinates, as we can
        # collapse them. Also, we will work with coordinates
        # as row vectors. Cartesians transform as:
        #   c' = (c - center) @ T + center + shift
        #      = c @ T + center @ (I - T) + shift
        # and
        #   c = p @ u_cell,    or   p = c @ u_inv
        # Thus
        #   p' = c' @ u_inv = p @ (u_cell @ T @ u_inv) + frac_offset
        # where
        #   frac_offset = (center @ (I - T) + shift) @ u_inv
        #
        # Important note: All the operations above are in principle
        # meant to be done using the full 3D unit cells, fractional
        # and Cartesian positions. However we can stick to using only
        # the in-plane components for the frac_offset, as both center
        # and shift are only in 2D.
        ab_inv = np.linalg.inv(self.ab_cell.T)
        transform = transform.T  # Coordinates as row vectors
        frac_offset = np.dot(center, np.identity(2) - transform[:2, :2])
        frac_offset += shift
        frac_offset = frac_offset.dot(ab_inv)

        # Now the part that needs to be done in 3D. This is critical
        # for systems that have c axis not orthogonal to the surface,
        # as both the frac_transform and the fractional-to-Cartesian
        # conversions mix the third component into the in-plane ones.
        ucell = self.ucell.T
        frac_transform = ucell.dot(transform).dot(np.linalg.inv(ucell))
        for atom in self:
            atom.pos[:2] = atom.pos.dot(frac_transform)[:2] + frac_offset
            collapse_fractional(atom.pos[:2], in_place=True)
            atom.cartpos[:2] = atom.pos.dot(self.ucell.T)[:2]

    # ----------------- SYMMETRY UPON TRANSFORMATION ------------------

    def _is_2d_transform_symmetric(self, matrix, center, translation, eps):
        """Return whether `self` is identical under a 2D symmetry transform.

        Parameters
        ----------
        matrix : numpy.ndarray
            Shape (2, 2). The point operation `matrix` to be applied
            to all Cartesian coordinates (as row vectors, on the right)
        center : Sequence
            Shape (2,). The point of application of the operation
            in `matrix`.
        translation : Sequence
            Shape (2,). Translation vector to be applied.
        eps : float
            Tolerance (Cartesian) for position equivalence

        Returns
        -------
        is_symmetric : bool
            Whether the slab is self-similar under the operation.

        Raises
        ------
        NeedsSublayersError
            If called before sublayers were created
        """
        # Use the version of the unit cell with unit vectors as rows
        ab_cell = self.ab_cell.T
        ab_inv = np.linalg.inv(ab_cell)
        releps = eps / np.linalg.norm(ab_cell, axis=1)

        # Run the comparison sublayer-wise
        if not self.sublayers:
            raise NeedsSublayersError(
                '2d-transform invariance check requires sublayers. '
                'Call create_sublayers(epsz), then try again.'
                )

        for layer in self.sublayers:
            # Collapse fractional coordinates before transforming
            cart_coords = np.array([atom.cartpos[:2] for atom in layer])
            cart_coords, frac_coords = collapse(cart_coords, ab_cell, ab_inv)

            # Create a transformed copy: shift, transform,
            # shift back, and apply rigid translation
            transf_coords = cart_coords - center
            transf_coords = transf_coords.dot(matrix) + center + translation

            # Collapse again
            transf_coords, _ = collapse(transf_coords, ab_cell, ab_inv)

            # Add extra atoms close to edges/corners for comparing
            cart_coords, _ = add_edges_and_corners(cart_coords,
                                                   frac_coords,
                                                   releps, ab_cell)
            # Finally compare interatomic distances
            distances = euclid_distance(transf_coords, cart_coords)
            if any(distances.min(axis=1) > eps):
                return False
        return True

    def is_mirror_symmetric(self, symplane, eps, glide=False):
        """Return if this slab is unchanged when applying a 2D mirror/glide.

        Parameters
        ----------
        symplane : SymPlane
            The plane across which a mirror/glide should be
            performed. If `symplane.is_glide`, a glide operation
            is applied.
        eps : float
            Cartesian tolerance for position equivalence.
        glide : bool, optional
            Whether a glide translation is applied. This value
            is ignored if `symplane.is_glide`.

        Returns
        -------
        symmetric : bool
            Whether slab is unchanged when 'mirrored' across `symplane`
        """
        matrix = symplane.point_operation(n_dim=2)
        glide_vec = np.zeros(2)
        if symplane.is_glide or glide:
            glide_vec = symplane.glide_vector
        return self._is_2d_transform_symmetric(matrix, symplane.pos,
                                               glide_vec, eps)

    def is_rotation_symmetric(self, axis, order, eps):
        """Return if this slab is unchanged when applying a 2D rotation.

        Parameters
        ----------
        axis : Sequence
            The 2D Cartesian coordinates of the axis around
            which the rotation should be tested.
        order : int
            The order of the rotation to be tested.
        eps : float
            Cartesian tolerance for position equivalence.

        Returns
        -------
        symmetric : bool
            Whether slab is unchanged when rotated around `axis`.
        """
        matrix = rotation_matrix_order(order, dim=2)
        return self._is_2d_transform_symmetric(matrix, axis, 0, eps)

    def isTranslationSymmetric(self, tv, eps, z_periodic=True, z_range=None):
        """
        Evaluates whether the slab is equivalent to itself when translated
        along the given cartesian translation vector tv.

        Parameters
        ----------
        tv : numpy array
            2- or 3-dimensional translation vectors are accepted.
        eps : float
            Error tolerance for positions (cartesian)
        z_periodic : bool, optional
            True for checking periodicity of a bulk slab, in which the c vector
            is a true unit cell vector. False otherwise.
        z_range : tuple of floats, optional
            Limit check to only atoms within a given range of cartesian
            coordinates. The default is None.

        Returns
        -------
        bool
            True if translation symmetric, else False.

        """
        if len(tv) == 2:  # two-dimensional displacement. append zero for z
            tv = np.append(tv, 0.)
        uc = np.copy(self.ucell)
        uc[:, 2] *= -1   # mirror c vector down
        uct = np.transpose(uc)
        releps = [eps / np.linalg.norm(uct[j]) for j in range(0, 3)]
        shiftv = tv.reshape(3, 1)
        # unlike in-plane operations, this one cannot be done sublayer-internal
        coordlist = [at.cartpos for at in self]
        shiftm = np.tile(shiftv, len(coordlist))
        oricm = np.array(coordlist)  # original cartesian coordinate matrix
        oricm[:, 2] *= -1
        shiftm[2] *= -1
        oripm = np.dot(np.linalg.inv(uc), oricm.transpose()) % 1.0
        # collapse (relative) coordinates to base unit cell
        oricm = np.dot(uc, oripm).transpose()
        # original cartesian coordinates collapsed to base unit cell
        tmpcoords = np.copy(oricm).transpose()
        # copy of coordinate matrix to be manipulated
        # determine which z to check
        if z_range is None:
            min_z = np.min(tmpcoords[2]) - eps
            max_z = np.max(tmpcoords[2]) + eps
        else:
            z_range = tuple(-v for v in z_range)
            min_z, max_z = min(z_range) - eps, max(z_range) + eps
        tmpcoords += shiftm
        if not z_periodic:
            # discard atoms that moved out of range in z
            tmpcoords = tmpcoords[:, tmpcoords[2] >= min_z]
            tmpcoords = tmpcoords[:, tmpcoords[2] <= max_z]
        tmpcoords = np.dot(uc, (np.dot(np.linalg.inv(uc), tmpcoords) % 1.0))
        # collapse coordinates to base unit cell
        # for every point in matrix, check whether is equal:
        for (i, p) in enumerate(oripm.transpose()):
            # get extended comparison list for edges/corners:
            addlist = []
            for j in range(0, 3):
                if abs(p[j]) < releps[j]:
                    addlist.append(oricm[i]+uct[j])
                if abs(p[j]-1) < releps[j]:
                    addlist.append(oricm[i]-uct[j])
            if len(addlist) == 2:
                # 2D coner - add the diagonally opposed point
                addlist.append(addlist[0]+addlist[1]-oricm[i])
            elif len(addlist) == 3:
                # 3D corner - add all diagonally opposed points
                addlist.extend([(p1 + p2 - oricm[i]) for (p1, p2) in
                                itertools.combinations(addlist, 2)])
                addlist.append(addlist[0] + addlist[1] + addlist[2]
                               - 2*oricm[i])
            for v in addlist:
                oricm = np.concatenate((oricm, v.reshape(1, 3)))
        distances = euclid_distance(tmpcoords.transpose(), oricm)
        # print(oricm)
        if any(min(sublist) > eps for sublist in distances):
            return False
        return True
