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

from abc import abstractmethod
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
from viperleed.tleedmlib.base import add_edges_and_corners, collapse
from viperleed.tleedmlib.base import collapse_fractional, pairwise
from viperleed.tleedmlib.base import rotation_matrix_order
from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.atom_containers import AtomContainer, AtomList
from viperleed.tleedmlib.classes.layer import Layer, SubLayer
from viperleed.tleedmlib.classes.sitetype import Sitetype

from .slab_errors import AlreadyMinimalError, InvalidUnitCellError
from .slab_errors import MissingLayersError, MissingSublayersError, SlabError
from .slab_utils import _left_handed, _z_distance


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
class BaseSlab(AtomContainer):
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
        self.atlist = AtomList()                                                # base
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
        if isinstance(item, SubLayer):
            return item in self.sublayers
        if isinstance(item, Layer):
            return item in self.layers
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
        return [lay for lay in self.layers if lay.is_bulk]

    @property
    def elements(self):
        """Return a tuple of elements in this slab, as originally read."""
        return tuple(self.n_per_elem.keys())

    @property
    def fewest_atoms_sublayer(self):
        """Return the sublayer with fewest atoms."""
        if not self.sublayers:
            raise MissingSublayersError(f'{type(self).__name__} has '
                                        'no sublayers defined')
        return min(self.sublayers, key=attrgetter('n_atoms'))

    @property
    @abstractmethod
    def is_bulk(self):
        """Return whether this is a bulk slab."""
        return False

    @property
    def n_atoms(self):
        """Return the number of atoms in this slab."""
        return self.atlist.n_atoms

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
        MissingLayersError
            If no layers are available
        """
        if not self.layers:
            raise MissingLayersError

        if self.n_layers == 1:
            return 0.                                                           # TODO: I don't think it's right that it is zero if there's only one layer. Think about it.
        # self.update_cartesian_from_fractional()                               # TODO: I don't think this is needed. It does not update anything for layers; only makes sense if we also .update_layer_coordinates.

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

    def _add_one_bulk_cell(self, bulk_layers, bulkc_par,
                           bulkc_par_atoms, bulkc_perp_atoms,
                           new_atoms_start_index=None):
        """Add one bulk unit cell and return the new atoms.

        This method is intended to be used only internally in
        Slab objects. Use with_extra_bulk_units (for non-bulk
        slabs) or with_double_thickness (for bulk ones).

        Parameters
        ----------
        bulk_layers : list of Layer
            The layers of the bottommost bulk cell to be repeated.
        bulkc_par : numpy.ndarray
            Component of the bulk repeat vector parallel to the
            c axis of the slab. This is used to expand the unit
            cell. Shape (3,).
        bulkc_par_atoms : numpy.ndarray
            Component of the bulk repeat vector parallel to the c axis,
            to be used for the Cartesian positions of the atoms. The
            difference with `bulkc_par` is only due to the fact that we
            store Cartesian coordinates of atoms with z going INTO THE
            SOLID, while the unit cell reference has z going TOWARD
            THE SURFACE. Shape (3,).
        bulkc_perp_atoms : numpy.ndarray
            Component of the bulk repeat vector perpendicular to the
            c axis, to be used for the Cartesian positions of the
            atoms. Shape (3,).
        new_atoms_start_index : int or None, optional
            The initial index to be used as num for the atoms that
            are added as a consequence of the duplication.

        Returns
        -------
        new_atoms : list
            The atoms added to this slab.
        new_layers : list
            The Layer objects added to this slab.
        """
        if new_atoms_start_index is None:
            new_atoms_start_index = self.n_atoms + 1
        # The mess with needing three different components of the
        # bulk repeat vector comes from the fact that we want to
        # add atoms (and layers) below, but not shift the current
        # ones relative to the current unit cell. This is useful
        # because the 'old atoms' appear then at the same in-plane
        # position in, e.g., VESTA.
        self.ucell.T[2] += bulkc_par  # Expand unit cell
        new_atoms = []
        new_layers = []
        for layer in bulk_layers.copy():  # .copy avoids infinite loop
            # Add a new bulk layer and duplicate of all its atoms
            new_layer = Layer(self, self.n_layers, True)
            new_layers.append(new_layer)
            self.layers.append(new_layer)
            for atom in layer:
                new_atom = atom.duplicate(add_to_atlists=False,
                                          num=new_atoms_start_index)
                new_atoms_start_index += 1
                new_layer.atlist.append(new_atom)
                new_atom.layer = new_layer
                self.atlist.append(new_atom)
                self.n_per_elem[new_atom.el] += 1
                new_atoms.append(new_atom)
                # Shift new atoms only perpendicular to unit-cell c
                new_atom.cartpos -= bulkc_perp_atoms

        # Shift old atoms up along unit-cell c
        for old_atom in self:
            if old_atom not in new_atoms:
                old_atom.cartpos += bulkc_par_atoms

        # Finally, recalculate the positions of all layers
        self.update_layer_coordinates()
        return new_atoms, new_layers

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
        tmplist = AtomList(*self.atlist)
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
            if layer.is_bulk:
                self.layers[layer.num+1].is_bulk = True
            self.layers.remove(layer)
            del layer
        self.layers.reverse()
        for i, layer in enumerate(self.layers):
            layer.update_position()
            layer.num = i
        self.atlist = tmplist
        return ct

    def _get_sublayers_for_el(self, element, eps):
        """Yield SubLayer objects for atoms of `element` within `eps`."""
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
            new_layer = SubLayer(self, 0)
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
            same_z.sort(key=attrgetter('element'))                              # TODO: is this even necessary? .sort should be stable, so element order should be preserved
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

    def full_update(self, rpars):
        """Update atoms and layers from `rpars`.

        This method is typically used after information of a
        slab is read from file (e.g., via `poscar.read`).

        The method ensures that all atoms are within the (0, 0) unit
        cell, then, if needed, uses the information from `rpars` to
        create layers calculate Cartesian coordinates (absolute and
        per layer), and to update elements and sites.
        Notice that atom displacements are NOT cleared, unless
        rpars.fileLoaded['VIBROCC'] is True-thy.

        Parameters
        ----------
        rpars : Rparams
            Information read from a PARAMETERS file.

        Returns
        -------
        None.
        """
        self.collapse_fractional_coordinates()
        self.update_cartesian_from_fractional()
        if not self.layers:
            self.createLayers(rpars)
        else:
            # It is only needed if there were layers already and
            # either (i) fractional coordinates had pos[2] out of
            # range, or (ii) topat_ori_z had not been calculated yet
            self.update_layer_coordinates()
        self._update_chem_elements(rpars)
        if not self.sitelist:
            self.initSites(rpars)
        if rpars.fileLoaded['VIBROCC']:
            for atom in self:
                atom.initDisp()

    @abstractmethod
    def get_bulk_repeat(self, rpars, only_z_distance=False):
        """Return the bulk repeat vector (with positive z).

        Notice that this method does **not attempt to identify an
        unknown bulk-repeat vector**.

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS to be interpreted.
        only_z_distance : bool, optional
            Whether a distance in the direction perpendicular to the
            surface (i.e., not necessarily along the c axis) should
            be returned rather than a full vector. This is ignored
            if `rpars.BULK_REPEAT` is a vector. Default is False.

        Returns
        -------
        bulk_repeat_vector : numpy.ndarray or float
            Bulk repeat vector pointing from the bulk to the surface,
            or its component along z. If `rpars.BULK_REPEAT` is a
            vector, a copy is returned. Otherwise, the vector is taken
            to be parallel to the c axis of this slab. Its length is
            calculated from a z distance (i.e., perpendicular to the
            surface). The z distance is taken as either the value of
            `rpars.BULK_REPEAT` (if it is not-None) or the z distance
            between the bottommost points of the lowest bulk and lowest
            non-bulk layers.
        """

    def get_minimal_ab_cell(self, eps, epsz=None, warn_convention=False):       # TODO: write a test case for the reduction of POSCAR Sb on Si(111)  # too-many-locals
        """Check if there is a 2D unit cell smaller than the current one.

        Parameters
        ----------
        eps : float
            Cartesian tolerance for in-plane comparisons (Angstrom)
        epsz : float or None, optional
            Cartesian tolerance for comparisons in the direction
            perpendicular to the surface. If not given on None, it
            is take equal to `eps`. Default is None.
        warn_convention : bool, optional
            If True, warnings are added to the current logger in case
            making the reduced unit cell stick to the conventions would
            result in a sub-optimal `SUPERLATTICE` matrix. Default is
            False.

        Returns
        -------
        mincell : np.ndarray
            The reduced 2D unit cell. Notice that `mincell` is such that
            `a`, `b` = `mincell`, i.e., it is transposed with respect
            to `slab.ab_cell`.

        Raises
        ------
        AlreadyMinimalError
            If the current 2D unit cell cannot be reduced.
        """
        if epsz is None:
            epsz = eps

        # Create a test slab: C projected to Z
        ts = copy.deepcopy(self)
        ts.project_c_to_z()
        ts.sort_by_z()
        ts.create_sublayers(epsz)

        # Use the lowest-occupancy sublayer (the one
        # with fewest atoms of the same chemical element)
        lowocclayer = ts.fewest_atoms_sublayer
        n_atoms = lowocclayer.n_atoms
        if n_atoms < 2:
            # Cannot be smaller if there's only 1 atom
            raise AlreadyMinimalError(
                f'{type(self).__name__}.get_minimal_ab_cell: fewest-atom '
                f'sublayer ({lowocclayer.num}) contains only one atom '
                f'({lowocclayer.atlist[0]}).'
                )

        # Create a list of candidate unit vectors as those connecting
        # atom pairs. Notice that it is enough to take any arbitrary
        # atom as a 'reference' as, if the unit cell can be reduced,
        # all atoms must have a 'copy'. We take the first one.
        plist = [at.cartpos[:2] for at in lowocclayer]
        candidate_unit_vectors = (p - plist[0] for p in plist[1:])              # TODO: since is_translation_symmetric is somewhat expensive, would it make sense to pre-filter the vectors keeping only the shortest ones among those parallel to one another?
        candidate_unit_vectors = (vec for vec in candidate_unit_vectors
                                  if ts.is_translation_symmetric(vec, eps))
        # Try to reduce the cell
        smaller, mincell = self._minimize_ab_area(ts.ab_cell.T, n_atoms,
                                                  candidate_unit_vectors, eps)
        if not smaller:
            raise AlreadyMinimalError(
                f'{type(self).__name__}.get_minimal_ab_cell: none '
                f'of the translation vectors tested gave a smaller area'
                )
        return self._reduce_mincell(mincell, eps, warn_convention)

    @staticmethod
    def _minimize_ab_area(ab_cell, n_atoms, candidate_unit_vectors, eps):
        """Reduce an ab_cell with n_atoms to the smallest area possible.

        Parameters
        ----------
        ab_cell : numpy.ndarray
            The cell to be reduced
        n_atoms : int
            The number of atoms in ab_cell. Used to determine the
            maximum downscaling that can be expected.
        candidate_unit_vectors : Iterable
            2D Vectors to be tested for reduction of ab_cell.
        eps : float
            Linear Cartesian tolerance to discern whether areas
            differ from one another. Its square is used for areas.

        Returns
        -------
        smaller : bool
            Whether ab_cell could be reduced.
        mincell : numpy.ndarray
            The reduced unit cell.
        """
        # Try to reduce the cell by using a pair of vectors from
        # [a, b, *candidate_unit_vectors]. Keep in mind that with
        # n_atoms, we cannot reduce the area by more than a factor
        # 1 / n_atoms (which would give one atom per reduced cell)
        ab_cell_area = abs(np.linalg.det(ab_cell))
        smallest_area = ab_cell_area / n_atoms
        smaller = False
        for vec in candidate_unit_vectors:
            # Try first replacing the current second unit vector
            tcell = np.array([ab_cell[0], vec])
            tcell_area = abs(np.linalg.det(tcell))
            if smallest_area < tcell_area + eps**2 < ab_cell_area:
                ab_cell, ab_cell_area = tcell, tcell_area
                smaller = True
                continue
            # Try replacing the current first unit vector instead
            tcell = np.array([ab_cell[1], vec])
            tcell_area = abs(np.linalg.det(tcell))
            if smallest_area < tcell_area + eps**2 < ab_cell_area:
                ab_cell, ab_cell_area = tcell, tcell_area
                smaller = True
        return smaller, ab_cell

    @staticmethod
    def _reduce_mincell(mincell, eps, warn_convention):
        """Ensure `mincell` is Minkowski-reduced; enforce conventions."""
        # Use Minkowski reduction to make mincell high symmetry
        mincell, *_ = leedbase.reduceUnitCell(mincell)

        # Cosmetic corrections
        if all(abs(np.diag(mincell)) < eps):
            # Swap a and b when matrix is off-diagonal
            mincell[[0, 1]] = mincell[[1, 0]]

        _is_diagonal = abs(mincell[1, 0]) < eps and abs(mincell[0, 1]) < eps
        if _is_diagonal:
            # If matrix is diagonal, make elements positive
            mincell = abs(mincell)
        # By convention, make the shorter vector the first one
        if np.linalg.norm(mincell[0]) > np.linalg.norm(mincell[1]) + eps:
            # If matrix is diagonal, DO NOT make it off-diagonal
            if not _is_diagonal:
                mincell = np.dot([[0, 1], [-1, 0]], mincell)
            elif _is_diagonal and warn_convention:
                _LOGGER.warning('The unit cell orientation does not follow '
                                'standard convention: to keep SUPERLATTICE '
                                'matrix diagonal, the first bulk vector must '
                                'be larger than the second. Consider swapping '
                                'the unit cell vectors.')
        # Finally, make sure it's right-handed
        if _left_handed(mincell):
            mincell = np.dot([[1, 0], [0, -1]], mincell)
        return mincell

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
        atlist.sort(key=attrgetter('num'))
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

    def is_equivalent(self, other, eps=0.001):
        """Return whether this slab is equivalent to another.

        Parameters
        ----------
        other : Slab
            The other slab to compare to.
        eps : float, optional
            Tolerance on Cartesian atomic coordinates.

        Returns
        -------
        equivalent : bool
            True if the Cartesian coordinates of all atoms match
            (with at least one other atom, if there are duplicates),
            False otherwise. Both slabs are copied and collapsed to
            the (0, 0) cell.
        """
        if not isinstance(other, BaseSlab):
            return False

        slabs = copy.deepcopy(self), copy.deepcopy(other)
        for slab in slabs:
            slab.collapse_cartesian_coordinates()
            if not slab.sublayers:
                slab.create_sublayers(eps)
            # Reorder sublayers by Z to then compare by index                   # TODO: is this necessary?
            slab.sublayers.sort(key=attrgetter('cartbotz'))

        if self.n_sublayers != other.n_sublayers:
            return False

        ab_cell = self.ab_cell.T
        releps = eps / np.linalg.norm(ab_cell, axis=1)

        for this_lay, other_lay in zip(*(s.sublayers for s in slabs)):
            if (this_lay.n_atoms != other_lay.n_atoms
                    or abs(this_lay.cartbotz - other_lay.cartbotz) > eps
                    or this_lay.element != other_lay.element):
                return False
            # Prepare Cartesian coordinates for comparisons
            this_coords = np.array([at.cartpos[:2] for at in this_lay])
            other_coords = np.array([at.cartpos[:2] for at in other_lay])

            # Add extra atoms at edges/corners
            frac_coords = [at.pos[:2] for at in this_lay]
            this_coords, _ = add_edges_and_corners(this_coords,
                                                   frac_coords,
                                                   releps, ab_cell)
            # Finally compare interatomic distances
            distances = euclid_distance(this_coords, other_coords)
            if any(distances.min(axis=1) > eps):
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
        Also assigns the duplicate_of variable for all atoms in self.atlist.
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
        newatlist = AtomList()
        for subl in ssl.sublayers:
            i = 0
            while i < subl.n_atoms:
                atom_i = subl.atlist[i]
                baseat = self.atlist.get(atom_i.num)
                j = i + 1
                while j < subl.n_atoms:
                    atom_j = subl.atlist[j]
                    if not atom_i.is_same_xy(atom_j, eps=rp.SYMMETRY_EPS):
                        j += 1
                        continue
                    subl.atlist.pop(j)
                    duplicate_atom = self.atlist.get(atom_j.num, None)
                    if duplicate_atom:
                        duplicate_atom.duplicate_of = baseat
                i += 1
            newatlist.extend(subl)
        ssl.atlist = newatlist
        ssl.update_element_count()   # update number of atoms per element again
        # update the layers. Don't use Slab.createLayers here to keep it
        #   consistent with the slab layers
        for i, layer in enumerate(ssl.layers):
            layer.slab = ssl
            layer.update_position()
            layer.num = i
            layer.atlist = [at for at in layer if at in ssl]
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
        self.atlist.sort(key=attrgetter('num'))

    def update_atom_numbers(self):
        """Assign new incremental numbers to atoms in the slab.

        If a `bulkslab` is defined, also update the numbers there
        to keep the two consistent. This method is especially
        useful when atoms are added and/or removed from the slab.

        Returns
        -------
        None.
        """
        self.sort_original()
        self.sort_by_element()
        try:
            bulk_atoms = self.bulkslab.atlist
        except AttributeError:
            bulk_atoms = AtomList()

        for i, atom in enumerate(self, start=1):
            bulk_atom = bulk_atoms.get(atom.num, None)
            if bulk_atom:
                bulk_atom.num = i
            atom.num = i
        bulk_atoms.update_atoms_map()
        self.atlist.update_atoms_map()

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
            layer.update_position()

    # ----------------------  TRANSFORMATIONS  ------------------------

    def apply_matrix_transformation(self, trafo_matrix):
        """Apply an orthogonal transformation to the unit cell and all atoms.

        The transformation is given as an orthogonal transformation
        matrix (O) which is applied to BOTH the unit cell and all
        Cartesian atomic coordinates. The unit cell (U, unit vectors
        as columns) is transformed to U' = O @ U. Atomic coordinates
        (v, as column vectors) are transformed to v' = O @ v. This
        transformation is essentially equivalent to a change of basis.

        This method differs from `rotate_atoms`, `mirror_atoms`,
        `rotate_unit_cell`, and `transform_unit_cell_2d` in that
        the former two only cause a transformation of the atoms
        but not of the unit cell, whereas the latter two transform
        the unit cell but not the atoms. Here both unit cell and
        atoms are transformed.

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
            glide_vec = symplane.par.dot(self.ab_cell.T) / 2
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

    def rotate_unit_cell(self, order, append_ucell_mod=True):
        """Rotate the unit cell around the origin.

        All atomic Cartesian coordinates are left unchanged. Note that
        this rotates the fractional atomic coordinates in the opposite
        direction as `rotate_atoms`.

        Parameters
        ----------
        order : int
            Order of rotation. A rotation by 2pi/`order` radians
            around the z axis will be performed.
        append_ucell_mod : bool, optional
            Whether the rotation applied should be registered (to
            allow reverting). Default is True.

        Returns
        -------
        None.
        """
        transform = rotation_matrix_order(order, dim=3)
        self.transform_unit_cell_2d(transform, on_row_vectors=False,
                                    append_ucell_mod=append_ucell_mod)

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

    def transform_unit_cell_2d(self, transform, on_row_vectors=True,
                               append_ucell_mod=True):
        """Apply a 2D-transformation matrix to the unit cell.

        All atomic Cartesian coordinates are left unchanged.
        Fractional coordinates are recalculated for the new
        unit cell, but are **not** collapsed.

        Parameters
        ----------
        transform : Sequence
            Shape (2, 2) or (3, 3). The matrix transformation to be
            applied to the unit cell. The transform is always applied
            'from the left' to the unit cell. `on_row_vectors` decides
            how the unit cell vectors are to be transformed.
        on_row_vectors : bool, optional
            Whether `transform` should be applied to a unit cell where
            the unit vectors are rows, or whether it should be applied
            on one with unit vectors as columns. Default is True.
        append_ucell_mod : bool, optional
            Whether the transformation applied should be registered
            (to allow reverting). Default is True.

        Notes
        -----
        Typical applications include transforming the unit cell with
        SUPERLATTICE/SYMMETRY_CELL_TRANSFORM (`on_row_vectors`=True),
        or applying a rotation or mirror (i.e., Hausholder) matrix to
        the unit cell (`on_row_vectors`=False).

        Returns
        -------
        None.
        """
        if self.topat_ori_z is None:
            # It makes sense to run this only if we do not yet
            # have the information about the Cartesian origin,
            # as the Cartesian coordinates stay unaltered.
            self.update_cartesian_from_fractional()
        if np.shape(transform) == (2, 2):
            transform_3d = np.identity(3)
            transform_3d[:2, :2] = transform
        else:
            transform_3d = np.asarray(transform)

        if on_row_vectors:                                                      # TODO: this has to be turned around when switching ucell.T
            # transform_3d is to be multiplied to the left for a unit
            # cell with unit vectors as rows. However, we use columns:
            #    ucell.T = transform_3d @ ucell.T
            #    ucell = (transform_3d @ ucell.T).T
            #          = ucell @ transform_3d.T
            side = 'rmul'
            transform_3d = transform_3d.T
            self.ucell = np.dot(self.ucell, transform_3d)
        else:
            side = 'lmul'
            self.ucell = np.dot(transform_3d, self.ucell)

        if append_ucell_mod:
            self.ucell_mod.append((side, transform_3d))
        self.update_fractional_from_cartesian()

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
        MissingSublayersError
            If called before sublayers were created
        """
        # Use the version of the unit cell with unit vectors as rows
        ab_cell = self.ab_cell.T
        ab_inv = np.linalg.inv(ab_cell)
        releps = eps / np.linalg.norm(ab_cell, axis=1)

        # Run the comparison sublayer-wise
        if not self.sublayers:
            raise MissingSublayersError(
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
            glide_vec = symplane.par.dot(self.ab_cell.T) / 2
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

    def is_translation_symmetric(self, translation, eps,                        # too-many-locals
                                 z_periodic=True, z_range=None):
        """Return if this slab is unchanged when translated.

        Parameters
        ----------
        translation : Sequence
            2- or 3-dimensional Cartesian translation vector
            for which translational symmetry is to be tested.
            The coordinate frame is the same as for atoms.
        eps : float
            Error tolerance for positions (Cartesian).
        z_periodic : bool, optional
            Whether the c unit vector should be considered as a repeat
            vector or not. Should be True for checking periodicity of
            a bulk slab, False otherwise. Default is True.
        z_range : Sequence or None, optional
            When a Sequence, `min(z_range)` and `max(z_range`) are used
            as Cartesian limits. Only atoms with min_z <= z <= max_z
            are compared. `z_range` is used only if `z_periodic` is
            False. Default is None.

        Returns
        -------
        is_symmetric : bool
            Whether this slab is unchanged (within `eps`) upon
            `translation`.

        Raises
        ------
        ValueError
            If `translation` is not a 2- or 3-D vector.
        """
        translation = np.array(translation)
        if translation.shape not in ((2,), (3,)):
            raise ValueError(
                f'{type(self).__name__}.is_translation_symmetric: '
                'not a 2D or 3D vector'
                )
        if translation.shape == (2,):  # 2D. Use zero for z component.
            translation = np.append(translation, 0.)

        # Let's work in the coordinate frame of the unit cell,                  # TODO: change when flipping cartpos[2]
        # i.e., with z axis directed from bulk to surface. This
        # requires updating the translation vector and z_range
        translation[2] *= -1
        if z_range is not None:
            z_range = (self.topat_ori_z - min(z_range),
                       self.topat_ori_z - max(z_range))

        # Use the better version of the unit cell, with a,b,c = uc
        ucell = self.ucell.T
        ucell_inv = np.linalg.inv(ucell)
        releps = eps / np.linalg.norm(ucell, axis=1)

        # Collapse the translation vector toward the origin. This is
        # useful to shortcut the check if the vector is close to zero
        translation, _ = collapse(translation, ucell, ucell_inv, 'round')
        if np.linalg.norm(translation) < eps:
            return True

        # Prepare Cartesian coordinates, collapsed to
        # the base cell, and their translated version
        frac_coords = collapse_fractional(np.array([at.pos for at in self]))
        cart_coords = frac_coords.dot(ucell)
        shifted_coords = cart_coords + translation

        # Unlike in-plane operations, the comparison cannot be
        # done one sublayer at a time. We will compare element
        # by element. Build atom-to-element index mappings for
        # both unshifted and shifted coordinates. Use a list
        # for element_per_atom for now, convert it to array later
        # to use as mask. The shifted one can be an array already
        element_per_atom = [at.el for at in self]
        element_per_atom_shifted = np.array(element_per_atom.copy())

        # Discard shifted atoms that moved out of range in z...
        if not z_periodic:
            if z_range is None:
                z_range = np.min(cart_coords[:, 2]), np.max(cart_coords[:, 2])
            in_range = np.all((shifted_coords[:, 2] >= min(z_range) - eps,
                               shifted_coords[:, 2] <= max(z_range) + eps),
                              axis=0)
            shifted_coords = shifted_coords[in_range]
            element_per_atom_shifted = element_per_atom_shifted[in_range]

        # ...and collapse also the shifted coordinates
        shifted_coords, _ = collapse(shifted_coords, ucell, ucell_inv)

        # Include replicas of atoms close to edges/corners, then
        # convert element_per_atom to array for use it as a mask
        (cart_coords,
         element_per_atom) = add_edges_and_corners(cart_coords,
                                                   frac_coords,
                                                   releps, ucell,
                                                   props=element_per_atom)
        element_per_atom = np.array(element_per_atom)

        # Finally, perform the element-wise comparison
        for element in self.elements:
            distances = euclid_distance(
                shifted_coords[element_per_atom_shifted == element],
                cart_coords[element_per_atom == element],
                )
            if any(distances.min(axis=1) > eps):
                return False
        return True
