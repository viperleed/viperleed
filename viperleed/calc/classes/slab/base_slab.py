"""Module base_slab of viperleed.calc.classes.slab.

Defines the BaseSlab class, useful to describe collections of atoms in
crystalline form. This is the abstract base class at the basis of both
BulkSlab and SurfaceSlab classes (for 3D- and 2D-periodic systems,
respectively), and contains generic functionality that is common to
both.
This module contains refactored and modified functionality that
used to be contained in the original slab module (2019-06-13) by
F. Kraushofer.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-02-21'
__license__ = 'GPLv3+'

from abc import abstractmethod
from collections import Counter
import copy
import itertools
import logging
from numbers import Real
from operator import attrgetter, itemgetter

import numpy as np
from scipy.spatial.distance import cdist as euclid_distance

from viperleed.calc.lib import leedbase
from viperleed.calc.lib.base import COLLAPSE_EPS
from viperleed.calc.lib.base import add_edges_and_corners
from viperleed.calc.lib.base import collapse
from viperleed.calc.lib.base import collapse_fractional
from viperleed.calc.lib.base import pairwise
from viperleed.calc.lib.base import rotation_matrix_order
from viperleed.calc.classes.atom import Atom
from viperleed.calc.classes.atom_containers import AtomContainer, AtomList
from viperleed.calc.classes.layer import Layer, SubLayer
from viperleed.calc.classes.sitetype import Sitetype

from .errors import AlreadyMinimalError
from .errors import EmptySlabError
from .errors import InvalidUnitCellError
from .errors import MissingElementsError
from .errors import MissingLayersError
from .errors import MissingSublayersError
from .errors import SlabError
from .errors import TooFewLayersError
from .utils import _left_handed, _z_distance


_LOGGER_NAME, _ = __name__.rsplit('.', maxsplit=1)
_LOGGER = logging.getLogger(_LOGGER_NAME)


# TODO: .cartpos[2] Issue #174
# TODO: too-many-instance-attributes
#   Perhaps there's a way to group attributes and methods
#   into some logical mix-in base classes?
# TODO: layer coordinates may not be up to date after we update_origin
#   This is also partly related to Issue #174. It would be solved if
#   we (i) use the same frame for atoms, layers, and slab, (ii) we
#   store topat_ori_z only right before a refcalc (with a dedicated
#   method. Perhaps save_refcalc_state), (iii) we remove all of the
#   instances of update_origin.
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
    """

    def __init__(self):
        """Initialize instance."""
        self.ucell = np.array([])
        self.poscar_scaling = 1.
        self.chemelem = set()
        self.n_per_elem = {}
        self.atlist = AtomList()
        self.layers = ()
        self.sublayers = ()
        self.sitelist = []
        self.ucell_mod = []
        self.ucell_ori = np.array([])
        self.topat_ori_z = None                                                 # TODO: SurfaceSlab only after we fix the .cartpos[2] -- Issue #174
        self.celltype = 'unknown'
        self.planegroup = 'unknown'
        self.foundplanegroup = 'unknown'
        self.orisymplane = None
        self.linklists = []
        self.bulkslab = None  # Deleted in BulkSlab.__init__

        # Remember the last value of the ELEMENT_MIX parameter that
        # was applied. Prevents repeated applications
        self._last_element_mix = None

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
                f'{type(self).__name__} has no unit cell defined'
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
    def bottom_atom(self):
        """Return the atom currently at the bottom of this slab."""
        # Do it the most efficient way possible, i.e., with the
        # atom container that should have the fewest atoms. We
        # rely on the z-sorting of layers and sublayers.
        if self.sublayers:
            container = self.sublayers[-1]                                      # TODO: this may not be right, as sublayers are also sorted by element for same-z!
        elif self.layers:
            container = self.layers[-1]
        else:
            container = self
        return min(container, key=lambda at: at.pos[2])

    @property
    def bulk_atoms(self):
        """Return atoms in this Slab that belong to its bulk."""
        return [at for at in self if at.is_bulk]

    @property
    def bulk_layers(self):
        """Return the layers of self that are bulk."""
        return tuple(lay for lay in self.layers if lay.is_bulk)

    @property
    def c_vector(self):
        """Return the out-of-plane unit-cell vector."""
        try:
            return self.ucell.T[2]
        except IndexError:  # Uninitialized
            raise InvalidUnitCellError(
                f'{type(self).__name__} has no unit cell defined'
                ) from None

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
    def interlayer_gaps(self):
        """Return the z distances between pairs of adjacent layers.

        Returns
        -------
        distances : generator of floats
            When iterated over, it returns the z distances between
            the bottommost atom of the higher layer and the topmost
            atom of the lower one. Distances are in the same order
            as self.layers. This means that the first element is the
            gap between the topmost layer and the second layer from
            the top. Notice that there are self.n_layers - 1 items
            in `distances`.

        Raises
        ------
        MissingLayersError
            If no layers are available.
        TooFewLayersError
            If only one layer is present.
        """
        if not self.layers:
            raise MissingLayersError
        if self.n_layers == 1:
            raise TooFewLayersError(f'{type(self).__name__} '
                                    'has only one layer')
        return (lay_below.cartori[2] - lay_above.cartbotz                       # TODO: change when flipping .cartpos[2] -- Issue #174
                for lay_above, lay_below in pairwise(self.layers))

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

    @property
    def non_bulk_layers(self):
        """Return a tuple of the non-bulk layers in this slab."""
        return tuple(lay for lay in self.layers if not lay.is_bulk)

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
        TooFewLayersError
            If only one layer is present.
        """
        return min(self.interlayer_gaps)

    @property
    def top_atom(self):
        """Return the atom currently at the bottom of this slab."""
        # Do it the most efficient way possible, i.e., with the
        # atom container that should have the fewest atoms. We
        # rely on the z-sorting of layers and sublayers.
        if self.sublayers:
            container = self.sublayers[0]                                       # TODO: this may not be right, as sublayers are also sorted by element for same-z!
        elif self.layers:
            container = self.layers[0]
        else:
            container = self
        return max(container, key=lambda at: at.pos[2])

    @classmethod
    def from_slab(cls, other):
        """Return a `cls` instance with attributes deep-copied from `other`."""
        if not isinstance(other, BaseSlab):
            raise TypeError(f'{cls.__name__}.from_slab: other is not a slab.')
        if other.__class__ is cls:
            return copy.deepcopy(other)

        instance = cls()
        memo = {id(other): instance}
        for attr, value in other.__dict__.items():
            if not hasattr(instance, attr):
                # Skip attributes that do not belong to cls instances
                continue
            setattr(instance, attr, copy.deepcopy(value, memo))
        return instance

    def _add_one_bulk_cell(self, bulk_layers, bulkc_par,                        # TODO: Issue #174
                           bulkc_par_atoms, bulkc_perp_atoms,
                           new_atoms_start_index=None):
        """Add one bulk unit cell and return the new atoms.

        This method is intended to be used only internally in
        Slab objects. Use with_extra_bulk_units (for non-bulk
        slabs) or with_double_thickness (for bulk ones).

        Parameters
        ----------
        bulk_layers : tuple of Layer
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
        new_layers : tuple
            The Layer objects added to this slab.

        Raises
        ------
        ValueError
            If `bulk_layers` is empty.
        """
        if not bulk_layers:
            raise ValueError('_add_one_bulk_cell: no bulk_layers to add.')
        if new_atoms_start_index is None:
            new_atoms_start_index = self.n_atoms + 1
        # The mess with needing three different components of the
        # bulk repeat vector comes from the fact that we want to
        # add atoms (and layers) below, but not shift the current
        # ones relative to the current unit cell. This is useful
        # because the 'old atoms' appear then at the same in-plane
        # position in, e.g., VESTA.
        self.c_vector[:] += bulkc_par  # Expand unit cell
        new_atoms = []
        new_layers = []
        for layer in bulk_layers:
            # Add a new bulk layer and duplicate of all its atoms
            new_layer = Layer(self, self.n_layers, is_bulk=True)
            new_layers.append(new_layer)
            self.layers += (new_layer,)
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
        return new_atoms, tuple(new_layers)

    def check_a_b_in_plane(self):
        """Raise InvalidUnitCellError if a, b have out-of-plane components."""
        try:
            a_b_z_component = self.ucell[2, :2]
        except IndexError:
            raise InvalidUnitCellError(f'{type(self).__name__} has no '
                                       'unit cell defined') from None
        if any(a_b_z_component):
            _err = ('Unit cell a and b vectors must not '
                    'have an out-of-surface (Z) component!')
            _LOGGER.error(_err)
            raise InvalidUnitCellError(_err)

    def clear_symmetry_and_ucell_history(self):
        """Set all symmetry information back to default values.

        This method also completely erases the history of unit-cell
        modifications, and sets the current unit cell as the 'original'
        one. This means that, if the unit cell was modified prior to
        a call to this method, the original unit cell **cannot** be
        recovered by a call to `revert_unit_cell`.

        Returns
        -------
        None.
        """
        self.ucell_mod = []
        self.ucell_ori = self.ucell.copy()
        self.celltype = 'unknown'
        self.foundplanegroup = self.planegroup = 'unknown'
        self.orisymplane = None

    def _get_layer_cut_positions(self, rpars, bulk_cuts):
        """Return crystal-cut positions from rpars.LAYER_CUTS."""
        # rpars.LAYER_CUTS is a 'list' of LayerCutTokens. Each either
        # .is_numeric (fractional positions) or .is_auto_cut (dc/dz).
        # We can collect the fractional ones immediately.
        cuts = [token.value for token in rpars.LAYER_CUTS if token.is_numeric]

        # The dc/dz need a bit of work. We use as cut
        # position the midpoints between atoms that are
        # further than the desired z-distance cut-off
        auto_cuts = (token for token in rpars.LAYER_CUTS if token.is_auto_cut)
        _c_to_z = self.c_vector[2] / np.linalg.norm(self.c_vector)
        atoms = sorted(self, key=lambda atom: atom.pos[2])
        for token in auto_cuts:
            cutoff = token.value
            cutoff *= _c_to_z if token.is_auto_cut_dc else 1
            try:
                lower, upper = token.get_bounds(bulk_cuts)
            except ValueError:
                # bulk_cuts > upper. This auto-cut is useless
                continue
            if lower < COLLAPSE_EPS:
                # Avoid mistreating atoms that are very close to c==0
                lower = -2*COLLAPSE_EPS
            cuts.extend(
                sum(at.pos[2] for at in pair) / 2
                for pair in pairwise(atoms)
                if (all(lower < at.pos[2] < upper for at in pair)
                    and _z_distance(*pair) > cutoff)
                )
        if bulk_cuts:
            cuts = [v for v in cuts if v > max(bulk_cuts) + 1e-6] + bulk_cuts
        cuts.sort()
        return cuts

    def create_layers(self, rpars, bulk_cuts=()):
        """Create a list of Layer objects based on `rpars`.

        After this call, the `layers` attribute contains a list of
        the layers created. If layers were already defined, they
        are overwritten.

        Notice that this method assumes that fractional coordinates
        of the atoms are **collapsed to the base cell**.

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS to be used to create layers. The
            N_BULK_LAYERS and LAYER_CUTS attributes are used.
        bulk_cuts : Sequence, optional
            Sequence of fractional positions representing
            automatically detected bulk layer cuts. Default
            is an empty tuple.

        Returns
        -------
        cuts : list
            Fractional coordinates (sorted) of the cuts performed.
            Notice that the number of layers created may be smaller
            than `len(cuts)`, as layers without atoms are deleted.
        """
        self.check_a_b_in_plane()
        if self.is_bulk:
            bulk_cuts = ()

        # Get a sorted list of fractional cut positions
        cuts = self._get_layer_cut_positions(rpars, bulk_cuts)
        n_cuts = len(cuts)  # For checking empty layers later
        cuts_iter = iter(cuts)

        # Remember current sorting to re-sort at the end
        self.atlist.save_sorting()
        self.sort_by_z()           # Bottommost atom first
        next_cut = next(cuts_iter, 1)
        newlayer = Layer(self, 0)  # Layer numbers are fixed later
        self.layers = ()
        for atom in self:
            while atom.pos[2] > next_cut:
                # May run multiple times if the user defined multiple
                # cuts between two atoms that would give empty layers
                next_cut = next(cuts_iter, 1)
                newlayer = Layer(self, 0)
            if newlayer not in self:
                newlayer.is_bulk = (self.is_bulk
                                    or self.n_layers < rpars.N_BULK_LAYERS)
                self.layers += (newlayer,)
            atom.layer = newlayer
            newlayer.atlist.append(atom)

        n_empty = (n_cuts + 1) - self.n_layers
        if n_empty:
            _plural = 's' if n_empty > 1 else ''
            _LOGGER.warning(
                f'Found {n_empty} layer{_plural} containing no atoms. '
                f'Layer{_plural} will be deleted. Check LAYER_CUTS parameter.'
                )
            rpars.setHaltingLevel(2)

        # Sort layers in the final order: top to bottom
        self.layers = tuple(reversed(self.layers))
        for i, layer in enumerate(self.layers):
            layer.update_position()
            layer.num = i
        self.atlist.restore_sorting()
        return cuts

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
            # pylint: disable-next=magic-value-comparison  # The '2'
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
            new_layer.atlist.extend(atoms)
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

        Raises
        ------
        InvalidUnitCellError
            If the (a,b) unit-cell vectors have out-of-plane components
        EmptySlabError
            If there are no atoms in this slab.
        MissingElementsError
            If this slab has no elements. Normally, this means that
            `update_element_count` was not called yet.
        """
        self.check_a_b_in_plane()
        if not self.n_atoms:
            raise EmptySlabError(f'{type(self).__name__} has no atoms')
        if not self.elements:
            raise MissingElementsError(
                f'{type(self).__name__} has no elements. '
                'Did you forget to update_element_count()?'
                )

        self.sort_by_z()
        subl = []
        for element in self.elements:
            subl.extend(self._get_sublayers_for_el(element, eps))

        # subl is sorted element-first. Re-sort it as in the __doc__
        sorted_sublayers = []

        # Work with subl sorted from bottom to top, and pop the last
        # element each time (i.e., the topmost layer to be processed)
        subl.sort(key=attrgetter('cartbotz'), reverse=True)                     # TODO: Issue #174
        while subl:
            same_z = [subl.pop()]  # accumulate sublayers with same z
            this_z = same_z[0].cartbotz
            while subl:
                if abs(subl[-1].cartbotz - this_z) < eps:
                    same_z.append(subl.pop())
                else:
                    break
            same_z.sort(key=attrgetter('element'))
            sorted_sublayers.extend(same_z)
        self.sublayers = tuple(sorted_sublayers)
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

    def collapse_fractional_coordinates(self, releps=None):                     # TODO: would be nicer to have the public method update everything correctly, and rather make this a private one.  USED ONLY IN SLABs and in iosearch
        """Ensure all atoms are inside the unit cell.

        Atoms whose fractional positions are outside the unit cell
        are 'back-folded' inside.

        Notice that calling this method DOES NOT UPDATE the Cartesian
        coordinates of the atoms. It must be preceded or followed by
        a call to `update_cartesian_from_fractional` to ensure that
        fractional and Cartesian coordinates are up to date.

        Parameters
        ----------
        releps : float or Sequence or None, optional
            Fractional tolerance for collapsing coordinates. If
            not None, coordinates are collapsed to the (fractional)
            interval [-releps, 1-releps]. If a sequence, it should
            have three items, corresponding to the fractional
            tolerances along the three unit vectors. Default is None.

        Returns
        -------
        None.

        See Also
        --------
        collapse_cartesian_coordinates
        """
        kwargs = {'in_place': True}
        if releps is not None:
            kwargs['eps'] = releps
        for atom in self:
            collapse_fractional(atom.pos, **kwargs)

    def full_update(self, rpars):
        """Update atoms and layers from `rpars`.

        This method is typically used after information of a
        slab is read from file (e.g., via `poscar.read`).

        The method ensures that all atoms are within the (0, 0) unit
        cell, then, if needed, uses the information from `rpars` to
        create layers and calculate Cartesian coordinates (absolute
        and per layer), and to update elements and sites.
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
        if not self.layers or len(self.bulk_layers) != rpars.N_BULK_LAYERS:
            self.create_layers(rpars)
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

    def get_minimal_ab_cell(self, eps, epsz=None, warn_convention=False):
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
        self_copy = copy.deepcopy(self)
        self_copy.project_c_to_z()
        self_copy.sort_by_z()
        self_copy.create_sublayers(epsz)

        # Use the lowest-occupancy sublayer (the one
        # with fewest atoms of the same chemical element)
        lowocclayer = self_copy.fewest_atoms_sublayer
        n_atoms = lowocclayer.n_atoms
        if n_atoms < 2:  # pylint: disable=magic-value-comparison
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
        candidate_unit_vectors = (p - plist[0] for p in plist[1:])              # TODO: since is_translation_symmetric is somewhat expensive, would it make sense to pre-filter the vectors keeping only the shortest ones among those parallel to one another? Perhaps it would even make sense to filter the candidates keeping only those that give an integer area reduction.
        candidate_unit_vectors = (
            vec for vec in candidate_unit_vectors
            if self_copy.is_translation_symmetric(vec, eps)
            )
        # Try to reduce the cell
        smaller, mincell = self._minimize_ab_area(self.ab_cell.T, n_atoms,
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

        candidate_unit_vectors = itertools.chain(ab_cell,
                                                 candidate_unit_vectors)
        candidate_ab_cells = itertools.combinations(candidate_unit_vectors, 2)
        candidate_ab_cells_and_areas = ((ab, abs(np.linalg.det(ab)))
                                        for ab in candidate_ab_cells)
        candidate_ab_cells_and_areas = (
            (ab, area) for ab, area in candidate_ab_cells_and_areas
            if smallest_area < area + eps**2 < ab_cell_area
            )

        # Finally, consider that the new unit cell area
        # must be an integer divisor of the current one
        candidate_ab_cells_and_area_ratios = (
            (ab, ab_cell_area / area)
            for ab, area in candidate_ab_cells_and_areas
            )
        acceptable_ab_cells_and_ratios = (
            (ab, round(ratio))
            for ab, ratio in candidate_ab_cells_and_area_ratios
            if abs(round(ratio) - ratio) < 1e-6
            )

        try:
            # Pick the one that realizes the largest area reduction
            mincell, _ = max(acceptable_ab_cells_and_ratios, key=itemgetter(1))
        except ValueError:  # No smaller cell
            return False, ab_cell
        return True, mincell

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

    def initSites(self, rp):                                                    # BUG: No sites if called twice
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
            the (0, 0) cell. This means that this slab is considered
            equivalent to other only if the atom positions differ
            by a translation that is a multiple of the unit-cell
            vectors.
        """
        if not isinstance(other, BaseSlab):
            return False

        slabs = copy.deepcopy(self), copy.deepcopy(other)
        for slab in slabs:
            slab.collapse_cartesian_coordinates()
            # Always re-create sublayers as eps may differ
            slab.create_sublayers(eps)
            # Reorder sublayers by Z to then compare by index
            slab.sublayers = tuple(sorted(slab.sublayers,
                                          key=attrgetter('cartbotz')))
        if slabs[0].n_sublayers != slabs[1].n_sublayers:
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
            distances = euclid_distance(other_coords, this_coords)
            if any(distances.min(axis=1) > eps):
                # Notice that we use axis=1 as axis=0 contains
                # this_coords, which has been expanded to add
                # atoms at edges/corners
                return False
        return True

    def project_c_to_z(self):
        """Make the c vector of the unit cell perpendicular to the surface.

        All atom coordinates are updated to fit the new basis.

        Returns
        -------
        None.
        """
        c_vec_xy = self.c_vector[:2]
        if any(c_vec_xy):  # Non-zero components
            self.update_cartesian_from_fractional()
            c_vec_xy[:] = 0
            self.collapse_cartesian_coordinates()  # Also updates fractional

    def remove_duplicate_atoms(self, eps=1e-3, epsz=1e-3, other_slab=None):
        """Remove atoms with same positions and chemical element.

        Should this slab have `layers` defined, they are updated
        accordingly without recreating them. Sublayers are instead
        created from scratch.

        This method explicitly recalculates the correct number of atoms
        per each (POSCAR) element, but does not recalculate Atom.num
        values. Explicitly call self.update_atom_numbers() afterwards,
        if needed.

        Parameters
        ----------
        eps : float, optional
            Cartesian tolerance (Angstrom) for in-plane
            equivalence. Default is 1e-3.
        epsz : float, optional
            Cartesian tolerance (Angstrom) for out-of-plane
            equivalence. Default is 1e-3.
        other_slab : Slab or None, optional
            If given and not None, atoms in `other_slab` that are
            removed from this slab are labelled as being duplicates
            of those that remain in this one. Atoms between the two
            slabs are compared by `num`. Default is None.

        Raises
        ------
        SlabError
            If this slab has layers defined, but removing atoms
            led to an inconsistency of which atoms survived in
            self and which ones survived in its layers.
        """
        try:
            other_slab_atlist = other_slab.atlist
        except AttributeError:
            other_slab_atlist = AtomList()

        # The easiest way to remove duplicates is going one sublayer
        # at a time, and check only the in-plane positions. Ensure we
        # have sublayers, and, for that, we also make sure that the
        # number-of-atoms-per-element is up to date.
        self.update_element_count()
        self.create_sublayers(epsz)
        new_atlist = AtomList()
        for sublayer in self.sublayers:
            i = 0
            while i < sublayer.n_atoms:
                atom_i = sublayer.atlist[i]
                other_base_atom = other_slab_atlist.get(atom_i.num, None)
                j = i + 1
                while j < sublayer.n_atoms:
                    atom_j = sublayer.atlist[j]
                    if not atom_i.is_same_xy(atom_j, eps=eps):
                        j += 1
                        continue
                    sublayer.atlist.pop(j)
                    duplicate_atom = other_slab_atlist.get(atom_j.num, None)
                    if duplicate_atom:
                        duplicate_atom.duplicate_of = other_base_atom
                i += 1
            new_atlist.extend(sublayer.atlist)
        self.atlist = new_atlist

        # Update number of atoms per element and sublayers again
        self.update_element_count()
        self.create_sublayers(epsz)

        # Finally, update the layers, without recreating them
        if not self.layers:
            return
        self._delete_removed_atoms_from_layers()

    def _delete_removed_atoms_from_layers(self):
        """Remove atoms not in self from layers. Keep only filled layers."""
        # First remove the atoms that are not in self any longer
        for layer in self.layers:
            atoms = [at for at in layer if at in self]
            layer.atlist.clear()
            layer.atlist.extend(atoms)

        # Then keep only non-empty layers. Some layers may
        # be completely empty, e.g., when detecting bulk
        self.layers = tuple(layer for layer in self.layers
                            if layer.n_atoms)
        for i, layer in enumerate(self.layers):
            layer.update_position()
            layer.num = i

        # Sanity check: we should now have the same atoms everywhere
        if set(self) != set(at for lay in self.layers for at in lay):
            raise SlabError(
                'Inconsistent atoms in slab and its layers. Did you forget '
                'to appropriately create_layers and update rparams before?'
                )

    def remove_vacuum_at_bottom(self, rpars):
        """Move all atoms along c to remove any vacuum at the bottom.

        This method **will update the origin of the z coordinates**.
        It should be called **before** a reference calculation.

        Parameters
        ----------
        rpars : Rparams
            The current run parameters. Attributes read: SYMMETRY_EPS

        Returns
        -------
        fractional_c_shift : float
            How much atoms have been moved down along c.

        Notes
        -----
        The atoms are moved down slightly less than what would be
            necessary to bring the bottom-most atom(s) at c == 0
            to avoid possible floating-point rounding issues when
            collapsing coordinates.
        This method does not check whether the coordinates are
            collapsed to the base cell. Make sure they are before
            calling this.
        """
        # Skip movements if the bottom atom is closer to zero than eps
        eps = (rpars.SYMMETRY_EPS, rpars.SYMMETRY_EPS, rpars.SYMMETRY_EPS.z)
        releps = np.dot(eps, np.linalg.inv(self.ucell.T))
        bottom_atom = self.bottom_atom
        if bottom_atom.pos[2] < releps[2]:
            return 0.

        delta_c = bottom_atom.pos[2] - 0.1 * releps[2]
        for atom in self:
            atom.pos[2] -= delta_c
        self.update_cartesian_from_fractional(update_origin=True)
        return delta_c

    def revert_unit_cell(self, restore_to=None):
        """Revert unit-cell and coordinate modifications.

        Notice that calling this method always collapses the atom
        coordinates to the base unit cell. This is independent of
        whether any unit-cell modification was actually reverted
        or not.

        This method will not revert the state of, e.g., layers.
        If an operation that can (potentially) change layers is
        reverted, the layers and sublayers attributes are cleared.

        Parameters
        ----------
        restore_to : Sequence or None, optional
            If given and not None, keep the oldest `len(restore_to)`
            modifications of the unit cell and coordinates. Default
            is None.

        Raises
        -------
        RuntimeError
            If any of the to-be-reverted operations stored in
            `slab.ucell_mod` has an unknown operation-type name.
        """
        n_keep = 0
        if restore_to is not None:
            n_keep = len(restore_to)

        self.update_cartesian_from_fractional()
        operations_to_undo = self.ucell_mod[n_keep:]
        for op_type, op_array in reversed(operations_to_undo):
            if op_type == 'add':
                frac_shift = np.dot(op_array, np.linalg.inv(self.ab_cell.T))
                for atom in self:
                    atom.translate_2d(-op_array, -frac_shift)
            elif op_type == 'c_shift':
                for atom in self:
                    atom.pos[2] -= op_array
                    collapse_fractional(atom.pos, in_place=True)
                self.update_cartesian_from_fractional(update_origin=True)
                self.layers = ()     # They're probably wrong
                self.sublayers = ()  # Probably wrong sorting
            elif op_type == 'lmul':
                self.ucell = np.dot(np.linalg.inv(op_array), self.ucell)
            elif op_type == 'rmul':
                self.ucell = np.dot(self.ucell, np.linalg.inv(op_array))
            else:
                raise RuntimeError(f'Invalid unit-cell modification {op_type}')
        self.collapse_cartesian_coordinates()
        self.ucell_mod = self.ucell_mod[:n_keep]

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
                f'update_element_count()? elements={self.elements}'
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

    def update_cartesian_from_fractional(self, update_origin=False):            # TODO: Issue #174
        """Assign absolute Cartesian coordinates to all atoms.

        The frame of reference has (x, y) as defined by the a and b
        unit-cell vectors, z == 0 for the topmost atom and positive
        going down through the slab (from the surface into the solid).

        This method can be used, e.g., to (re)calculate all Cartesian
        coordinates when distortions have been applied to the unit
        cell (but atoms should retain their position relative to the
        modified cell).

        Notice that this method does not touch the atomic coordinates
        in this slab's `bulkslab`.

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
            atom.cartpos[2] = self.topat_ori_z - atom.cartpos[2]                # TODO: Remove with Issue #174

    def _update_chem_elements(self, rpars):                                     # TODO: Why aren't we also taking into account ELEMENT_RENAME here? Issue #108
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
        try:
            self.bulkslab.update_element_count()
        except AttributeError:  # BulkSlab or no bulkslab available
            pass

    def update_fractional_from_cartesian(self):
        """Calculate atoms' fractional coordinates from Cartesian ones.

        This method is commonly used when `ucell` is modified, and
        the atoms should retain their absolute Cartesian positions
        while recomputing their fractional coordinates relative to
        the new unit cell.

        Notice that this method does not touch the atomic coordinates
        in this slab's `bulkslab`.

        Returns
        -------
        None.
        """
        uci = np.linalg.inv(self.ucell)
        for atom in self:
            # Flip over the z Cartesian coordinate, as
            # we store it as 'positive going down'
            cartpos = atom.cartpos.copy()
            cartpos[2] = self.topat_ori_z - cartpos[2]                          # TODO: Remove when flipping .cartpos[2] -- Issue #174
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
            self.layers = ()
            self.sublayers = ()

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
                            f'one or three, got {len(scaling)}')
        if not all(isinstance(s, Real) for s in scaling):
            raise TypeError(f'{type(self).__name__}.apply_scaling: '
                            f'invalid scaling factor. Expected one '
                            'or three numbers')
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

        if on_row_vectors:                                                      # TODO: this has to be turned around when switching ucell.T -- Issue #175
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

    def translate_atoms_2d(self, shift):
        """Add a 2D Cartesian shift to all atomic coordinates.

        Parameters
        ----------
        shift : Sequence
            Two elements. The 2D Cartesian rigid shift to apply
            to all atoms.

        Notes
        -----
        This method **does not collapse** the atom coordinates to
            the base unit cell. Explicitly follow a call to this
            method with .collapse_cartesian_coordinates() if the
            atoms should be back-folded into the base unit cell.
        This operation can be undone with a call to `revert_unit_cell`,
            as the shift is registered as one of the `ucell_mod`.
        """
        shift = np.asarray(shift)
        # Since fractional positions and Cartesian positions are
        #    frac = cart @ ab_inv,
        # we can move Cartesian by shift, and fractional by:
        frac_shift = shift.dot(np.linalg.inv(self.ab_cell.T))
        for atom in self:
            atom.translate_2d(shift, frac_shift)
        self.ucell_mod.append(('add', shift))

    def translate_atoms_c(self, c_fraction):
        """Move all atoms up by `c_fraction` along the c vector.

        Parameters
        ----------
        c_fraction : float
            The rigid shift in fractional coordinates to be
            added to all the atoms.

        Notes
        -----
        At variance with other atom- and unit-cell-transforming
            operations, this method **does collapse** the atom
            coordinates to the base unit cell and **updates the
            stored position of the topmost atom**.
        As translating atoms along c has the potential to modify
            layers and sublayers, the `.layers` ans `.sublayers`
            attributes are always cleared. Call `create_(sub)layers`
            again to recreate them.
        This operation can be undone with a call to `revert_unit_cell`,
            as the shift is registered as one of the `ucell_mod`. The
            clearing of layers cannot be automatically undone.
        """
        for atom in self:
            atom.pos[2] += c_fraction
            collapse_fractional(atom.pos, in_place=True)
        self.update_cartesian_from_fractional(update_origin=True)
        self.ucell_mod.append(('c_shift', c_fraction))
        self.layers = ()
        self.sublayers = ()

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

        # Let's work in the coordinate frame of the unit cell,                  # TODO: change when flipping .cartpos[2] -- Issue #174
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
        # convert element_per_atom to array to use it as a mask
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
