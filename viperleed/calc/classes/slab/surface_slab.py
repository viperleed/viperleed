"""Module surface_slab of viperleed.calc.classes.slab.

Defines the SurfaceSlab class, a BaseSlab describing a 2D-periodic
crystal. This is available in viperleed.calc as the Slab alias for
backwards compatibility with the former module-based structure of
slab.py by @fkraushofer (originally created on 2019-06-13).
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)'
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-02-21'
__license__ = 'GPLv3+'

from collections import Counter
import copy
import itertools
import logging
from operator import itemgetter

import numpy as np
from scipy.spatial import KDTree, distance as sp_distance

from viperleed.calc.lib import leedbase
from viperleed.calc.lib.base import NonIntegerMatrixError
from viperleed.calc.lib.base import SingularMatrixError
from viperleed.calc.lib.base import add_edges_and_corners
from viperleed.calc.lib.base import collapse
from viperleed.calc.lib.base import collapse_fractional
from viperleed.calc.lib.base import ensure_integer_matrix
from viperleed.calc.lib.base import pairwise
from viperleed.calc.classes.atom import Atom
from viperleed.calc.classes.atom_containers import AtomList
from viperleed.calc.files.parameters.errors import (
    InconsistentParameterError
    )
from viperleed.calc.lib.periodic_table import PERIODIC_TABLE, COVALENT_RADIUS

from .base_slab import BaseSlab
from .bulk_slab import BulkSlab
from .errors import AlreadyMinimalError
from .errors import AtomsTooCloseError
from .errors import EmptySlabError
from .errors import MissingBulkSlabError
from .errors import MissingLayersError
from .errors import NoBulkRepeatError
from .errors import NotEnoughVacuumError
from .errors import NoVacuumError
from .errors import TooFewLayersError
from .errors import SlabError
from .errors import WrongVacuumPositionError

try:
    import ase
except ImportError:
    _HAS_ASE = False
else:
    _HAS_ASE = True


_LOGGER_NAME, _ = __name__.rsplit('.', maxsplit=1)
_LOGGER = logging.getLogger(_LOGGER_NAME)
_MIN_VACUUM = 5.0   # Angstrom
_VACUUM_EPS = 1e-4  # In fractional c coordinates


class SurfaceSlab(BaseSlab):
    """A class representing a semi-infinite solid.

    Contains unit cell, element information and atom coordinates.
    Also has a variety of convenience functions for manipulating
    and updating the atoms. Slabs can be created from an ase.Atoms
    object using the `from_ase` class method. Another common way
    is using `viperleed.calc.files.poscar.read`.

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
    displists : list of list of Atom
        List of lists of atoms which are displaced together. This
        differs from `linklists` in that while `linklists` stores
        the linking based on the current symmetry, the `displists`
        store actual displacement linking based on the symmetry set
        at the time the displacement is assigned.
    deltas_initialized : bool
        Set by Rparams.generateSearchPars
    preprocessed : bool
        True if the POSCAR that this slab was read from had the
        'Plane group = XY' comment in the header, indicating that
        it was processed by viperleed.calc before.
    symbaseslab : Slab or None
        Slab with the smallest in-plane unit-cell area that shows
        the full symmetry of the slab.
    bulkslab : Slab or None
        Slab object containing only bulk layers
    """

    def __init__(self):
        """Initialize instance."""
        super().__init__()
        self.displists = []
        self.preprocessed = False
        self.deltas_initialized = False
        self.symbaseslab = None

    @property
    def is_bulk(self):
        """Return whether this is a bulk slab."""
        return False

    @property
    def thickness(self):
        """Return the thickness of this slab along z."""
        top_atom_z = min(at.cartpos[2] for at in self)                          # TODO: swap min and max with .cartpos[2] -- Issue #174
        bottom_atom_z = max(at.cartpos[2] for at in self)
        return abs(top_atom_z - bottom_atom_z)

    @property
    def vacuum_gap(self):
        """Return the z distance that does not not contain atoms."""
        return self.c_vector[2] - self.thickness

    @classmethod
    def from_ase(cls, ase_atoms, sort_elements=True):
        """Return a slab initialized from an ase.Atoms object.

        Parameters
        ----------
        ase_atoms : ase.Atoms
            The Atomic-simulation-environment Atoms object to
            initialize from.
        sort_elements : bool, optional
            Whether chemical elements should be sorted in
            alphabetical order before storing. Useful to
            prevent incorrect ordering when an existing
            PHASESHIFTS file is used. Default is True.

        Returns
        -------
        slab : Slab
            The Slab, initialized from ase_atoms. May be empty
            if the ase module is not available.

        Raises
        ------
        ModuleNotFoundError
            If the ase module is not available.
        TypeError
            If `ase_atoms` is not an ase.Atoms object.
        """
        if not _HAS_ASE:
            _err_msg = ('Slab creation from ase.Atoms not '
                        'supported: Module ase not found')
            _LOGGER.error(_err_msg)
            raise ModuleNotFoundError(_err_msg, name='ase', path=__file__)
        if not isinstance(ase_atoms, ase.Atoms):
            raise TypeError(f'{cls.__name__}.from_ase requires an ase.Atoms '
                            f'object, found {type(ase_atoms)!r} instead.')
        self = cls()
        self.ucell = ase_atoms.cell.copy().T
        elems = ase_atoms.symbols
        n_per_el = dict(Counter(elems))
        if sort_elements:
            n_per_el = {k: n_per_el[k] for k in sorted(n_per_el.keys())}
        self.n_per_elem = n_per_el
        for elem, pos in zip(elems, ase_atoms.get_scaled_positions()):
            self.atlist.append(Atom(elem, pos, self.n_atoms + 1, self))
        self.update_cartesian_from_fractional()
        return self

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
            A new slab instance with attributes copied from `other`.
            `new_slab` will not have .layers and .sublayers defined
            if other.is_bulk. This is to prevent `new_slab` from
            having too many (i.e., more than RParams.N_BULK_LAYERS)
            bulk layers.
        """
        new_slab = super().from_slab(other)
        if other.is_bulk:
            new_slab.layers = ()
            new_slab.sublayers = ()
        return new_slab

    # Used only in ensure_minimal_bulk_ab_cell
    def _change_bulk_cell(self, rpars, new_ab_cell):
        """Assign a new 2D unit cell to the bulk, if possible.

        The new cell is applied only if it gives an acceptable
        `SUPERLATTICE` matrix. Otherwise exceptions are raised.

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS to be used and updated. Attributes accessed:
            (read/write) SUPERLATTICE; (read) BULK_REPEAT, SYMMETRY_EPS,
            superlattice_defined.
        new_ab_cell : Sequence
            Shape (2, 2), representing the new 2D unit cell in
            Cartesian coordinates. The new unit vectors should
            be rows. A tentative SUPERCELL is calculated from it.

        Raises
        ------
        NonIntegerMatrixError
            If `new_ab_cell` gives a superlattice matrix that is not
            integer-valued.
        calc.files.parameters.errors.InconsistentParameterError
            If `new_ab_cell` gives a different superlattice matrix
            than the one in `rpars`.
        """
        # Calculate new SUPERLATTICE matrix
        superlattice = self.ab_cell.T.dot(np.linalg.inv(new_ab_cell))
        try:
            superlattice = ensure_integer_matrix(superlattice, eps=5e-2)
        except NonIntegerMatrixError:
            raise NonIntegerMatrixError(
                'Automatically detected bulk SUPERLATTICE is '
                f'not integer-valued: \n{superlattice}'                         # TODO: guilib array formatter
                ) from None
        if (rpars.superlattice_defined
                and not np.allclose(superlattice, rpars.SUPERLATTICE)):
            raise InconsistentParameterError(
                'Automatically detected minimum-area bulk unit cell differs '
                'from the cell defined by the SUPERLATTICE parameter. '
                'Consider changing the SUPERLATTICE parameter. Found matrix: '
                f'\n{superlattice.astype(int)}'                                 # TODO: guilib array formatter
                )
        rpars.SUPERLATTICE = superlattice

        if not self.bulkslab:
            raise MissingBulkSlabError
        # Notice that there is no need to re-center along bulk c, as
        # we're only modifying the in-plane cell. There's also no
        # reason to modify the c coordinates, hence z_periodic=False
        self.bulkslab.apply_bulk_cell_reduction(rpars.SYMMETRY_EPS,
                                                epsz=rpars.SYMMETRY_EPS.z,
                                                new_ab_cell=new_ab_cell,
                                                recenter=False,
                                                z_periodic=False)

    def check_atom_collisions(self, eps=0.05):
        """Raise if any pairs of Atoms are too close to one another.

        Parameters
        ----------
        eps : float, optional
            Minimum Cartesian distance (Angstrom) for Atoms to be
            considered colliding. Default is 0.05.

        Raises
        ------
        MissingLayersError
            If this method is called before this Slab has layers.
        AtomsTooCloseError
            If any pair of atoms, including 2D replicas, are too
            close to one another.
        """
        if not self.layers:
            raise MissingLayersError('Cannot check_atom_collisions '
                                     'without layers')

        ucell = self.ucell.T
        ucell_inv = np.linalg.inv(ucell)
        releps = eps / np.linalg.norm(self.ab_cell.T, axis=1)

        # We will check atoms including 2D replicas. We should take
        # atoms from two adjacent layers if the layers are close
        # enough (most of the times this will not be the case).
        layers = (*self.layers, None)
        for upper_layer, lower_layer in pairwise(layers):
            cartesians = [at.cartpos for at in upper_layer]
            try:
                layer_dist = abs(lower_layer.cartori[2] - upper_layer.cartbotz)
            except AttributeError:  # lower_layer is None
                layer_dist = 2*eps  # Just to skip the next check
            if layer_dist < eps:
                cartesians.extend(at.cartpos for at in lower_layer)
            cartesians = np.array(cartesians)
            cartesians[:, 2] = self.topat_ori_z - cartesians[:, 2]              # TODO: .cartpos[2]. We can entirely skip this when we flip -- Issue #174
            cartesians, fractionals = collapse(cartesians, ucell, ucell_inv)
            cartesians, _ = add_edges_and_corners(cartesians, fractionals,
                                                  releps, ucell)
            too_close = sp_distance.pdist(cartesians, 'euclidean') <= eps
            if any(too_close):
                raise AtomsTooCloseError(
                    f'{type(self).__name__} contains atoms closer than {eps} '
                    'angstrom. This is problematic for a LEED calculation. '
                    'Are there duplicate atoms? Remove duplicates, e.g., '
                    f'via slab.remove_duplicate_atoms(eps={eps}, epsz={eps})'
                    )

    def check_vacuum_gap(self):
        """Complain if there's not enough vacuum or it is not at the top.

        An interatomic z distance of at least 5 angstrom is considered
        an acceptable vacuum gap if it is between the topmost atom and
        the (periodic replica of the) bottommost atom.

        This slab is modified only if there is a large-enough vacuum
        gap very close to the top, i.e., if the atom right above the
        vacuum gap is closer to c=1 than 1e-4.

        Raises
        ------
        NoVacuumError
            If the vacuum gap found is extremely small, i.e., if
            the top and bottom atoms are less than 2e-4 apart along
            the c vector.
        NotEnoughVacuumError
            If a small vacuum gap is found (< 5 angstrom).
        WrongVacuumPositionError
            If there is enough vacuum, but it is somewhere in
            the middle of the slab rather than at the top.
        """
        # Work with a copy not to mess around sorting and positions.
        # This will also become the vacuum-fixed version of self.
        self_copy = copy.deepcopy(self)

        # Ensure the c fractional coordinates are in [0, 1).
        # Also keep the Cartesians up to date, as we will look
        # for the pair of atoms that are furthest along z.
        self_copy.collapse_cartesian_coordinates(update_origin=True)
        self_copy.sort_by_z()  # From bottom to top

        # Collect all the Cartesian coordinates. Also include
        # a periodic replica of the bottom-most atom at the top
        z_cartpos = [self_copy.topat_ori_z - at.cartpos[2] for at in self_copy]
        replica_z = self_copy.topat_ori_z - self_copy.atlist[0].cartpos[2]        # TODO: .cartpos[2] -- Issue #174
        replica_z += self_copy.c_vector[2]  # Don't care about in-plane
        z_cartpos.append(replica_z)

        # Notice that taking the max on both the distance and the index
        # ensures that ind_above is always the largest index (i.e., the
        # highest atom) that realizes a certain vacuum gap. This is
        # also the reason for the round(..., 5): avoids float errors
        vacuum_gap, ind_above = max(
            ((round(c_above - c_below, 5), i+1)
             for i, (c_below, c_above) in enumerate(pairwise(z_cartpos)))
            )
        already_at_top = ind_above == self.n_atoms
        if vacuum_gap >= _MIN_VACUUM:
            bottom_atom = self_copy.atlist[0 if already_at_top else ind_above]
            info_msg = (
                f'The slab has a {vacuum_gap:.2f} A vacuum gap above '
                f'{self_copy.atlist[ind_above-1]}, which is not at the top. '
                )
            self._move_vacuum_to_the_top(bottom_atom, already_at_top,
                                         self_copy, info_msg)
            return

        # Not enough vacuum. Assume top/bottom atoms are at top/bottom
        current_vacuum = (z_cartpos[-1]     # Replica of bottom atom
                          - z_cartpos[-2])  # The actual top one
        self._complain_about_small_vacuum(self_copy, current_vacuum)

    @staticmethod
    def _complain_about_small_vacuum(fixed_slab, current_vacuum):
        """Raise VacumError(s) for a slab with too little vacuum."""
        # Assume that the top/bottom atoms are the actual top/bottom.
        # Remember whether the gap is very very small.
        top_atom = fixed_slab.atlist[-1]
        bot_atom = fixed_slab.atlist[0]
        no_vacuum = (bot_atom.pos[2] - top_atom.pos[2]) % 1.0 < 2*_VACUUM_EPS

        # Produce a fixed version: Expand c vector, then ensure
        # that the bottom atom is at least at c==_VACUUM_EPS. Also,
        # give a bit more vacuum than needed to avoid floating-point
        # precision issues. We don't fixed_slab to cause complaints.
        extra_vacuum = _MIN_VACUUM + _VACUUM_EPS - current_vacuum
        c_extra = fixed_slab.c_vector * extra_vacuum / fixed_slab.c_vector[2]
        fixed_slab.c_vector[:] += c_extra
        fixed_slab.update_fractional_from_cartesian()
        if bot_atom.pos[2] < _VACUUM_EPS:
            fixed_slab.translate_atoms_c(_VACUUM_EPS - bot_atom.pos[2])
        err_msg = (
            f'The slab has a {current_vacuum:.2f} A vacuum gap, which is '
            f'smaller than the minimum ({_MIN_VACUUM} A). Will assume '
            f'that {bot_atom} is at the bottom and {top_atom} at the top'
            )
        exc = NoVacuumError if no_vacuum else NotEnoughVacuumError
        raise exc(err_msg, fixed_slab)

    def _move_vacuum_to_the_top(self, bottom_atom, already_at_top,
                                fixed_slab, info_msg):
        """Move atoms along c, if possible, to have vacuum at the top.

        After a successful call to this method, self has its bottom
        atom at c >= 1e-4. This method **does not** check that the
        vacuum gap is large enough!

        Parameters
        ----------
        bottom_atom : Atom
            The atom that was detected as the one right above vacuum.
        already_at_top : bool
            Whether the vacuum gap of self is already in the top
            part of the slab, i.e., there is no atom above it.
        fixed_slab : SurfaceSlab
            A copy of self to be used for producing a vacuum-
            corrected version. Used for error reporting.
        info_msg : str
            Partial error message to be reported in case this slab
            cannot be automatically corrected (i.e., the vacuum is
            somewhere in the middle).

        Raises
        ------
        WrongVacuumPositionError
            If the vacuum gap is somewhere in the middle of self.
        """
        c_shift = collapse_fractional(-bottom_atom.pos[2], eps=_VACUUM_EPS)
        if abs(c_shift) > _VACUUM_EPS and already_at_top:
            # Vacuum is already at the top, and the bottom atom is
            # far enough from the bottom edge. Leave all as is.
            return

        info_msg += f' The bottom atom is {bottom_atom}, at c={-c_shift:.5f}'
        if abs(c_shift) < _VACUUM_EPS:
            # Vacuum is close to the top, but the bottom atom is very
            # close to the c edge. Move it up, at c == _VACUUM_EPS
            c_shift += _VACUUM_EPS
            _LOGGER.info(info_msg + '. Atoms were shifted along c (by '
                         f'{c_shift:.5f}) so that {bottom_atom} is at the '
                         f'bottom (at c={_VACUUM_EPS:.5f}). PARAMETERS (e.g., '
                         'LAYER_CUTS, BULK_LIKE_BELOW, BULK_REPEAT, ...) '
                         'may need adjustment.')
            self.translate_atoms_c(c_shift)
            return

        # The vacuum is somewhere in the middle. Complain.
        fixed_slab.translate_atoms_c(c_shift + _VACUUM_EPS)
        raise WrongVacuumPositionError(info_msg, fixed_slab)

    def detect_bulk(self, rpars, second_cut_min_spacing=1.2):
        """Determine the minimal bulk repeat vector from BULK_LIKE_BELOW.

        Notice that this method modifies `rpars.BULK_REPEAT` and
        `rpars.SUPERLATTICE` so they contain the detected (minimal)
        bulk-repeat vector and the superlattice matrix. It also
        modifies `rpars.LAYER_CUTS` and `rpars.N_BULK_LAYERS` to
        reflect the newly detected bulk. Finally, this slab will
        have its `.layers` and `bulkslab` set accordingly.

        Parameters
        ----------
        rpars : Rparams
            Run parameters object. Attributes accessed:
            (read)  BULK_LIKE_BELOW, SYMMETRY_EPS,
                    SUPERLATTICE, superlattice_defined
            (write) BULK_REPEAT, LAYER_CUTS, N_BULK_LAYERS,
                    SUPERLATTICE.
        second_cut_min_spacing : float, optional
            Minimal z spacing in Angstrom between sublayers of the
            detected bulk unit to make an additional cut through
            the bulk. If no two sublayers are further apart then
            this threshold, the bulk will be treated as a single
            layer. Default is 1.2.

        Returns
        -------
        bulk_cuts : list of float
            Layer cuts for self that generate the bulk.
        bulk_dist : float
            The largest distance found between bulk sublayers. The
            second bulk cut position is placed there if `bulk_dist`
            is larger than `second_cut_min_spacing`.

        Raises
        ------
        ValueError
            If no `BULK_LIKE_BELOW` was defined in `rpars` or if
            `BULK_REPEAT` is already defined.
        NoBulkRepeatError
            If no repeat vector is found.
        """
        if rpars.BULK_LIKE_BELOW <= 0:
            raise ValueError(f'{type(self).__name__}.detect_bulk: rpars '
                             'must have a positive BULK_LIKE_BELOW defined')
        if rpars.BULK_REPEAT is not None:
            raise ValueError(f'{type(self).__name__}.detect_bulk: rpars '
                             'already has a BULK_REPEAT defined')

        # Work with a temporary copy of this slab to mess up layers,
        # and a temporary copy of rpars to modify LAYER_CUTS, and
        # N_BULK_LAYERS. This copy will also store the correct
        # BULK_REPEAT and SUPERLATTICE till we're sure that the
        # procedure is successful
        self_copy = copy.deepcopy(self)
        rpars_copy = copy.deepcopy(rpars)
        rpars_copy.LAYER_CUTS.update_from_sequence([rpars.BULK_LIKE_BELOW])
        rpars_copy.N_BULK_LAYERS = 1
        self_copy.create_layers(rpars_copy)

        # Make sure that the bottom-most atom is at c=0. This is
        # important to prevent the topmost bulk atoms from being
        # back-folded to the bottom when creating the 'thick' bulk
        # slab below.
        bottom_atom_c = self_copy.remove_vacuum_at_bottom(rpars_copy)

        # Create a pseudo-bulk slab to determine the correct repeat
        # c vector: the c vector now is very likely to be wrong (as
        # it is just chopped off from the one of self). Notice that
        # here we should not re-center the coordinates, as we know
        # that the bulk repeat is incorrect.
        self_copy.make_bulk_slab(rpars_copy, recenter=False)

        # Reduce in-plane bulk. This also updates .SUPERLATTICE
        self_copy.ensure_minimal_bulk_ab_cell(rpars_copy)

        # Detect new c vector (z_periodic=False because current c is
        # the one of self, and is most likely wrong), and collapse cell
        # in the process to get the proper bulk slab. Notice that the
        # next call also saves the new c vector in .BULK_REPEAT
        try:
            self_copy.bulkslab.ensure_minimal_c_vector(rpars_copy,
                                                       z_periodic=False)
        except AlreadyMinimalError as exc:
            _LOGGER.error('Automatic bulk detection failed: Found no bulk '
                          'repeat vector below the specified cut-off.')
            raise NoBulkRepeatError('Failed to detect '
                                    'bulk repeat vector') from exc

        # Identify the cut positions that give bulk layers
        # pylint: disable-next=protected-access
        bulk_cuts, bulk_dist = self_copy._get_bulk_cuts(rpars.SYMMETRY_EPS.z,
                                                        second_cut_min_spacing)
        bulk_cuts = [b + bottom_atom_c for b in bulk_cuts]

        # Finally, update self and rpars with the new information
        rpars.BULK_REPEAT = rpars_copy.BULK_REPEAT
        rpars.SUPERLATTICE = rpars_copy.SUPERLATTICE
        rpars.N_BULK_LAYERS = len(bulk_cuts)
        new_cuts = self.create_layers(rpars, bulk_cuts=bulk_cuts)               # TODO: Issue #121
        rpars.LAYER_CUTS.update_from_sequence(new_cuts)
        self.make_bulk_slab(rpars)
        return bulk_cuts, bulk_dist

    def _get_bulk_cuts(self, epsz, second_cut_min_spacing):
        """Return cut positions for self to give bulk layers."""
        # Calculate cut plane (fractional) between bulk and non-bulk
        # portion of self: will take atoms in the bottommost non-bulk
        # layer, and those in the topmost bulk layer. However, take
        # into account that we may not yet have self.layers defined
        # (nor their bulk or non-bulk nature)
        frac_atoms_z = [at.pos[2] for at in self]

        bulk_height = abs(self.bulkslab.c_vector[2])
        frac_bulk_thickness = bulk_height / abs(self.c_vector[2])
        frac_lowest_pos = min(frac_atoms_z)
        frac_bulk_onset = (frac_lowest_pos + frac_bulk_thickness
                           - (epsz / self.c_vector[2]))
        slab_cuts = [(max(f for f in frac_atoms_z if f < frac_bulk_onset)
                     + min(f for f in frac_atoms_z if f > frac_bulk_onset))
                     / 2]

        if self.bulkslab.n_sublayers == 1:
            return slab_cuts, 0.

        # Look for second cut, where bulk sublayers are furthest apart
        lays_and_distances = (
            (lay_above, abs(lay_above.cartbotz - lay_below.cartbotz))
            for lay_above, lay_below in pairwise(self.bulkslab.sublayers)
            )
        lay_above_cut, maxdist = max(lays_and_distances, key=itemgetter(1))
        if maxdist >= second_cut_min_spacing:
            # Get fractional position, in bulk cell, of the cut
            bulk_frac_cut_from_lowest = (
                lay_above_cut.pos[2]
                - (0.5 * maxdist) / bulk_height
                - min(lay.pos[2] for lay in self.bulkslab.sublayers)
                )
            slab_cuts.append(
                frac_lowest_pos
                + bulk_frac_cut_from_lowest * frac_bulk_thickness
                )
        return slab_cuts, maxdist

    def ensure_minimal_bulk_ab_cell(self, rpars, warn_convention=False):
        """Make sure the `bulkslab` has the smallest possible in-plane cell.

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS to be used and updated. Attributes accessed:
            (read/write) SUPERLATTICE; (read) BULK_REPEAT, SYMMETRY_EPS,
            superlattice_defined; setHaltingLevel method
        warn_convention : bool, optional
            Whether warnings should be logged in case the reduced cell
            could not be made to follow standard conventions.

        Raises
        ------
        NonIntegerMatrixError
            If reduction to a minimal cell would give a SUPERLATTICE
            matrix with non-integer values.
        calc.files.parameters.errors.InconsistentParameterError
            If reduction to a minimal cell would give a SUPERLATTICE
            different from the one in `rpars`.
        """
        bulk = self.bulkslab or self.make_bulk_slab(rpars)
        eps = rpars.SYMMETRY_EPS
        try:
            mincell = bulk.get_minimal_ab_cell(eps, eps.z, warn_convention)
        except AlreadyMinimalError:
            return
        self._change_bulk_cell(rpars, mincell)

    def get_bulk_repeat(self, rpars, only_z_distance=False):
        """Return the bulk repeat vector (with positive z).

        Notice that this method does **not attempt to identify an
        unknown bulk-repeat vector**. If that's what you're after,
        use `identify_bulk_repeat()` instead.

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

        Raises
        ------
        TooFewLayersError
            If rpars.BULK_REPEAT is None, and this slab has either
            no bulk layers or no non-bulk layers.
        """
        if isinstance(rpars.BULK_REPEAT, np.ndarray):
            bulkc = rpars.BULK_REPEAT.copy()
            if bulkc[2] < 0:
                bulkc *= -1
            return bulkc

        zdiff = (rpars.BULK_REPEAT if rpars.BULK_REPEAT is not None
                 else self._get_bulk_to_non_bulk_distance())
        if only_z_distance:
            return zdiff
        return self.c_vector * zdiff / self.c_vector[2]

    def _get_bulk_to_non_bulk_distance(self):
        """Return the z distance from bottommost bulk and non-bulk layers."""
        bulk_layers = self.bulk_layers
        non_bulk_layers = self.non_bulk_layers
        if not bulk_layers:
            raise TooFewLayersError(
                f'{type(self).__name__} has no bulk layers. '
                'Did you forget to .create_layers(rpars)?'
                )
        if not non_bulk_layers:
            raise TooFewLayersError(f'{type(self).__name__} has only bulk '
                                    'layers. Check LAYER_CUTS or explicitly '
                                    'define a BULK_REPEAT.')
        # Use the distance between the bottommost bulk layer and the
        # bottommost non-bulk layer, i.e., the current 'thickness of
        # bulk', plus the gap between 'bulk' and 'non-bulk' parts.
        return bulk_layers[-1].cartbotz - non_bulk_layers[-1].cartbotz          # TODO: Issue #174?

    def get_nearest_neighbours(self):
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
        distances, _ = tree.query(atom_coords_center, k=2)
        distances = distances[:, 1]  # First column is zeros
        return dict(zip(self, distances))

    def getSurfaceAtoms(self, rp):
        """Checks which atoms are 'at the surface', returns them as a set."""

        _PTL = set(el.lower() for el in PERIODIC_TABLE)
        _RADII = COVALENT_RADIUS

        atoms = copy.deepcopy(self.atlist)
        # run from top to bottom of slab
        atoms.sort(key=lambda atom: -atom.pos[2])
        covered = set()
        surfats = set()
        for atom in atoms:                                                      # TODO: make this radius into an Atom property
            if atom.el in rp.ELEMENT_MIX:
                # radius as weighted average
                totalocc = 0.0
                r = 0.0
                for chemel in rp.ELEMENT_MIX[atom.el]:
                    if chemel.lower() in _PTL:
                        if chemel in atom.site.occ:
                            r += (_RADII[chemel.capitalize()]
                                  * atom.site.occ[chemel])
                            totalocc += atom.site.occ[chemel]
                    else:
                        _LOGGER.error(
                            'Error identifying surface atoms: Could '
                            f'not identify {chemel} as a chemical element.')
                        rp.setHaltingLevel(2)
                        return []
                if totalocc == 0:
                    _LOGGER.error('Unexpected point encountered in '
                                  'generateSearchInput: GS001')
                else:
                    r /= totalocc
            else:
                # radius of atom
                chemical_element = rp.ELEMENT_RENAME.get(atom.el, atom.el)
                try:
                    r = _RADII[chemical_element.capitalize()]
                except KeyError:
                    _LOGGER.error('Error identifying surface atoms: Could not '
                                  f'identify {atom.el} as a chemical element.')
                    rp.setHaltingLevel(2)
                    return []
            r *= 1.2
            surfats.update(a for a in self
                           if a.pos[2] >= atom.pos[2] and a not in covered)
            covered.update(
                a for a in self
                if a.pos[2] < atom.pos[2] and a.is_same_xy(atom, eps=r)
                )
            if len(covered) + len(surfats) >= len(atoms):
                break   # that's all of them
        return surfats

    def has_atoms_in_multiple_c_cells(self):
        """Return whether atoms appear in multiple cells along c.

        It makes sense to call this method only before the
        first call to collapse_cartesian_coordinates() or
        collapse_fractional_coordinates(). After that, this
        method will always return True.

        Returns
        -------
        has_multiple_cells : bool
            Whether the fractional coordinates of the atoms in this
            slab are found in more than one base unit cell along c.
            Atoms very close (<1e-4) to cell boundaries are not taken
            into consideration for the test.

        Raises
        ------
        EmptySlabError
            If this method is called before this slab has any atom.
        """
        if not self.n_atoms:
            raise EmptySlabError
        cells = (at.pos[2] // 1
                 for at in self
                 if _VACUUM_EPS < at.pos[2] % 1.0 < 1 - _VACUUM_EPS)
        try:
            first_cell = next(cells)
        except StopIteration:  # All atoms close to the edges
            return False
        return any(cell != first_cell for cell in cells)

    def identify_bulk_repeat(self, eps, epsz=None):
        """Find a repeat vector for which the bulk matches the slab above.

        Notice that this may not be the shortest of such vectors: the
        vector returned is the the shortest among those that bring
        `bulkslab` to match the bottommost part of the non-bulk part
        of this slab. Should `bulkslab` be 'too thick' (e.g., the
        LAYER_CUTS select composite layers with an exaggerated number
        of sublayers), the vector returned will not be the shortest
        possible. If you are after the absolutely shortest vector, use
        `slab.bulkslab.get_minimal_c_vector(eps, epsz)`, or
        `slab.bulkslab.ensure_minimal_c_vector(rpars)`.

        Parameters
        ----------
        eps : float
            Cartesian tolerance for in-plane comparisons (Angstrom)
        epsz : float or None, optional
            Cartesian tolerance for comparisons along the direction
            perpendicular to the surface. If not given or None, take
            the same value as `eps`. This argument is used only if
            this slab or its `bulkslab` do not have yet `sublayers`.
            Default is None.

        Returns
        -------
        repeat_vec : numpy.ndarray
            Repeat vector in Cartesian coordinates, with (x, y) in the
            surface plane, and z directed from the solid towards the
            surface, i.e., opposite to the usual LEED convention.

        Raises
        ------
        MissingBulkSlabError
            If this slab has no `bulkslab`.
        NoBulkRepeatError
            If a bulk repeat vector could not be found.
        """
        if self.bulkslab is None:
            raise MissingBulkSlabError(
                f'{type(self).__name__}.identify_bulk_repeat: missing '
                'bulkslab. Call make_bulk_slab, then try again'
                )
        if epsz is None:
            epsz = eps
        if not self.sublayers:
            self.create_sublayers(epsz)
        if not self.bulkslab.sublayers:
            self.bulkslab.create_sublayers(epsz)

        n_bulk_lay = self.bulkslab.n_sublayers
        if self.n_sublayers < 2 * n_bulk_lay:
            raise NoBulkRepeatError(
                f'{type(self).__name__}.identify_bulk_repeat: failed. '
                f'Too few sublayers in slab ({self.n_sublayers}) with '
                f'respect to bulkslab ({n_bulk_lay}). Should be at least '
                f'{2*n_bulk_lay}'
                )

        # Pick the region of the slab that should be considered for
        # comparison: from the topmost bulk layer all the way down.
        # Notice that we take the topmost bulk from self rather than
        # from self.bulkslab, as that one may have a different frame
        # if it was re-centred along z when creating it.
        z_range = (
            self.sublayers[-n_bulk_lay].cartbotz,  # topmost bulk               # TODO: .cartpos[2] -- Issue #174
            self.sublayers[-1].cartbotz            # bottommost
            )

        # Construct candidate repeat vectors as all those connecting
        # the first non-bulk layer to all atoms of the bottommost layer
        first_non_bulk_layer = self.sublayers[-1-n_bulk_lay]
        if first_non_bulk_layer.element != self.sublayers[-1].element:
            # Can't detect a repeat vector if elements don't match
            raise NoBulkRepeatError(
                f'{type(self).__name__}.identify_bulk_repeat: failed. '
                'Chemical elements mismatched in first non-bulk sublayer '
                f'({first_non_bulk_layer.num}, {first_non_bulk_layer.element})'
                f' and bottommost sublayer ({self.sublayers[-1].element})'
                )

        ori = first_non_bulk_layer.cartpos                                      # TODO: need to flip with .cartpos[2]? Issue #174
        repeat_vectors = []
        for atom in self.sublayers[-1]:
            test_v = atom.cartpos - ori          # From slab to bulk
            if self.is_translation_symmetric(test_v, eps, z_periodic=False,
                                             z_range=z_range):
                repeat_vectors.append(-test_v)   # From bulk to slab
        if not repeat_vectors:
            raise NoBulkRepeatError(
                f'{type(self).__name__}.identify_bulk_repeat: failed. '
                'No repeat vectors'
                )

        # Pick the shortest repeat vector, and return it with the
        # 'conventional' z direction (positive from bulk to surface)
        shortest = min(repeat_vectors, key=np.linalg.norm)
        shortest[2] *= -1                                                       # TODO: .cartpos[2] -- Issue #174
        return shortest

    def make_bulk_slab(self, rpars, recenter=True):
        """Return a copy of this slab with only bulk layers.

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS used for making the new slab. Attributes
            used (read-only): SUPERLATTICE, superlattice_defined,
            BULK_REPEAT, SYMMETRY_EPS.
        recenter : bool, optional
            Whether the fractional coordinates of the bulk slab
            that is created should be adjusted so that the topmost
            bulk atom is at the same position. This can safely be
            left to the default as long as the bulk repeat vector
            is known and correct. Default is True.

        Returns
        -------
        bulk_slab : Slab
            Copy of self with the correct bulk unit cell. The same
            slab can then be accessed as the `bulkslab` attribute
            of this slab. Notice that no attempt is made on making
            the returned slab have a 'minimal' unit cell in either
            in- or out-of-plane directions. If you want that, you
            can later call `self.ensure_minimal_bulk_ab_cell(rpars)`
            (makes `self.bulkslab` have smallest in-plane cell) and
            `self.bulkslab.ensure_minimal_c_vector(rpars)` (also
            makes the repeat vector shortest). If you already know
            the minimal bulk c vector and/or in-plane cell, you can
            also use `self.bulkslab.apply_bulk_cell_reduction(...)`.

        Raises
        ------
        MissingLayersError
            If this method is called on a slab that has no layers.
            This normally means that `create_layers` was not called
            before.
        TooFewLayersError
            If this method is called on a slab that has only one
            layer. Normally this means that there is a problem with
            the LAYER_CUTS.
        """
        if self.n_layers < 2:
            err_txt = ('Less than two layers detected. Check POSCAR '
                       'and consider modifying LAYER_CUTS.')
            _LOGGER.error(err_txt)
            exc = TooFewLayersError(err_txt)
            if not self.layers:
                exc = MissingLayersError('No layers detected. May need to '
                                         'create_layers or to full_update')
            raise exc

        # Construct bulk slab
        bulk_slab = BulkSlab.from_slab(self)
        bulk_slab.clear_symmetry_and_ucell_history()
        bulk_slab.atlist = AtomList(bulk_slab.bulk_atoms)
        bulk_slab.layers = bulk_slab.bulk_layers
        bulk_slab.sublayers = ()

        kwargs = {
            'eps': rpars.SYMMETRY_EPS,
            'epsz': rpars.SYMMETRY_EPS.z,
            'new_c_vec': self.get_bulk_repeat(rpars),
            'new_ab_cell': np.dot(np.linalg.inv(rpars.SUPERLATTICE),
                                  self.ab_cell.T),
            'recenter': recenter,
            'z_periodic': False,  # We cannot guarantee that it will be
            }
        bulk_slab.apply_bulk_cell_reduction(**kwargs)
        a_norm, b_norm = np.linalg.norm(bulk_slab.ab_cell.T, axis=1)
        if rpars.superlattice_defined and a_norm > b_norm + 1e-4:
            _LOGGER.warning(
                'The bulk unit cell defined by SUPERLATTICE does not '
                'follow standard convention: the first vector is longer '
                'than the second. Make sure beams are labelled correctly.'
                )
            rpars.setHaltingLevel(1)
        self.bulkslab = bulk_slab
        return bulk_slab

    def make_supercell(self, transform):
        """Return a copy of the slab replicated according to `transform`.

        The 'inverse' (i.e., leading to a size reduction) of this
        operation can be obtained by calling `make_subcell` with
        the same `transform`. Atoms in the duplicated replicas are
        marked as duplicates of those in this slab.

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
            repeat_ranges = [range(r) for r in repeats]
            for atom in super_slab.atlist.copy():
                for i, j in itertools.product(*repeat_ranges):
                    # pylint: disable=compare-to-zero
                    if i == j == 0:  # Skip already existing atom
                        continue
                    # Duplicate the atom, and implicitly store it
                    # in super_slab. Then change its fractional pos
                    duplicate_atom = atom.duplicate()
                    duplicate_atom.pos[:2] += (i, j)

        # We now may have added new atoms that will have fractional
        # coordinates outside the original cell. Make sure we have
        # the correct Cartesian coordinates before transforming the
        # unit cell.
        super_slab.update_cartesian_from_fractional()
        super_slab.ab_cell[:] = super_slab.ab_cell.dot(transform.T)

        # Starting from the stored Cartesian coordinates, collapse
        # all atoms (incl. fractional coordinates) to the new cell
        super_slab.collapse_cartesian_coordinates()
        return super_slab

    def make_subcell(self, rpars, transform):
        """Return a subcell of this slab.

        This is the inverse of `make_supercell(transform)` with
        the same `transform` matrix.

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS used to determine how the copy is to be
            prepared. Attributes used: SYMMETRY_EPS.
            Both attributes are used as Cartesian tolerances for
            removal of duplicate atoms resulting from the reduction
            of size.
        transform : Sequence
            Shape (2, 2). The transformation matrix defining the
            subcell-to-supercell relationship. It should be an
            integer-valued matrix. The `transform` is inverted,
            then applied by left multiplication with the in-plane
            unit cell, when the unit cell has vectors on the rows.
            Default is None.

        Returns
        -------
        subcell_slab : Slab
            A new slab with reduced in-plane size according
            to `transform`. Duplicate atoms resulting from
            the size reduction are removed. Atoms that have
            been removed from the original slab as a result
            of the size reduction are marked as duplicates
            of those that remain in `subcell_slab`.

        Raises
        ------
        ValueError
            If `transform` is not a (2x2) matrix.
        ValueError
            If `transform` is an invalid transformation for
            this slab, that is, if slab is not a supercell
            with `transform` as its superlattice matrix.
        NonIntegerMatrixError
            If `transform` is not integer-valued.
        SlabError
            If reducing the current slab was unsuccessful.
            This is normally because either `slab.layers`
            or `rpars` were out of date.

        Notes
        -----
        This method is especially useful in case the slab
            is a supercell containing multiple copies of
            the same slab, i.e., if `get_minimal_ab_cell`
            gives a possible reduction. In that case,
            checking the symmetry on the smaller slab
            is more efficient.
        """
        subcell_slab = copy.deepcopy(self)
        subcell_slab.clear_symmetry_and_ucell_history()
        subcell_slab.update_cartesian_from_fractional()

        transform = ensure_integer_matrix(transform)
        if transform.shape != (2, 2):
            raise ValueError(f'Invalid transform.shape={transform.shape}. '
                             'Expected transform.shape=(2, 2)')
        if np.allclose(transform, np.identity(2)):
            return subcell_slab

        # Reduce dimensions in (x, y), then remove duplicates
        try:
            subcell_slab.ab_cell[:] = np.dot(self.ab_cell,
                                             np.linalg.inv(transform).T)
        except np.linalg.LinAlgError as exc:  # singular
            raise SingularMatrixError(f'transform={transform} '
                                      'is singular') from exc
        subcell_slab.collapse_cartesian_coordinates(update_origin=True)
        subcell_slab.remove_duplicate_atoms(rpars.SYMMETRY_EPS,
                                            rpars.SYMMETRY_EPS.z,
                                            other_slab=self)
        # Now make sure that transform actually gave a subcell
        supercell_atoms = subcell_slab.n_atoms * abs(np.linalg.det(transform))
        if self.n_atoms != supercell_atoms:
            raise ValueError(
                f'Slab is not a supercell with transform={transform} '
                'as its superlattice matrix'.replace('\n', ',')
                )
        return subcell_slab

    def restoreOriState(self, keepDisp=False):
        """Resets the atom positions and site vibrational amplitudes to the
        original state, and stores the deviations as offset instead."""
        for site in self.sitelist:
            if site.oriState is None:
                continue
            siteats = [at for at in self if at.site == site]
            for el in site.occ:
                for at in siteats:
                    o = site.vibamp[el] - site.oriState.vibamp[el]
                    if el in at.offset_vib:
                        o += at.offset_vib[el]
                    at.offset_vib[el] = o
                    if abs(at.offset_vib[el]) < 1e-6:
                        del at.offset_vib[el]
                site.vibamp[el] = site.oriState.vibamp[el]
                site.occ[el] = site.oriState.occ[el]
        uci = np.linalg.inv(self.ucell)
        self.update_cartesian_from_fractional()
        for at in self:
            if at.oriState is None:
                continue
            for el in at.disp_occ.keys():
                o = np.array([0., 0., 0.])
                if el in at.offset_geo:
                    o += at.offset_geo[el]
                if el in at.disp_geo_offset:
                    o += at.disp_geo_offset[el][0]
                else:
                    o += at.disp_geo_offset['all'][0]
                o[2] *= -1
                ro = np.dot(uci, o)
                nro = at.pos - at.oriState.pos + ro
                for i in range(0, len(nro)):  # minimize abs(nro)
                    while nro[i] > 0.5:
                        nro[i] -= 1
                    while nro[i] < -0.5:
                        nro[i] += 1
                no = np.dot(self.ucell, nro)
                no[2] *= -1
                at.offset_geo[el] = no
            at.pos = at.oriState.pos
            at.disp_geo_offset = {'all': [np.zeros(3)]}
        self.collapse_fractional_coordinates()
        self.update_cartesian_from_fractional()
        self.update_layer_coordinates()
        if keepDisp:
            return
        for at in self:
            at.known_deltas = []
            at.initDisp(force=True)
            at.constraints = {1: {}, 2: {}, 3: {}}
        return

    def with_extra_bulk_units(self, rpars, n_cells):
        """Return a copy with extra bulk unit cells added at the bottom.

        Parameters
        ----------
        rpars : Rparams
            The PARAMETERS  from which BULK_REPEAT is taken.
        n_cells : int
            The number of bulk unit cells to be added to the
            copy returned.

        Returns
        -------
        bulk_appended : Slab
            A copy of self with `n_cells` bulk cells added at the
            bottom.
        new_bulk_atoms : list
            The Atoms that were added to append the new bulk cells.

        Raises
        ------
        MissingLayersError
            If this method is called before any layer was identified,
            or if none of the layers is labelled as bulk.
        SlabError
            If the bulk repeat vector would give an unreasonably small
            interlayer distance. This may mean that (1) there are too
            few layers (e.g., only one thick layer), or (2) the bulk
            was not identified correctly (e.g., wrong BULK_REPEAT).
        ValueError
            If `n_cells` is not a positive number.

        Note
        ----
        Layers in the `bulk_appended` slab retuned are re-labelled
        so that there are exactly `rpars.N_BULK_LAYERS` bulk layers
        at the bottom. This means that layers of this slab that are
        labelled as bulk are turned into non-bulk layers in the
        `bulk_appended` slab returned (as are atoms in those layers).
        """
        if n_cells <= 0:
            raise ValueError('with_extra_bulk_units: n_cells must be positive')
        eps = 0.1
        bulk_appended = copy.deepcopy(self)
        bulk_layers = bulk_appended.bulk_layers
        if not bulk_layers:
            raise MissingLayersError('No bulk layers to duplicate')

        # First get the bulk repeat vector (with positive z).
        bulk_c = self.get_bulk_repeat(rpars)
        if abs(bulk_c[2]) < eps:
            raise SlabError('Bulk interlayer distance is too small: '
                            f'{bulk_c[2]}. Check LAYER_CUTS.')

        # Now take the component of bulk_c parallel to unit-cell c.
        # This will be used for expanding the unit cell.
        slab_c = bulk_appended.c_vector
        slab_c_direction = slab_c / np.linalg.norm(slab_c)
        bulk_c_par = bulk_c.dot(slab_c_direction) * slab_c_direction

        # Since we use the opposite-z convention for atoms, we also
        # need the two components with flipped z. We will move old
        # atoms up along unit-cell c, while new atoms will be shifted
        # only perpendicular to it.
        bulk_c[2] *= -1                                                         # TODO: .cartpos[2] -- Issue #174
        bulk_c_par_atoms = bulk_c.dot(slab_c_direction) * slab_c_direction
        bulk_c_perp_atoms = bulk_c - bulk_c_par_atoms

        bulk_appended.update_cartesian_from_fractional()
        new_bulk_atoms = []
        for _ in range(n_cells):
            # The trick to get it right is always adding atoms as
            # duplicates of the bottommost bulk layers, which may
            # either be the ones that this slab originally had, or
            # those that we have just added.
            # pylint: disable=protected-access
            (added_atoms,
             bulk_layers) = bulk_appended._add_one_bulk_cell(bulk_layers,
                                                             bulk_c_par,
                                                             bulk_c_par_atoms,
                                                             bulk_c_perp_atoms)
            new_bulk_atoms.extend(added_atoms)
        bulk_appended.collapse_cartesian_coordinates(update_origin=True)

        # Make sure that the BULK_REPEAT is not too short.
        # We can only cover the case in which the BULK_REPEAT
        # gives negative interlayer distances. See discussion
        # at https://github.com/viperleed/viperleed/issues/187
        if any(d < eps for d in bulk_appended.interlayer_gaps):
            raise SlabError('BULK_REPEAT is too small. Gives layers '
                            'with negative interlayer distances')

        # Make sure the number of bulk layers in the new slab is
        # still consistent with rpars.N_BULK_LAYERS. This means
        # turning all layers that were bulk before into 'non-bulk'
        for non_bulk_layer in bulk_appended.layers[:-rpars.N_BULK_LAYERS]:
            non_bulk_layer.is_bulk = False

        bulk_appended.sort_original()
        # Invalidate outdated information
        bulk_appended.sublayers = ()

        if not bulk_appended.is_bulk:
            bulk_appended.bulkslab = None
        return bulk_appended, new_bulk_atoms
