# -*- coding: utf-8 -*-
"""Module surface_slab of viperleed.tleedmlib.classes.slab.

Created on 2023-02-21, originally Jun 13 2019

@author: Florian Kraushofer (@fkraushofer)
@author: Michele Riva (@michele-riva)

Defines the SurfaceSlab class, a BaseSlab describing a 2D-periodic crystal.
This is available in tleedmlib as the Slab alias for backwards compatibility
with the former module-based structure of slab.py.
"""

from collections import Counter
import copy
import logging
from math import remainder as round_remainder

import numpy as np

from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.classes.atom_containers import AtomList
from viperleed.tleedmlib.periodic_table import PERIODIC_TABLE, COVALENT_RADIUS

from .base_slab import BaseSlab
from .bulk_slab import BulkSlab
from .slab_errors import AlreadyMinimalError, MissingBulkSlabError
from .slab_errors import NoBulkRepeatError, SlabError

try:
    import ase
except ImportError:
    _HAS_ASE = False
else:
    _HAS_ASE = True


_LOGGER = logging.getLogger('tleedm.slab')


class SurfaceSlab(BaseSlab):
    """A class representing a semi-infinite solid.

    Contains unit cell, element information and atom coordinates.
    Also has a variety of convenience functions for manipulating
    and updating the atoms. Slabs can be created from an ase.Atoms
    object using the `from_ase` class method. Another common way
    is using `viperleed.tleedmlib.files.poscar.read`.

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
        it was processed by tleedm before.
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
        self.symbaseslab = None                                                 # TODO: needed?

    @property
    def is_bulk(self):
        """Return whether this is a bulk slab."""
        return False

    @property
    def thickness(self):
        """Return the thickness of this slab along z."""
        top_atom_z = max([at.cartpos[2] for at in self])
        bottom_atom_z = min([at.cartpos[2] for at in self])
        return top_atom_z - bottom_atom_z
        # return abs(top_atom_z - bottom_atom_z)

    @property
    def vacuum_gap(self):
        """Return the z distance that does not not contain atoms."""
        return self.ucell[2, 2] - self.thickness

    @classmethod
    def from_ase(cls, ase_atoms, sort_elements=True):                           # TODO: ensure atoms are not too close to edges in c (could use poscar function, but beware of cyclic imports)
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
        TypeError
            If `ase_atoms` is not an ase.Atoms object.
        """
        if not _HAS_ASE:
            _LOGGER.warning('Slab creation from ase.Atoms not '                 # TODO: probably we should raise (ModuleNotFoundError?)
                            'supported: Module ase not found.')
            return cls()

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

    def changeBulkCell(self, rp, newcell):
        """Takes a unit cell (a,b), calculates a new SUPERLATTICE parameter
        from it and creates a new bulk slab with that cell. If a different
        SUPERLATTICE was defined by the user, outputs an error and returns."""
        # calculate new SUPERLATTICE matrix
        newSL = self.ab_cell.T.dot(np.linalg.inv(newcell))
        if not np.all(abs(newSL - newSL.round()) < 5e-2):
            _LOGGER.error(
                'Automatically detected bulk SUPERLATTICE is '
                'not integer-valued: \n'+str(newSL)+'\n'
                '# User-defined SUPERLATTICE will be used instead.')
            rp.setHaltingLevel(2)
            return
        else:
            newSL = newSL.round()
        if rp.superlattice_defined:
            _LOGGER.warning(
                'Automatically detected minimum-area bulk unit cell differs '
                'from the cell defined by the SUPERLATTICE parameter. '
                'Consider changing the SUPERLATTICE parameter. Found matrix: '
                '\n' + str(newSL.astype(int)))
            rp.setHaltingLevel(2)
            return
        else:
            rp.SUPERLATTICE = newSL
            self.bulkslab = self.makeBulkSlab(rp)

    def detectBulk(self, rp, second_cut_min_spacing=1.2):
        """
        Determine the minimal bulk repeat vector from BULK_LIKE_BELOW.

        Parameters
        ----------
        rp : Rparams
            Run parameters object. Contains BULK_LIKE_BELOW, not modified.
        second_cut_min_spacing : float, optional
            Minimal z spacing in Angstrom between sublayers of the detected
            bulk unit to make an additional cut through the bulk. If no two
            sublayers are below this threshold, the bulk will be treated as
            a single layer.

        Raises
        ------
        RuntimeError
            Raised if no repeat vector is found.

        Returns
        -------
        bulk_repeat : np.array
            The new repeat vector.
        slab_cuts : list of float
            Layer cuts needed for the bulk.

        """
        tsl = copy.deepcopy(self)   # temporary copy to mess up layers in
        rp_dummy = copy.deepcopy(rp)
        rp_dummy.LAYER_CUTS = [rp.BULK_LIKE_BELOW]
        rp_dummy.N_BULK_LAYERS = 1
        rp_dummy.BULK_LIKE_BELOW = 0.
        tsl.createLayers(rp_dummy)
        # create a pseudo-bulkslab
        tsl.bulkslab = tsl.makeBulkSlab(rp_dummy)
        bsl = tsl.bulkslab
        bsl.create_sublayers(rp.SYMMETRY_EPS_Z)
        # reduce unit cell in xy
        try:
            mincell = bsl.get_minimal_ab_cell(rp.SYMMETRY_EPS,
                                              rp.SYMMETRY_EPS_Z)
        except AlreadyMinimalError:
            pass
        else:
            tsl.changeBulkCell(rp, mincell)
            bsl = tsl.bulkslab
        bsl.update_cartesian_from_fractional()
        # detect new C vector
        newC = bsl.getMinC(rp, z_periodic=False)
        if newC is None:
            _LOGGER.error('Automatic bulk detection failed: Found no bulk '
                          'repeat vector below the specified cutoff.')
            raise RuntimeError('Failed to detect bulk repeat vector.')

        # reduce the bulk slab
        rp_dummy.BULK_REPEAT = -newC
        rp_dummy.SUPERLATTICE = np.eye(2)
        tsl.bulkslab = tsl.makeBulkSlab(rp_dummy)
        bsl = tsl.bulkslab
        bsl.create_sublayers(rp.SYMMETRY_EPS_Z)

        # calculate cut plane
        frac_bulk_thickness = abs(newC[2]) / abs(self.ucell[2, 2])
        frac_lowest_pos = min([at.pos[2] for at in self])
        frac_bulk_onset = (frac_lowest_pos + frac_bulk_thickness
                           - (rp.SYMMETRY_EPS_Z / self.ucell[2, 2]))
        slab_cuts = [(max([at.pos[2] for at in self
                          if at.pos[2] < frac_bulk_onset])
                     + min([at.pos[2] for at in self
                            if at.pos[2] > frac_bulk_onset]))
                     / 2]

        # now look for potential second cut
        if bsl.n_sublayers > 1:
            maxdist = abs(bsl.sublayers[1].cartbotz
                          - bsl.sublayers[0].cartbotz)
            cutlayer = 0
            for i in range(1, bsl.n_sublayers - 1):
                d = abs(bsl.sublayers[i+1].cartbotz
                        - bsl.sublayers[i].cartbotz)
                if d > maxdist:
                    maxdist = d
                    cutlayer = i
            if maxdist >= second_cut_min_spacing:
                bulkcut_frac_from_lowest = (
                    bsl.sublayers[cutlayer].pos[2]
                    - (maxdist / (2 * abs(bsl.ucell[2, 2])))
                    - min([at.pos[2] for at in bsl]))
                slab_cuts.append(frac_lowest_pos
                                 + (bulkcut_frac_from_lowest
                                    * bsl.ucell[2, 2] / self.ucell[2, 2]))
        return (-newC, slab_cuts)

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
        """
        if isinstance(rpars.BULK_REPEAT, np.ndarray):
            bulkc = rpars.BULK_REPEAT.copy()
            if bulkc[2] < 0:
                bulkc *= -1
            return bulkc

        cvec = self.ucell.T[2]
        if rpars.BULK_REPEAT is None:
            # Use the distance between the bottommost bulk layer
            # and the bottommost non-bulk layer, i.e., the current
            # 'thickness of bulk', plus the gap between 'bulk' and
            # 'non-bulk' parts.
            blayers = self.bulk_layers
            zdiff = (blayers[-1].cartbotz
                     - self.layers[blayers[0].num - 1].cartbotz)
        else:  # rpars.BULK_REPEAT is float
            zdiff = rpars.BULK_REPEAT
        if only_z_distance:
            return zdiff
        return cvec * zdiff / cvec[2]

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
            r *= 1.2    # TODO: !!! test if this is enough
            surfats.update(a for a in self
                           if (a.pos[2] >= atom.pos[2]
                               and a not in covered))
            covered.update(
                a for a in self
                if a.pos[2] < atom.pos[2] and a.is_same_xy(atom, eps=r)
                )
            if len(covered) + len(surfats) >= len(atoms):
                break   # that's all of them
        return surfats

    def identify_bulk_repeat(self, eps, epsz=None):
        """Find a repeat vector for which the bulk matches the slab above.

        Notice that this may not be the shortest of such vectors: the
        vector returned is the the shortest among those that bring
        `bulkslab` to match the bottommost part of the non-bulk part
        of this slab. Should `bulkslab` be 'too thick' (e.g., the
        LAYER_CUTS select composite layers with an exaggerated number
        of sublayers), the vector returned will not be the shortest
        possible.

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
                'bulkslab. Call makeBulkSlab, then try again'
                )
        if epsz is None:
            epsz = eps
        if not self.sublayers:
            self.create_sublayers(epsz)
        if not self.bulkslab.sublayers:
            self.bulkslab.create_sublayers(epsz)

        n_bulk_lay = self.bulkslab.n_sublayers
        if self.n_sublayers < 2*n_bulk_lay:
            raise NoBulkRepeatError(
                f'{type(self).__name__}.identify_bulk_repeat: failed. '
                f'Too few sublayers in slab ({self.n_sublayers}) with '
                f'respect to bulkslab ({n_bulk_lay}). Should be at least '
                f'{2*n_bulk_lay}'
                )

        # Pick the region of the slab that should be considered for
        # comparison: from the topmost bulk layer all the way down
        z_range = (self.bulkslab.sublayers[0].cartbotz,                         # TODO: @fkraushofer used to be self.sublayers[-n_bulk_lay], and both used .cartpos[2] instead of .cartbotz. I think that using self.bulkslab.sublayers[0] is equivalent, and that using .cartbotz is a bit better as it may loose atoms at the bottom. Did I misinterpret something?
                   self.sublayers[-1].cartbotz)

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

        ori = first_non_bulk_layer.cartpos                                      # TODO: need to flip with .cartpos[2]?
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
        shortest[2] *= -1                                                       # TODO: .cartpos[2]
        return shortest

    def makeBulkSlab(self, rp):
        """Copies self to create a bulk slab, in which everything apart from the
        bulk layers is deleted. Returns that bulk slab."""
        if self.n_layers < 2:
            err_txt = 'Less than two layers detected. Check POSCAR and consider modifying LAYER_CUTS.'
            _LOGGER.error(err_txt)
            raise RuntimeError(err_txt)

        # construct bulk slab
        bsl = BulkSlab.from_slab(self)
        bsl.clear_symmetry_and_ucell_history()
        bsl.atlist = AtomList(at for at in bsl if at.is_bulk)
        bsl.layers = bsl.bulk_layers
        bsl.update_cartesian_from_fractional()
        al = bsl.atlist[:]     # temporary copy
        al.sort(key=lambda atom: atom.pos[2])
        topat = al[-1]
        if type(rp.BULK_REPEAT) == np.ndarray:
            bulkc = np.copy(rp.BULK_REPEAT)
            if bulkc[2] < 0:
                bulkc = -bulkc
        else:
            cvec = bsl.ucell[:, 2]
            zdiff = 0
            if rp.BULK_REPEAT is None:
                # assume that interlayer vector from bottom non-bulk to top
                #  bulk layer is the same as between bulk units
                zdiff = (bsl.layers[-1].cartbotz
                         - self.layers[bsl.layers[0].num-1].cartbotz)
            elif isinstance(rp.BULK_REPEAT, (float, np.floating)):
                zdiff = rp.BULK_REPEAT
            if  abs(zdiff) < 1e-5:
                err_txt = 'Unable to detect bulk interlayer vector. Check POSCAR and consider explicitly setting BULK_REPEAT.'
                _LOGGER.error(err_txt)
                raise RuntimeError(err_txt)
            bulkc = cvec * zdiff / cvec[2]
        bsl.ucell[:, 2] = bulkc
        # reduce dimensions in xy
        superlattice = np.identity(3, dtype=float)
        superlattice[:2, :2] = rp.SUPERLATTICE.T
        bsl.ucell = np.dot(bsl.ucell, np.linalg.inv(superlattice))
        a_vec, b_vec = bsl.ab_cell.T
        if (rp.superlattice_defined
                and np.linalg.norm(a_vec) > np.linalg.norm(b_vec) + 1e-4):
            _LOGGER.warning(
                'The bulk unit cell defined by SUPERLATTICE does '
                'not follow standard convention: the first vector is larger '
                'than the second. Make sure beams are labelled correctly.')
            rp.setHaltingLevel(1)
        bsl.collapse_cartesian_coordinates(update_origin=True)
        # if self.ucell_mod is not empty, don't drag that into the bulk slab.
        bsl.ucell_mod = []
        # position in c is now random; make topmost bulk atom top
        topc = topat.pos[2]
        for at in bsl:
            at.pos[2] = (at.pos[2] - topc + 0.9999) % 1.
        # now center around cell midpoint
        cl = [at.pos[2] for at in bsl]
        midpos = (max(cl)+min(cl))/2
        for at in bsl:
            at.pos[2] = (at.pos[2] - midpos + 0.5) % 1.
        bsl.update_cartesian_from_fractional(update_origin=True)
        bsl.update_element_count()   # update the number of atoms per element
        # remove duplicates
        bsl.create_sublayers(rp.SYMMETRY_EPS_Z)
        newatlist = AtomList()
        for subl in bsl.sublayers:
            i = 0
            while i < subl.n_atoms:
                j = i+1
                while j < subl.n_atoms:
                    if subl.atlist[i].is_same_xy(subl.atlist[j],
                                                 eps=rp.SYMMETRY_EPS):
                        subl.atlist.pop(j)
                    else:
                        j += 1
                i += 1
            newatlist.extend(subl)
        bsl.atlist = newatlist
        bsl.update_element_count()   # update number of atoms per element again
        # update the layers. Don't use Slab.createLayers here to keep it
        #   consistent with the slab layers
        for i, layer in enumerate(bsl.layers):
            layer.slab = bsl
            layer.update_position()
            layer.num = i
            layer.atlist = [at for at in layer if at in bsl]
        return bsl

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

        It is important to realize that the `bulk_appended` slab
        returned will have **more** than `rpars.N_BULK_LAYERS` bulk
        layers: There will be n_cells * len(self.bulk_layers).

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
        SlabError
            If this method is called before any layer was identified,
            or if none of the layers is labelled as bulk.
        SlabError
            If the bulk repeat vector would give an unreasonably small
            interlayer distance. This may mean that (1) there are too
            few layers (e.g., only one thick layer), or (2) the bulk
            was not identified correctly (e.g., wrong BULK_REPEAT).
        """
        bulk_appended = copy.deepcopy(self)
        bulk_layers = bulk_appended.bulk_layers
        if not bulk_layers:
            raise SlabError('No bulk layers to duplicate')

        # First get the bulk repeat vector (with positive z). Take into
        # account that there already may be multiple bulk layers at the
        # bottom. In that case the thickness of the bulk is larger than
        # the bulk repeat, and it must be a multiple of its z component
        bulk_c = self.get_bulk_repeat(rpars)
        if abs(bulk_c[2]) < 0.1:
            raise SlabError('Bulk interlayer distance is too small: '
                            f'{bulk_c[2]}. Check LAYER_CUTS.')
        bulk_thickness = abs(bulk_layers[0].cartori[2]
                             - bulk_layers[-1].cartbotz)
        if bulk_thickness > bulk_c[2]:
            assert abs(round_remainder(bulk_thickness, bulk_c[2])) < 1e-3
            bulk_c *= round(bulk_thickness/bulk_c[2]) + 1

        # Now take the component of bulk_c parallel to unit-cell c.
        # This will be used for expanding the unit cell.
        slab_c = bulk_appended.ucell.T[2]
        slab_c_direction = slab_c / np.linalg.norm(slab_c)
        bulk_c_par = bulk_c.dot(slab_c_direction) * slab_c_direction

        # Since we use the opposite-z convention for atoms, we also
        # need the two components with flipped z: We will move old
        # atoms up along unit-cell c, while new atoms will be shifted
        # only in the ab plane.
        bulk_c[2] *= -1                                                          # TODO: .cartpos[2]
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

        bulk_appended.sort_original()
        # Invalidate outdated information
        bulk_appended.sublayers = []

        if not bulk_appended.is_bulk:
            bulk_appended.bulkslab = None
        return bulk_appended, new_bulk_atoms
