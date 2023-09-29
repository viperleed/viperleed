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

import numpy as np

from viperleed.tleedmlib.classes.atom import Atom
from viperleed.tleedmlib.periodic_table import PERIODIC_TABLE, COVALENT_RADIUS

from .base_slab import BaseSlab
from .bulk_slab import BulkSlab


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
    chemelem : list of str
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
        self.symbaseslab = None

    @property
    def is_bulk(self):
        """Return whether this is a bulk slab."""
        return False

    @property
    def thickness(self):
        """Return the thickness of this slab along z."""
        top_atom_z = max([at.cartpos[2] for at in self.atlist])
        bottom_atom_z = min([at.cartpos[2] for at in self.atlist])
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
            n_atoms = len(self.atlist)
            self.atlist.append(Atom(elem, pos, n_atoms + 1, self))
        self.getCartesianCoordinates()
        return self

    def changeBulkCell(self, rp, newcell):
        """Takes a unit cell (a,b), calculates a new SUPERLATTICE parameter
        from it and creates a new bulk slab with that cell. If a different
        SUPERLATTICE was defined by the user, outputs an error and returns."""
        # calculate new SUPERLATTICE matrix
        abst = np.transpose(self.ucell[:2, :2])
        newSL = np.dot(abst, np.linalg.inv(newcell))
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
        bsl.createSublayers(rp.SYMMETRY_EPS_Z)
        # reduce unit cell in xy
        changecell, mincell = bsl.getMinUnitCell(rp)
        if changecell:
            tsl.changeBulkCell(rp, mincell)
            bsl = tsl.bulkslab
        bsl.getCartesianCoordinates()
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
        bsl.createSublayers(rp.SYMMETRY_EPS_Z)

        # calculate cut plane
        frac_bulk_thickness = abs(newC[2]) / abs(self.ucell[2, 2])
        frac_lowest_pos = min([at.pos[2] for at in self.atlist])
        frac_bulk_onset = (frac_lowest_pos + frac_bulk_thickness
                           - (rp.SYMMETRY_EPS_Z / self.ucell[2, 2]))
        slab_cuts = [(max([at.pos[2] for at in self.atlist
                          if at.pos[2] < frac_bulk_onset])
                     + min([at.pos[2] for at in self.atlist
                            if at.pos[2] > frac_bulk_onset]))
                     / 2]

        # now look for potential second cut
        if len(bsl.sublayers) > 1:
            maxdist = abs(bsl.sublayers[1].cartbotz
                          - bsl.sublayers[0].cartbotz)
            cutlayer = 0
            for i in range(1, len(bsl.sublayers) - 1):
                d = abs(bsl.sublayers[i+1].cartbotz
                        - bsl.sublayers[i].cartbotz)
                if d > maxdist:
                    maxdist = d
                    cutlayer = i
            if maxdist >= second_cut_min_spacing:
                bulkcut_frac_from_lowest = (
                    bsl.sublayers[cutlayer].atlist[0].pos[2]
                    - (maxdist / (2 * abs(bsl.ucell[2, 2])))
                    - min([at.pos[2] for at in bsl.atlist]))
                slab_cuts.append(frac_lowest_pos
                                 + (bulkcut_frac_from_lowest
                                    * bsl.ucell[2, 2] / self.ucell[2, 2]))
        return (-newC, slab_cuts)

    def getBulkRepeat(self, rp):
        """Based on a pre-existing definition of the bulk, tries to identify
        a repeat vector for which the bulk matches the slab above. Returns that
        vector in cartesian coordinates, or None if no match is found."""
        eps = rp.SYMMETRY_EPS
        if len(self.sublayers) == 0:
            self.createSublayers(rp.SYMMETRY_EPS_Z)
        if self.bulkslab is None:
            self.makeBulkSlab(rp)
        if len(self.bulkslab.sublayers) == 0:
            self.bulkslab.createSublayers(rp.SYMMETRY_EPS_Z)
        nsub = len(self.bulkslab.sublayers)
        if len(self.sublayers) < 2*nsub:
            return None
        # nonbulk_subl = self.sublayers[:-nsub]
        z_range = (self.sublayers[-nsub].atlist[0].cartpos[2],
                   self.sublayers[-1].atlist[0].cartpos[2])
        baseLayer = self.sublayers[-1-nsub]
        ori = baseLayer.atlist[0].cartpos  # compare displacements from here
        repeat_vectors = []
        for at in self.sublayers[-1].atlist:
            v = at.cartpos - ori
            if self.isTranslationSymmetric(v, eps, z_periodic=False,
                                           z_range=z_range):
                repeat_vectors.append(-v)
        if len(repeat_vectors) == 0:
            return None
        # pick the shortest repeat vector
        cv = min(repeat_vectors, key=lambda x: np.linalg.norm(x))
        cv[2] = -cv[2]  # leed coordinates to standard
        return cv

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
            surfats.update(a for a in self.atlist
                           if (a.pos[2] >= atom.pos[2]
                               and a not in covered))
            covered.update(a for a in self.atlist
                           if (a.pos[2] < atom.pos[2]
                               and a.isSameXY(atom.cartpos[:2], eps=r)))
            if len(covered) + len(surfats) >= len(atoms):
                break   # that's all of them
        return surfats

    def makeBulkSlab(self, rp):
        """Copies self to create a bulk slab, in which everything apart from the
        bulk layers is deleted. Returns that bulk slab."""
        if  len(self.layers) < 2:
            err_txt = 'Less than two layers detected. Check POSCAR and consider modifying LAYER_CUTS.'
            _LOGGER.error(err_txt)
            raise RuntimeError(err_txt)

        # construct bulk slab
        bsl = BulkSlab.from_slab(self)
        bsl.resetSymmetry()
        bsl.atlist = [at for at in bsl.atlist if at.layer.isBulk]
        bsl.layers = [lay for lay in bsl.layers if lay.isBulk]
        bsl.getCartesianCoordinates()
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
        if (rp.superlattice_defined and np.linalg.norm(bsl.ucell[:2, 0]) >
                np.linalg.norm(bsl.ucell[:2, 1]) + 1e-4):
            _LOGGER.warning(
                'The bulk unit cell defined by SUPERLATTICE does '
                'not follow standard convention: the first vector is larger '
                'than the second. Make sure beams are labelled correctly.')
            rp.setHaltingLevel(1)
        bsl.collapseCartesianCoordinates(updateOrigin=True)
        # if self.ucell_mod is not empty, don't drag that into the bulk slab.
        bsl.ucell_mod = []
        # position in c is now random; make topmost bulk atom top
        topc = topat.pos[2]
        for at in bsl.atlist:
            at.pos[2] = (at.pos[2] - topc + 0.9999) % 1.
        # now center around cell midpoint
        cl = [at.pos[2] for at in bsl.atlist]
        midpos = (max(cl)+min(cl))/2
        for at in bsl.atlist:
            at.pos[2] = (at.pos[2] - midpos + 0.5) % 1.
        bsl.getCartesianCoordinates(updateOrigin=True)
        bsl.updateElementCount()   # update the number of atoms per element
        # remove duplicates
        bsl.createSublayers(rp.SYMMETRY_EPS_Z)
        newatlist = []
        for subl in bsl.sublayers:
            i = 0
            while i < len(subl.atlist):
                j = i+1
                while j < len(subl.atlist):
                    if subl.atlist[i].isSameXY(subl.atlist[j].cartpos[:2],
                                               eps=rp.SYMMETRY_EPS):
                        subl.atlist.pop(j)
                    else:
                        j += 1
                i += 1
            newatlist.extend(subl.atlist)
        bsl.atlist = newatlist
        bsl.updateElementCount()   # update number of atoms per element again
        # update the layers. Don't use Slab.createLayers here to keep it
        #   consistent with the slab layers
        for i, layer in enumerate(bsl.layers):
            layer.slab = bsl
            layer.getLayerPos()
            layer.num = i
            layer.atlist = [at for at in layer.atlist if at in bsl.atlist]
        return bsl

    def restoreOriState(self, keepDisp=False):
        """Resets the atom positions and site vibrational amplitudes to the
        original state, and stores the deviations as offset instead."""
        for site in self.sitelist:
            if site.oriState is None:
                continue
            siteats = [at for at in self.atlist if at.site == site]
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
        self.getCartesianCoordinates()
        for at in self.atlist:
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
        self.collapseFractionalCoordinates()
        self.getCartesianCoordinates()
        self.updateLayerCoordinates()
        if keepDisp:
            return
        for at in self.atlist:
            at.deltasGenerated = []
            at.initDisp(force=True)
            at.constraints = {1: {}, 2: {}, 3: {}}
        return

    def updateAtomNumbers(self):
        """Updates atom oriN - should not happen normally, but necessary if
        atoms get deleted."""
        for (i, at) in enumerate(self.atlist):
            at.oriN = i+1

