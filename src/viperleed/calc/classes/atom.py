"""Module atom of viperleed.calc.classes.

Class storing position and other properties of individual atoms (to be
used with Slab, Layer, etc).
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2019-06-13'
__license__ = 'GPLv3+'

import copy
from dataclasses import field
import logging
from typing import Dict


import numpy as np

from viperleed.calc.lib.coordinates import add_edges_and_corners
from viperleed.calc.lib.dataclass_utils import frozen

_LOGGER = logging.getLogger(__name__)


class AtomError(Exception):
    """Base exception for Atom objects."""


class Atom:                                                                     # TODO: Issue #174 -- change description of cartpos when flipping .cartpos[2]
    """Class storing information about a single atom.

    Attributes
    ----------
    el : str
        Element type
    pos : numpy.ndarray
        Atom position as fractional coordinate
    num : int
        A progressive number identifying this atom. Normally, it would
        be the atom number in the POSCAR (i.e., the VESTA progressive,
        non-element-specific number).
    slab : Slab
        The Slab that the atom belongs to.
    layer : Layer or None
        Layer object that the atom belongs to (if any).
    site : Sitetype or None
        Site type of the atom, containing vibration amplitude and
        occupation. Supplied by the SITE_DEF parameter, assigned by
        the initSites method of Slab.
    cartpos : numpy.ndarray
        Position in Cartesian coordinates, with the highest atom
        at z = 0, positive z pointing into the surface
    linklist : list of Atom
        Defines to which other atoms this Atom is linked. Notice that
        if this Atom is linked to another one, the linklist of both
        are the same object, i.e., there is only one linklist shared
        by all atoms in the same group. This is used together with
        symrefm to determine symmetry-conserving displacements.
    displist : list of Atom
        Like linklist, but keeps track of which symmetry was active
        when displacement was defined.
    freedir : int or numpy.ndarray
        Defines whether the atom can be moved or is locked by symmetry.
        0: no movements, 1: completely free,
        np.array([0|1, 0|1|-1]): parallel movement to a, b, or diagonal
    symrefm : numpy.ndarray
        Defines how a translation of the atom at self.linklist[0]
        should be transformed to affect this atom.
    disp_vib, disp_geo, disp_occ : dict
        Keys are elements, values are lists of vibration/geometry/
        occupation offsets
    disp_geo_offset : dict
        Keys are elements; values are also formatted as lists
        for convenience, but should be one element long.
    disp_center_index : dict of dict
        Which index in the displacement range corresponds to
        'no change'
    dispInitialized : bool
        disp_* variables get initialized after readVIBROCC by Atom.initDisp
    known_deltas : list of str
        Filenames of delta files generated or found for this atom
    offset_geo, offset_vib, offset_occ : dict
        Offsets from self.cartpos, self.site.vib, self.site.occ
        per element
    constraints : dict
        Parameter constraints for restrict.f per element, where
        1, 2, 3 = geo, vib, occ. Can be integer-valued index in
        disp_* range or a tuple (atom, element)
    oriState : Atom
        Deep copy of self before a search is applied
    duplicate_of : Atom
        If this atom is identical to another one by translational
        symmetry (through SYMMETRY_CELL_TRANSFORM or domain supercell
        creation), this points to the atom in the base cell.
    """

    def __init__(self, el, pos, num, slab):
        """Initialize instance."""
        self.el = el
        self.pos = np.asarray(pos, dtype=float)
        self.cartpos = None                                                     # TODO: Issue #174 -- calculate right away
        self.num = num
        self.slab = slab
        self.layer = None
        self.site = None

        self.duplicate_of = None
        self.known_deltas = []
        self.oriState = None

        self.linklist = []
        self.freedir = 1
        self.symrefm = np.identity(2)

        self.displist = []
        self.disp_vib = {'all': [0.]}
        self.disp_geo = {'all': [np.zeros(3)]}
        self.disp_occ = {el: [1.]}
        self.disp_labels = {'geo': '',
                            'vib': '',
                            'occ': ''}
        self.disp_lin_steps = {'geo': [],
                               'vib': [],
                               'occ': []}
        self.disp_geo_offset = {'all': [np.zeros(3)]}
        self.disp_center_index = {'vib': {'all': 0},
                                  'geo': {'all': 0},
                                  'occ': {el: 0}}
        self.disp_ranges = DisplacementRanges(_atom=self)
        self.dispInitialized = False
        self.offset_geo = {}
        self.offset_vib = {}
        self.offset_occ = {}
        self.constraints = {1: {}, 2: {}, 3: {}}
        self.disp_ranges.update()

    def __repr__(self):
        """Return a representation string of this Atom."""
        return (f'{type(self).__name__}(el={self.el}, pos={self.pos}, '
                f'num={self.num}, slab={self.slab})')

    def __str__(self):
        """Return a string version of this Atom."""
        return f'{type(self).__name__}({self.num} {self.el})'

    def __format__(self, fmt_spec):
        """Return a formatted version of self."""
        return format(str(self), fmt_spec)

    @property
    def is_bulk(self):
        """Return whether this Atom is in a bulk layer."""
        try:
            return self.layer.is_bulk
        except AttributeError:
            return False

    def assignConstraint(self, mode, targetel='', value=None, linkAtEl=None,
                         index=None):
        """
        Assign a displacement constraint to this atom. Can be assigned for all
        elements or only one. Constraint is either a fixed value, or another
        (Atom, element) pair to link to.

        Parameters
        ----------
        mode : integer
            Which parameter to constrain. 1: geo, 2: vib, 3: occ
        targetel : string, optional
            If undefined, constrain for all elements. Otherwise,
            only constrain for the given element.
        value : float, optional
            The value to which the given parameter should be
            constrained. If the value is not in the disprange,
            assignment will be skipped. Leave at default (None)
            if 'linkAtEl' or 'index' argument is passed instead.
        linkAtEl : tuple, optional
            Format (Atom, str). The atom and element to which the
            parameter should be linked. If that atom's displist has
            a different length, assignment will be skipped. The
            element may be an empty string. In that case linking
            is considered to apply to all elements. Leave at default
            if 'value' or 'index' argument is passed instead. Default
            is None.
        index : int, optional
            The index to which the given parameter should be
            constrained. Leave at default if 'linkAtEl' or
            'value' argument is passed instead. Indices are
            1-based.

        Raises
        ------
        TypeError
            If no constrain type is given among
            `value`, `linkAtEl`, and `index`
        ValueError
            If more than one constraint type is given among
            `value`, `linkAtEl`, and `index`
        ValueError
            If `mode` is not one of the acceptable displacement modes
        """
        eps = 1e-5
        pars = len([p for p in [value, linkAtEl, index] if p is not None])
        if pars > 1:
            raise ValueError(
                f'{type(self).__name__}.assignConstraint: Can only constrain '
                'to either a fixed value or to another atom, not both'
                )
        if not pars:
            raise TypeError(
                f'{type(self).__name__}.assignConstraint: Exactly one '
                'constraint needed among "index", "value", or "linkAtEl"'
                )

        if mode == 1:
            td = self.disp_geo
        elif mode == 2:
            td = self.disp_vib
        elif mode == 3:
            td = self.disp_occ
        else:  # offset is not allowed here
            raise ValueError(f'{type(self).__name__}.assignConstraint: '
                             f'Unknown key {mode} for mode ({self})')

        if index is not None or value is not None:
            if targetel == '':
                els = list(td.keys())
            else:
                if targetel in td:
                    els = [targetel]
                elif 'all' in td:
                    els = ['all']
                else:
                    _LOGGER.warning(
                        f'Cannot assign constraint for {self}: Element '
                        f'{targetel} does not have displacements assigned.'
                        )
                    return
            for el in els:
                if index:
                    if index > len(td[el]) or index < 1:
                        _LOGGER.warning(
                            f'Cannot assign constraint for {self}, '
                            f'element {el}: index {index} is out of bounds'
                            )
                    else:
                        self.constraints[mode][el] = index - 1
                    continue
                # else value
                if mode == 1:
                    dirvec = td[el][-1] - td[el][0]  # dir of disp range
                    dirvec = dirvec / np.linalg.norm(dirvec)
                    match = [(np.linalg.norm(v - value*dirvec) < eps)
                             for v in td[el]]
                else:
                    match = [(abs(v - value) < eps) for v in td[el]]
                if any(match):
                    self.constraints[mode][el] = match.index(True)
                else:
                    _LOGGER.warning(
                        f'Cannot assign constraint for {self}, element {el}: '
                        f'value {value} is not in displacement list'
                        )
        else:  # linkAtEl
            (at, el2) = linkAtEl
            if mode == 1:
                td2 = at.disp_geo
            elif mode == 2:
                td2 = at.disp_vib
            elif mode == 3:
                td2 = at.disp_occ
            if el2 in td2:
                listlen = len(td2[el2])
            else:
                listlen = len(list(td2.values())[-1])
            if targetel == '':
                els = list(td.keys())
            else:
                els = [targetel]
            for el in els:
                if len(td[el]) != listlen:
                    _LOGGER.warning(f'Cannot constrain {self} to {at}: '
                                    'displacement list lengths differ.')
                    return
            for el in els:
                self.constraints[mode][el] = linkAtEl

    def assignDisp(self, mode, disprange, targetel='', primary=True,
                   displist=[], disp_label='', disp_lin_steps = []):
        """
        Assigns a list of displacements to this atom, for all or a given
        element.

        Parameters
        ----------
        mode : int
            What to displace. 1: geo, 2: vib, 3: occ, 4: geo offset
        disprange :
            The list of displacements. For geometric offsets, pass
            a list with only one element.
        targetel : str, optional
            If given, assignment is made only for that element,
            otherwise for all. Default is an empty string.
        primary : bool, optional
            Defines whether the assigned displacement should be passed
            along to linked atoms. If True, assignDisp is called for
            these atoms, with primary=False. FOR INTERNAL USE ONLY.
            Default is True.
        displist : list, optional
            Elements are Atom objects. Passed in secondary assignment
            to later link parameters (the 'linklist' defines how the
            'displist' is defined, but can change via the SYM_DELTA
            parameter).

        Raises
        ------
        ValueError
            If `mode` is not an acceptable displacement mode.
        """
        if targetel.lower() == 'vac':
            # don't write stuff for vacancies - they don't get displaced, and
            #   occupation can be determined implicitly
            return
        eps = 1e-5  # tolerance for comparing displacements
        dr = copy.copy(disprange)
        # to make sure multiple atoms do not get the same list object
        if mode == 1:
            td = self.disp_geo
            self.disp_labels['geo'] = disp_label
            self.disp_lin_steps['geo'] = disp_lin_steps
        elif mode == 2:
            td = self.disp_vib
            self.disp_labels['vib'] = 'N/A'  # direction not applicable for vib
            self.disp_lin_steps['vib'] = disp_lin_steps
        elif mode == 3:
            td = self.disp_occ
            self.disp_labels['occ'] = 'N/A'  # direction not applicabel for occ
            self.disp_lin_steps['occ'] = disp_lin_steps
        elif mode == 4:
            td = self.disp_geo_offset
        else:
            raise ValueError(f'{type(self).__name__}.assignDisp: '
                             f'Unknown key {mode} for mode ({self})')
        if targetel == '':
            els = list(td.keys())
        else:
            els = [targetel]
        for el in els:
            if mode == 4:  # special case: offset
                if el not in td:
                    td[el] = dr
                else:
                    match = True
                    if (abs(dr[0][2]) > eps and
                            abs(td[el][0][2] - dr[0][2]) > eps):
                        if abs(td[el][0][2]) < eps:
                            td[el][0][2] = dr[0][2]
                        else:
                            match = False
                    if (np.linalg.norm(dr[0][:2]) > eps and
                            np.linalg.norm(td[el][0][:2] - dr[0][:2]) > eps):
                        if all([(abs(f) < eps) for f in td[el][0][:2]]):
                            td[el][0][:2] = dr[0][:2]
                        else:
                            match = False
                if not match:
                    _LOGGER.warning(
                        'Atom.assignDisp: Trying to assign offset, '
                        f'but {self} already has an offset defined.'
                        ' Skipping second assignment.'
                        )
                continue
            if (el not in td
                    or (len(td[el]) == 1 and np.linalg.norm(td[el]) < 1e-5)
                    or (mode == 3 and len(td[el]) == 1)):
                # did not assign for this element yet -> OK, store
                td[el] = dr
                # also store center
                if mode != 3:
                    n = [np.linalg.norm(v) for v in dr]
                else:
                    n = [abs(v - self.site.occ[el]) for v in dr]
                smode = {1: 'geo', 2: 'vib', 3: 'occ'}
                self.disp_center_index[smode[mode]][el] = n.index(min(n))
                continue
            # is warning required: every value in td[el] also in disprange?
            match = True
            if el in td:
                if len(td[el]) != len(dr):
                    match = False
                else:
                    for v in td[el]:
                        # check whether all values from td[el] are in disprange
                        found = False
                        for v2 in dr:
                            if np.linalg.norm(v-v2) < eps:
                                found = True
                                break
                        if not found:
                            match = False
                            break
                if not match and len(td[el]) > 1:
                    _LOGGER.warning(
                        'Atom.assignDisp: Trying to assign displacement list, '
                        f'but {self} already has displacements assigned. '
                        'Skipping second assignment.'
                        )
                    return    # also don't assign to other atoms
                # in case of linking, base assignments on current dr
                dr = td[el][:]
        # assign to atoms in linklist:
        if primary:
            if not self.displist:
                self.displist = [self]
                self.slab.displists.append(self.displist)
            for at in [at for at in self.linklist if at != self]:
                if mode == 1 or mode == 4:
                    tm = np.identity(3)
                    tm[:2, :2] = np.dot(at.symrefm,
                                        np.linalg.inv(self.symrefm))
                    newdr = [np.dot(tm, v) for v in dr]
                else:
                    newdr = dr[:]
                at.assignDisp(mode, newdr, targetel, primary=False,
                              displist=self.displist,
                              disp_label=disp_label,
                              disp_lin_steps=disp_lin_steps)
            return
        if self.displist and displist != self.displist:
            _LOGGER.warning(
                f'{self} is being linked to different groups in the '
                'DISPLACEMENTS file. This will interfere with correct '
                'parameter linking in the search! Check SYM_DELTA settings.'
                )
        if self not in displist:
            displist.append(self)
        self.displist = displist

    def clearOffset(self, mode, targetel='', primary=True, displist=[]):
        """Revert an offset for self and all its symmetry-equivalent atoms.

        The offset restored is the one saved in .oriState
        (typically from POSCAR or VIBROCC).

        Parameters
        ----------
        mode : int
            Which offset to restore. 1: geo, 2: vib, 3: occ
        targetel : str, optional
            If passed, assignment is made only for that element,
            otherwise for all.
        primary : bool, optional
            Defines whether assignment should be passed along
            to linked atoms. This will call assignDisp for these
            atoms, with primary=False.
        displist : list, optional
            Elements are Atom objects. Passed in secondary assignment
            to later link parameters (the 'linklist' defines how the
            'displist' is defined, but can change via the SYM_DELTA
            parameter).

        Raises
        ------
        ValueError
            If `mode` is not one of the acceptable displacement modes.
        """
        if self.oriState is None or targetel.lower() == 'vac':
            return
        if mode == 1:
            td = self.offset_geo
            od = self.oriState.offset_geo
        elif mode == 2:
            td = self.offset_vib
            od = self.oriState.offset_vib
        elif mode == 3:
            td = self.offset_occ
            od = self.oriState.offset_occ
        else:
            raise ValueError(f'{type(self).__name__}.clearOffset: '
                             f'Unknown key {mode} for mode ({self})')
        if targetel == '':
            els = list(td.keys())
        else:
            els = [targetel]
        for el in els:
            if el not in od:
                del td[el]
            else:
                td[el] = od[el]
        # assign to atoms in linklist:
        if primary:
            if not self.displist:
                self.displist = [self]
                self.slab.displists.append(self.displist)
            for at in [at for at in self.linklist if at != self]:
                at.clearOffset(mode, targetel, primary=False,
                               displist=self.displist)
            return
        if self.displist and displist != self.displist:
            _LOGGER.warning(
                f'{self} is being linked to different groups in the '
                'DISPLACEMENTS file. This will interfere with correct '
                'parameter linking in the search! Check SYM_DELTA settings.'
                )
        if self not in displist:
            displist.append(self)
        self.displist = displist

    def copyOriState(self, other):
        """Deepcopy positions and offsets from another atom into oriState."""
        self.storeOriState()
        self.oriState.pos = other.pos.copy()
        self.oriState.cartpos = other.cartpos.copy()
        self.oriState.offset_geo = copy.deepcopy(other.offset_geo)
        self.oriState.offset_vib = copy.deepcopy(other.offset_vib)
        self.oriState.offset_occ = copy.deepcopy(other.offset_occ)

    def duplicate(self, add_to_atlists=True, num=None):
        """Return a somewhat lightweight copy of this Atom.

        Attributes position and elements are deep-copied, all others
        are instead references to those of this Atom. This includes
        in particular, site, displacements, slab, and layer. The new
        Atom can also be automatically added to the existing slab and
        layer.

        Parameters
        ----------
        add_to_atlists : bool, optional
            Whether the duplicate atom should be added to atom lists
            in the existing Slab and Layer objects. Default is True.
        num : int or None, optional
            Progressive number to give to the new atom returned.
            If not given or None, it is taken from the slab of
            this Atom. Default is None.

        Returns
        -------
        newat : Atom
            The duplicate atom that was created.
        """
        if num is None:
            num = self.slab.n_atoms + 1
        newat = Atom(self.el, self.pos.copy(), num, self.slab)
        newat.cartpos = self.cartpos.copy()
        newat.duplicate_of = self
        newat.site = self.site
        newat.dispInitialized = True
        newat.disp_vib = self.disp_vib
        newat.disp_geo = self.disp_geo
        newat.disp_occ = self.disp_occ

        if add_to_atlists:
            self.slab.atlist.append(newat)                                      # TODO: consider an AtomContainer.add_atom(atom) abstract method!
            self.slab.n_per_elem[self.el] += 1
            if self.layer is not None:
                self.layer.atlist.append(newat)
                newat.layer = self.layer
        return newat

    def initDisp(self, force=False):                                            # TODO: all the DISP stuff should complain if .is_bulk
        """Initialize displacement dictionaries based on self.site.

        This method should not be called before a site has been
        given to this Atom.

        This method must be called before displacements can be
        merged with their offsets, as this method prepares the
        correct entries in disp_occ.

        Parameters
        ----------
        force : bool, optional
            Whether displacements should be cleared also
            if they are already present. Default is False.

        Raises
        ------
        RuntimeError
            If this method is called before a site is available
        """
        if not self.site:
            raise RuntimeError('Cannot initialize displacements '
                               'before a site is defined')
        if self.dispInitialized and not force:
            return

        self.dispInitialized = True
        self.disp_vib = {'all': [0.]}
        self.disp_geo = {'all': [np.zeros(3)]}
        self.disp_occ = {}
        self.disp_center_index = {'vib': {'all': 0},
                                  'geo': {'all': 0},
                                  'occ': {}}
        for k, v in self.site.occ.items():
            if v > 0 or k in self.site.mixedEls:
                self.disp_occ[k] = [v]
                self.disp_center_index['occ'][k] = 0

    def is_same_xy(self, cartpos, eps=1e-3):
        """Return whether this atom is close to a 2D cartpos.

        If the atom is close to an edge or corner, its replicas
        are also considered.

        Parameters
        ----------
        cartpos : numpy.ndarray or Atom
            2D Cartesian coordinates to check against the position
            of this atom. If an Atom, its in-plane Cartesian position
            is used.
        eps : float, optional
            The precision to which positions are expected to match.
            The default is 1e-3.

        Returns
        -------
        bool
            True if positions match, else False.
        """
        if isinstance(cartpos, Atom):
            cartpos = cartpos.cartpos[:2]
        abt = self.slab.ab_cell.T
        releps = eps / np.linalg.norm(abt, axis=1)
        complist, _ = add_edges_and_corners([self.cartpos[:2]], (self.pos,),
                                            releps, abt)
        return any(np.linalg.norm(cartpos - complist, axis=1) < eps)

    def distance(self, cartpos, include_c_replicas=False):
        """Return the distance of this atom from a cartpos.

        2D or 3D replicas of this atom in other unit cells are also 
        considered, and the minimum distance is returned.

        Parameters
        ----------
        cartpos : numpy.ndarray or Atom
            3D Cartesian coordinates to check against the position
            of this atom. If an Atom, its Cartesian position is used.
        include_c_replicas : bool
            Whether replicas in the out-of-plane directions should be 
            considered. The default is False (2D-replicas only).

        Returns
        -------
        float
            The smallest absolute distance among the replicas.
        """
        if isinstance(cartpos, Atom):
            cartpos = cartpos.cartpos
        releps = (1.01, 1.01, 1.01 if include_c_replicas else 0)
        complist, _ = add_edges_and_corners([self.cartpos], (self.pos,),
                                            releps, self.slab.ucell.T)
        return min(np.linalg.norm(cartpos - complist, axis=1))

    def mergeDisp(self, el):
        """
        Merges the offsets from VIBROCC and DISPLACEMENTS into the
        displacements lists from DISPLACEMENTS for the given element.

        For vibration and occupational offsets a consistency check is
        performed. The offset lists will be emptied.

        Raises
        -------
        RuntimeError
            If this method is called before initDisp.
        """
        if not self.dispInitialized:
            raise RuntimeError('Need to .initDisp before displacements '
                               'can be merged with offsets')

        self.storeOriState()
        # geometric offsets from DISPLACEMENTS
        geo_d_offset = self.disp_geo_offset.get(el,
                                          self.disp_geo_offset['all'])[0]
        if el not in self.disp_geo:
            self.disp_geo[el] = copy.copy(list(self.disp_geo['all']))
        self.disp_geo[el] = list(self.disp_geo[el] + geo_d_offset)
        self.disp_geo_offset = {'all': [np.zeros(3)]}

        # geometric offsets from VIBROCC
        if el in self.offset_geo:
            geo_offset = self.offset_geo[el]
            self.disp_geo[el] = [geo_step + geo_offset for geo_step in self.disp_geo[el]]
            del self.offset_geo[el]

        # vibration offsets from VIBROCC
        if el not in self.disp_vib:
            self.disp_vib[el] = copy.copy(self.disp_vib['all'])
        if el in self.offset_vib:
            vib_offset = self.offset_vib[el]
            final_vib_steps = [vib_step + vib_offset for vib_step in self.disp_vib[el]]
            if any(np.array(final_vib_steps) + self.site.vibamp[el] < 0):
                _LOGGER.error(f'Vibration offset for {self} defined in '
                              'VIBROCC would result in negative vibration '
                              'amplitude. Offset will be ignored.')
            else:
                self.disp_vib[el] = final_vib_steps
            del self.offset_vib[el]

        # vibration offsets from VIBROCC
        if el in self.offset_occ:
            occ_offset = self.offset_occ[el]
            final_occ_steps = [occ_step + occ_offset for occ_step in self.disp_occ[el]]
            if any(np.array(final_occ_steps) < 0) or any(np.array(final_occ_steps) > 1):
                _LOGGER.error(
                    f'Occupational offset for {self} defined in '
                    'VIBROCC would result in unphysical concentration'
                    '(occupation <0 or >1). Offset will be ignored.'
                    )
            else:
                self.disp_occ[el] = final_occ_steps
            del self.offset_occ[el]
        self.disp_ranges.update()

    def storeOriState(self):
        """Stores the initial values from the input files for this atom."""
        if self.oriState is None:
            self.oriState = self.duplicate(add_to_atlists=False)                # TODO: potential problem: duplicate only shallow copies some attributes.

    def translate_2d(self, cart_shift, frac_shift):
        """Apply a 2D translation to this Atom."""
        self.cartpos[:2] += cart_shift
        self.pos[:2] += frac_shift


@frozen
class DisplacementRanges:
    """The minimum and maximum displacements for all elements of an atom."""

    _atom: Atom
    geo: Dict[str, tuple] = field(default_factory=dict)
    vib: Dict[str, tuple] = field(default_factory=dict)
    occ: Dict[str, tuple] = field(default_factory=dict)
    __hash__ = None  # Not hashable, as it has mutable attributes

    def __iter__(self):
        """Yield displacement-range dictionaries."""
        yield self.geo
        yield self.vib
        yield self.occ

    def __str__(self):
        """Return a string version of these displacement ranges."""
        lines = [f'Ranges for {self._atom}']
        names = 'geo', 'vib', 'occ'
        for _range, name in zip(self, names):
            lines.extend(f'{name} {el}: {minmax}'
                         for el, minmax in _range.items())
        return '\n'.join(lines)

    def update(self):
        """Update the minimum and maximum displacements."""
        displacements = (self._atom.disp_geo,
                         self._atom.disp_vib,
                         self._atom.disp_occ)
        for disp_type, _range in zip(displacements, self):
            for element, disp in disp_type.items():
                disp = np.asarray(disp)
                new_range = disp.min(axis=0), disp.max(axis=0)
                try:
                    old_range = _range[element]
                except KeyError:
                   _range[element] = new_range
                   continue
                _range[element] = (np.minimum(old_range[0], new_range[0]),
                                   np.maximum(old_range[1], new_range[1]))
