"""Module structdomains of viperleed.gui.leedsim.classes.

Defines the LEEDStructuralDomains class, used for constructing a LEED
pattern with one or more structures, i.e., domains characterized by
superlattice matrices that cannot be reduced to one another by one of
the symmetry operations of the bulk. Those that are related by bulk
symmetry operations are handled by the LEEDSymmetryDomains class.
The LEEDStructuralDomains class is essentially a list of
LEEDSymmetryDomains instances.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-03-21'
__license__ = 'GPLv3+'

from collections.abc import MutableSequence
import copy
import re

import numpy as np

from viperleed.gui.helpers import conventional_angles
from viperleed.gui.leedsim.classes.equivalent_beams import LEEDEquivalentBeams
from viperleed.gui.leedsim.classes.leedparameters import LEEDParametersList
from viperleed.gui.leedsim.classes.symdomains import LEEDSymmetryDomains

from viperleed.gui import decorators as dev_


DEFAULT_NAME = re.compile(r"^S\d+$")


# TODO: We would need to fiddle much less if only we passed
# the LEEDParameters to LEEDEquivalentBeams instead of just
# the superlattices. This would save a lot of potential cache
# clear operations. One would need LEEDParameters to then
# implement an appropriate __hash__. The tricky part there
# is that a == b --> hash(a) == hash(b), so we have to
# care about symmetry-equivalent domains. Also, one should
# have the hash care about the hash(-1) = -2 problem that
# can very likely happen on the superlattices.


class LEEDStructuralDomains(  # pylint: disable=too-many-ancestors
        MutableSequence):
    """Represent several structural domains in a LEED pattern.

    This is a container class that can be used to represent a list of
    distinct structural domains contributing to a LEED pattern. Each
    of the structural domains also has its symmetry-induced domains

    LEEDStructuralDomains(params)

    Parameters
    ----------
    params
        LEEDParametersList or anything that LEEDParametersList accepts
    """

    @dev_.exec_time
    def __init__(self, params, keep_duplicates=False):
        """Initialize LEEDStructuralDomains instance.

        Parameters
        ----------
        params : iterable or LEEDParametersList
            Each element is either dict, ConfigParser, or
            LEEDParameters. Will be fed to a LEEDParametersList.
        keep_duplicates : bool, default=False
            If True, parameters giving identical LEED
            patterns are kept. Otherwise only those
            giving 'unique' patterns are kept.

        Returns
        -------
        None.
        """
        self.__parameters = LEEDParametersList(params, keep_duplicates)
        self.__list = [LEEDSymmetryDomains(p) for p in self.parameters]

        # __last_eq is a reference to the last LEEDEquivalentBeams
        # instance used to combine all structural and symmetry domains.
        # Accessed (read-only) via .eq_beams_last_config
        # __last_domains a dictionary of the last domain indices used.
        # Accessed (read-only) via .last_domains_used
        # Both are set in self.equivalent_spots()
        self.__last_eq = None
        self.__last_domains = {}

    def __delitem__(self, elem):
        """Delete item at a given index or slice."""
        ids_deleted = self.domain_ids[elem]

        if isinstance(elem, slice):
            to_delete = range(*elem.indices(self.n_domains))
        else:
            to_delete = [elem]

        # Indices of domains without a given name, before
        # deletion, excluding those that will be deleted:
        old_unnamed = [i
                       for i, dom_id in enumerate(self.domain_ids)
                       if DEFAULT_NAME.mathch(dom_id) and i not in to_delete]

        # Actually delete stuff
        del self.__list[elem]
        del self.__parameters[elem]

        if self.eq_beams_last_config is None:
            return

        # Now take care of the LEEDEquivalentBeams:
        new_unnamed = [i
                       for i, dom_id in enumerate(self.domain_ids)
                       if DEFAULT_NAME.mathch(dom_id)]

        # If any of the unnamed has changed its 'place', we
        # have to clear the cache of LEEDEquivalentBeams, as
        # their index-based names will now be different.
        # TODO: remove from the cache only those that may have
        # gotten messed up
        if any(new_idx != old_idx
               for new_idx, old_idx in zip(new_unnamed, old_unnamed)):
            self.eq_beams_last_config.clear_cache('self')

        # If deletion affected the domains used
        # in the last .equivalent_spots() call,
        # that's no longer up to date.
        if any(old_id in self.last_domains_used for old_id in ids_deleted):
            self.__last_eq = None
            self.__last_domains = {}

    def __getitem__(self, elem):
        """Retrieve item at a given index or slice.

        Notice that passing a slice will return a list object, not a
        LEEDStructuralDomains object!
        """
        return self.__list[elem]

    def __setitem__(self, idx, data):
        """Set the value of element(s) at a given index or slice.

        The input is checked for consistency, the same way as in
        LEEDParametersList.__setitem__(). Replacing domains with
        others having the same name and superlattice, but having
        e.g., a different plane group, will trigger a clearing
        of the LEEDEquivalentBeams cache. This means that some
        of the old configurations will need to be recalculated.

        Parameters
        ----------
        idx : int or slice
        data : anything that can instantiate a LEEDParametersList
        """
        # Keep a few values that we will use later on
        old_parameters = copy.deepcopy(self.parameters)
        old_ids = self.domain_ids
        old_superlattices = self.superlattices

        # Leave all the processing, as well as consistency and
        # duplicates check to LEEDParametersList.__setitem__
        self.parameters[idx] = data

        # Then keep track of what was truly changed
        replaced_idx = [i
                        for i, (old, new) in enumerate(zip(old_parameters,
                                                           self.parameters))
                        if new != old]
        for i in replaced_idx:
            self.__list[i] = LEEDSymmetryDomains(self.parameters[i])

        # Now decide whether we have to fix the LEEDEquivalentBeams:
        # - if any of those removed has the same name and superlattice
        #   as one of those inserted --> hash corrupted. clear cache
        # - if any of those removed was used last time
        #   --> invalidate __last_eq
        if not self.eq_beams_last_config:
            return

        # Select only the replaced old/new ids and superlattices
        new_ids = self.domain_ids
        new_superlattices = self.superlattices
        old_ids, new_ids = zip(*((old_ids[i], new_ids[i])
                                 for i in replaced_idx))
        old_superlattices, new_superlattices = zip(*(
            (old_superlattices[i], new_superlattices[i])
            for i in replaced_idx
            ))

        for old_id, old_superlattice in zip(old_ids, old_superlattices):
            same_name_and_super = [np.all(sup == old_superlattice)
                                   for sup, new_id in zip(new_superlattices,
                                                          new_ids)
                                   if new_id == old_id]
            if any(same_name_and_super):
                self.eq_beams_last_config.clear_cache('all')
                break
        if any(old_id in self.last_domains_used for old_id in old_ids):
            self.__last_eq = None
            self.__last_domains = {}

    def __len__(self):
        """Return number of structural domains."""
        return len(self.__list)

    def __repr__(self):
        """Return string representation of self."""
        txt = 'LEEDStructuralDomains([\n'
        txt += '\n'.join(repr(d) for d in self.parameters)
        txt += '\n])'
        return txt

    @property
    def bulk_basis(self):
        """Return the reciprocal-space basis of the bulk lattice."""
        return self[0].bulk_basis

    @property
    def denominators(self):
        """Return a dictionary of the denominators for beam indices.

        Returns
        -------
        dict
            key : str
                Same values as in self.domain_ids
            value : int
                Denominator to be used for the beams coming
                from the structural domain at key
        """
        return {struct_id: struct.denominator_for_bulk_beams
                for struct_id, struct in zip(self.domain_ids, self)}

    @property
    def domain_ids(self):
        """Return an up-to-date list of domain 'ids'.

        Each id is either a name, if a name is present in the
        corresponding LEEDParameters dictionary, or a string
        "S"+number, with number-1 is the index in self

        Returns
        -------
        list of str
        """
        return [param['name'] if param['name'] else f"S{i+1}"
                for i, param in enumerate(self.parameters)]

    @property
    def eq_beams_last_config(self):
        """Return the LEEDEquivalentBeams of the last configuration."""
        return self.__last_eq

    @property
    def extinct(self):
        """Return beams that are extinct in the last configuration.

        If the configuration (i.e., domains selected, beam-incidence
        angles) is not up to date, call .equivalent_spots() before!

        Returns
        -------
        list
            List of extinct beams as either tuples or BeamIndex,
            depending on whether self.fractional is False or True

        Raises
        ------
        RuntimeError
            When this property is called before the .equivalent_spots()
            method has been invoked at least once.
        """
        if self.eq_beams_last_config is None:
            raise RuntimeError("LEEDStructuralDomains: Need to "
                               "call .equivalent_spots() once "
                               "before .extinct is available")
        return self.eq_beams_last_config.extinct

    @property
    def fractional(self):
        """Return whether beam indices should be expressed fractional.

        Return whether the beam indices need to be expressed as
        fractional, or if it is enough to use the numerators. In the
        latter case the processing of equivalent beams is faster.

        Return
        ------
        bool
            True if it requires fractional indices.
        """
        denominators = [s.denominator_for_bulk_beams for s in self]
        return any(den != denominators[0] for den in denominators[1:])

    @property
    def indexed_beams(self):
        """Return a list of indexed beams.

        If the configuration (i.e., domains selected, beam-incidence
        angles) is not up to date, call .equivalent_spots() before!

        Returns
        -------
        dict
            keys : {'beams', 'group_indices',
                    'overlap_domains', 'extinct_domains'}
            values : list
                One element per each one of the beams generated
                from the parameters passed at instantiation.

        Each element in the values is:
            'beams' : tuple or BeamIndex
                Tuple with only numerators if not self.fractional.
                Full fractional index of type BeamIndex otherwise.
            'group_indices' : int
                Progressive index that uniquely identifies the beam
                equivalence class to which beam belongs.  If the beam
                is fully extinct, its group index is negative.
            'overlap_domains' : list
                List of domain identifiers specifying which domains
                contribute to this beam. Each element is a tuple
                of the form (struct_id, (symm_id1, ...)), where
                struct_id is one of self.domain_ids, and symm_id
                one of LEEDSymmetryDomains.domain_ids
            'extinct_domains' : list
                List of domains, among those that overlap, that
                only contribute to the beam in an extinct fashion.
                Each entry has the same form as in 'overlap_domains'.
                If a structure contributes only with non-extinct
                contributions (i.e., all its symmetry domains are
                non-extinct), it will not appear. If a subset of
                the symmetry domains of one structure contributes
                in an extinct fashion, only these will be listed
                in the second element of (struct_id, (sym_id1, ...)).

        Raises
        ------
        RuntimeError
            When this property is called before the equivalent_spots()
            method has been invoked at least once.
        """
        if self.eq_beams_last_config is None:
            raise RuntimeError("LEEDStructuralDomains: Need to call "
                               ".equivalent_spots() once before "
                               ".indexed_beams is available")
        return self.eq_beams_last_config.indexed_beams

    @property
    def last_domains_used(self):
        """Return identifiers for the domains in the last configuration."""
        if self.eq_beams_last_config is None:
            raise RuntimeError("Need to call .equivalent_spots() once "
                               "before .last_domains_used is available")
        return self.__last_domains

    @property
    def n_domains(self):
        """Return number of distinct domains produced by bulk symmetry."""
        return len(self)

    @property
    def parameters(self):
        """Return the LEEDParametersList that underlies self."""
        return self.__parameters

    @property
    def superlattices(self):
        """Superlattice matrices of each domain.

        Returns
        -------
        dict
            Dictionary of the form {id: superlattices} for each of the
            structural domains, where superlattices is a numpy.ndarray
            of superlattice matrices
        """
        superlattices = (p.superlattices for p in self)
        return dict(zip(self.domain_ids, superlattices))

    def insert(self, index, value,  # pylint: disable=arguments-differ
               keep_duplicates=False):
        """Override list.insert allowing check for duplicates."""
        # the next raises errors if data is not compatible
        old_ids = self.domain_ids
        old_superlattices = self.superlattices
        old_len = len(self.parameters)
        self.parameters.insert(index, value, keep_duplicates)

        # do something only if something was actually inserted
        if len(self.parameters) == old_len:
            return
        self.__list.insert(index, value)

        # See if we should fiddle with LEEDEquivalentBeams
        if not self.eq_beams_last_config:
            return

        # If we shifted some unnamed structures we will have to
        # clear the 'self'. If the inserted structure
        # has the same name and superlattice of any of the
        # ones already present we will also need to clear
        # the 'dict' cache.
        if any(DEFAULT_NAME.match(dom_id) for dom_id in old_ids[index:]):
            self.eq_beams_last_config.clear_cache('self')
        new_id = self.domain_ids[index]
        new_superlattice = self.superlattices[index]
        for old_id, superlattice in zip(old_ids, old_superlattices):
            if new_id == old_id and np.all(superlattice == new_superlattice):
                self.eq_beams_last_config.clear_cache('dict')
                break

    def append(self, value,  # pylint: disable=arguments-differ
               keep_duplicates=False):
        """Override list.append allowing check for duplicates."""
        self.insert(len(self), value, keep_duplicates)

    def extend(self, values,  # pylint: disable=arguments-differ
               keep_duplicates=False):
        """Override list.extend allowing check for duplicates."""
        for value in values:
            self.append(value, keep_duplicates)

    @dev_.exec_time
    # @dev_.profile_lines
    def equivalent_spots(self, domains=None, theta=None, phi=None):
        """Return spots of selected domains that are equivalent.

        Return a list of beams grouped by equivalence, given a primary
        beam direction defined by polar and azimuthal angles, for the
        selected domains

        Parameters
        ----------
        domains : dict, list, int or None (default=None)
            When passing a dictionary, the keys are used to select
            the structural domains based on their name; names not
            found are skipped. When passing a list it should have as
            many elements as there are domains. In both cases each
            value/element can be either None (selects all domains),
            an integer, or a list of integers, selecting which
            symmetry-equivalent domains to use. When a single
            integer is passed, this is used as the positional index
            for all the domains.
        theta, phi : number or None
            polar and azimuthal directions in degrees of the primary
            beam.  If not given or None, those defined at instantiation
            are used.  theta should be in the [-90, 90] range.
            phi is positive counterclockwise, and measured from the x
            axis in the Cartesian reference of the real-space basis of
            the first domain.

        Returns
        -------
        equivalent_spots : list
            list of beams grouped by equivalence
        """
        # Figure out which domains we want:
        domains, structures, struct_ids = self.process_domains_input(domains)

        # Check and process the angles
        if theta is None:
            theta = self.parameters[0]['beamIncidence'][0]
        if phi is None:
            phi = self.parameters[0]['beamIncidence'][1]
        theta, phi = conventional_angles(float(theta), float(phi))

        # Call the equivalent_spots() method of each of the relevant
        # LEEDSymmetryDomains instances so that the LEEDEquivalentBeams
        # instances are correctly initialized
        for symmetry_domains, structure in zip(domains, structures):
            structure.equivalent_spots(symmetry_domains, theta, phi)

        # Now each structure.eq_beams_last_config contains a
        # reference to the correct LEEDEquivalentBeams instance.
        # Combine them into another LEEDEquivalentBeams instance
        # that includes all structures.  Also build ids.
        eq_beams = [s.eq_beams_last_config for s in structures]
        symm_ids = [s.last_domains_used for s in structures]
        combined_ids = tuple(zip(struct_ids, symm_ids))
        self.__last_domains = dict(combined_ids)

        # Check if we need to pass denominators together with
        # the numerators.  This is necessary if we have multiple
        # structures with differently sized cells.
        denominators = []
        if self.fractional:
            denominators = [self.denominators[struct_id]
                            for struct_id in struct_ids]

        self.__last_eq = LEEDEquivalentBeams(eq_beams,
                                             domain_ids=combined_ids,
                                             denominators=denominators,
                                             caller=self)
        if not self.eq_beams_last_config.was_cached:
            self.__correct_beam_index()
        elif len(combined_ids) == 1:
            # Single domain. indexed_beams['overlap_domains']
            # and indexed_beams['extinct_domains'] need fixing
            indexed = self.indexed_beams
            struct_id = combined_ids[0][0]
            indexed['overlap_domains'] = [
                [(struct_id, tuple(sym_ids))]
                for sym_ids in indexed['overlap_domains']
                ]
            indexed['extinct_domains'] = [
                [(struct_id, tuple(sym_ids))] if sym_ids else []
                for sym_ids in indexed['extinct_domains']
                ]
        return self.indexed_beams

    def from_domain_id(self, domain_id):
        """Return the LEEDSymmetryDomains with a certain id.

        Parameters
        ----------
        domain_id : str
            Identifier of the domain. If using the
            positional index, use self[idx] instead.

        Returns
        -------
        LEEDSymmetryDomains
            The LEEDSymmetryDomains instance with
            the given domain_id

        Raises
        ------
        ValueError
            If there is no domain with id == domain_id
        """
        try:
            idx = self.domain_ids.index(domain_id)
        except ValueError as err:
            raise ValueError(f"{domain_id!r} not found.") from err
        return self[idx]

    def set_beam_incidence(self, theta=None, phi=None):
        """Edit the beam incidence.

        Also triggers an update of the beam equivalence
        if self.equivalent_spots() was already called.

        Parameters
        ----------
        theta : float, optional
            Polar incidence angle in degrees, measured with
            respect to the direction perpendicular to the surface.
            Range [-90, 90].
        phi : float, optional
            Azimuthal angle of incidence in degrees, measured with
            respect to the positive x axis in the Cartesian reference
            of the bulk.  Positive counterclockwise when looking down
            on the surface.

        Returns
        -------
        None.
        """
        if theta is None and phi is None:
            return

        old_angles = self.parameters[0]['beamIncidence']

        # Update the underlying parameters, check if anything changed
        self.parameters.set_beam_incidence(theta, phi)

        new_angles = self.parameters[0]['beamIncidence']
        if new_angles == old_angles:
            return

        # Update the individual symmetry domains. This
        # triggers an update of the LEEDSymmetryDomains if
        # their .equivalent_spots() method was called before.
        for sym_domain in self:
            sym_domain.set_beam_incidence(*new_angles)
        # Finally update the full spot equivalence
        if self.eq_beams_last_config is not None:
            self.equivalent_spots(self.last_domains_used)

    def process_domains_input(self, domains):
        """Get correct structural domains. Used in equivalent_spots.

        Parameters
        ----------
        domains : dict, list, int or None (default=None)
            When passing a dictionary, the keys are used to select
            the structural domains based on their name; names not
            found are skipped. When passing a list it should have as
            many elements as there are domains. In both cases each
            value/element can be either None (selects all domains),
            an integer, or a list of integers, selecting which
            symmetry-equivalent domains to use. When a single
            integer is passed, it is used for all the domains.

        Returns
        -------
        domains : list
            Fixed list of domain indices
        structures : list
            Which structures were requested. Each element is a
            LEEDSymmetryDomains instance
        ids : list
            Which of the ids correspond to the structures

        Raises
        ------
        ValueError
            When the number of domain indices passed is inconsistent
            (i.e., smaller or larger) than the number of structural
            domains.
        TypeError
            When the entries in the domain indices are not integers
        """
        ids = self.domain_ids
        structures = self.__list
        if domains is None:
            domains = [None]*self.n_domains
        elif isinstance(domains, dict):
            domains_dict = domains
            ids = []
            domains = []
            structures = []
            for struct_id, symm_ids in domains_dict.items():
                try:
                    structure = self.from_domain_id(struct_id)
                except ValueError:
                    pass
                else:
                    ids.append(struct_id)
                    domains.append(symm_ids)
                    structures.append(structure)
        elif hasattr(domains, '__len__') and len(domains) != self.n_domains:
            raise ValueError("equivalent_spots: Inconsistent number "
                             "of structural domains. Expected "
                             f"{self.n_domains}, found {len(domains)}")
        elif isinstance(domains, int):
            domains = [domains]*self.n_domains
        else:
            raise TypeError("equivalent_spots: Invalid domains "
                            f"{type(domains).__name__}")
        return domains, structures, ids

    @dev_.exec_time
    def __correct_beam_index(self):
        """Correct symmetry domain indices in self.indexed_beams.

        When the LEEDEquivalentBeams combines different structural
        domains, it places in the 'overlapping list' and 'extinct list'
        the full domain ids used, e.g., [('S1', (0, 1, 2)), ...], even
        if the beam being considered only has contributions from domain
        1 of structure 'S1', as this information is not available.

        Here we sort out this mess by taking the LEEDSymmetryDomains
        that were used for the superposition, and selecting only
        those that contain each beam.

        This method is used only internally, and is called
        at the end of self.equivalent_spots() on an incorrect
        beam list.  One can safely assume that a beamlist is
        correct if the LEEDEquivalentBeams was loaded from cache
        """
        # NB: one could do this processing already as part of
        # LEEDEquivalentBeams__.index_beams(), but it does not
        # really give a significant improvement (1.8s tot vs.
        # 1.9s tot for a LEEDPattern with fishbone + 1x1).
        # This does not justifz the mess of extra code that is
        # needed there, which also makes LEEDEquivalentBeams
        # partly dependent on LEEDStructuralDomains.
        if self.eq_beams_last_config is None:
            raise RuntimeError("Cannot reindex beams if .equivalent_spots() "
                               "was never called before.")
        indexed = self.indexed_beams

        # for i, (beam, group, overlap, extinct) in enumerate(indexed):
        for i, (beam, _,
                overlap, extinct) in enumerate(zip(*indexed.values())):
            correct_overlap = []
            correct_extinct = []
            extinct = dict(extinct)
            for struct_id, _ in overlap:
                dom_beam = beam
                if self.fractional:
                    den = self.denominators[struct_id]
                    dom_beam = (round(beam[0]*den), round(beam[1]*den))
                struct = self.from_domain_id(struct_id)
                corrected = struct.domains_containing_beam(dom_beam)
                correct_overlap.append((struct_id, corrected))

                extinct_ids = extinct.get(struct_id, None)
                if extinct_ids:
                    correct_extinct.append(
                        (struct_id,
                         tuple(symm_id
                               for symm_id in extinct_ids
                               if symm_id in corrected))
                        )

            # indexed[i] = (beam, group, correct_overlap, correct_extinct)
            indexed['overlap_domains'][i] = correct_overlap
            indexed['extinct_domains'][i] = correct_extinct
